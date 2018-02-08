"""Classes for parsing and analyzis of DNA migration data."""
from collections import OrderedDict, Counter
import matplotlib.pyplot as plt
from io import BytesIO
from matplotlib.backends.backend_pdf import PdfPages

from bandwagon import BandsPattern, BandsPatternsSet
import flametree

try:
    from plateo.exporters.plate_to_matplotlib_plots import PlateColorsPlotter
    from plateo.parsers import plate_from_platemap_spreadsheet
    from plateo.containers import Plate96
    PLATEO_AVAILABLE = True
except:
    PLATEO_AVAILABLE = False


from ..bands_predictions import predict_digestion_bands
from ..plots import plot_cuts_map, plot_all_constructs_cuts_maps
from ..tools import load_genbank
from .Clone import Clone
from .BandsObservation import BandsObservation


class ClonesObservations:
    """All useful informations for a collection of clones to be validated.

    Parameters
    ----------
    clones
      Either a list of Clones (each with a unique name) or a dictionnary
      ``{name: Clone}``.

    constructs_records
      A dictionnary ``{construct_name: biopython_record}`` indicating the
      sequence of the different constructs. For each construct, set the
      attribute ``construct.linear = False`` if the construct is circular.

    """

    def __init__(self, clones, constructs_records):
        """Initialize"""
        if isinstance(clones, (list, tuple)):
            clones = {clone.name: clone for clone in clones}
        self.clones = clones
        self.constructs_records = constructs_records
        self.constructs_digestions = {ct: {} for ct in constructs_records}


    def get_digestion_bands_for_construct(self, construct_id, digestion):
        """Return the bands resulting from the digestion of the construct.

        This function enables some memoization (no digestion is computed
        twice).

        Parameters
        ----------

        construct_id
          ID of the construct as it appears in this object'
          ``constructs_records``.

        digestion
          For instance ``('BamHI', 'XbaI')``

        """
        construct_digestions = self.constructs_digestions[construct_id]
        if digestion not in construct_digestions:
            construct_record = self.constructs_records[construct_id]
            construct_digestions[digestion] = predict_digestion_bands(
                str(construct_record.seq),
                linear=construct_record.__dict__.get('linear', True),
                enzymes=digestion
            )
        return construct_digestions[digestion]

    def get_clone_digestion_bands(self, clone, construct_id=None):
        """Return ``{digestion: bands}`` for all digestions of this clone."""
        if construct_id is None:
            construct_id = clone.construct_id
        return {
            digestion: self.get_digestion_bands_for_construct(
                construct_id=clone.construct_id, digestion=digestion)
            for digestion in clone.digestions
        }

    def validate_all_clones(self, relative_tolerance=0.1,
                            min_band_cutoff=None, max_band_cutoff=None):
        """Return ``{clone: CloneValidation}`` for all clones."""
        return OrderedDict([
            (clone_name, clone.validate_bands(
                bands_by_digestion=self.get_clone_digestion_bands(clone),
                relative_tolerance=relative_tolerance,
                min_band_cutoff=min_band_cutoff,
                max_band_cutoff=max_band_cutoff
            ))
            for clone_name, clone in self.clones.items()
        ])

    def identify_all_clones(self, relative_tolerance=0.1,
                            min_band_cutoff=None, max_band_cutoff=None):
        """Return ``{clone: {construct_id: CloneValidation}}`` for all clones.

        Parameters
        ----------
        relative_tolerance
          Tolerance, as a ratio of the full ladder span. If =0.1, then the
          discrepancy will have a value of 1 when a band's nearest
          correspondent in the other pattern is more that 10% of the ladder
          span apart.

        min_band_cutoff
          Discrepancies involving at least one band below this minimal band
          size will be ignored. By default, it will be set to the smallest
          band size in the ladder.

        max_band_cutoff
          Discrepancies involving at least one band above this minimal band
          size will be ignored. By default, it will be set to the smallest
          band size in the ladder.

        """
        return OrderedDict([
            (clone_name, {
                construct_id: clone.validate_bands(
                    bands_by_digestion=self.get_clone_digestion_bands(
                        clone, construct_id=construct_id),
                    relative_tolerance=relative_tolerance,
                    min_band_cutoff=min_band_cutoff,
                    max_band_cutoff=max_band_cutoff
                )
                for construct_id in self.constructs_records
            })
            for clone_name, clone in self.clones.items()
        ])

    def validations_summary(self, validations, sort_clones_by_score=True):
        """Return ``{construct_id: [CloneValidation, ...]}``.

        To each construct corresponds a list of the validation of all clones
        associated with that construct, from the best-scoring to the
        least-scoring.
        """
        results = OrderedDict()
        for clone_name, validation in validations.items():
            construct_id = self.clones[clone_name].construct_id
            if construct_id not in results:
                results[construct_id] = []
            results[construct_id].append(validation)
        if sort_clones_by_score:
            for cst_id, validations in results.items():
                results[cst_id] = sorted(validations,
                                         key=lambda v: v.max_discrepancy)
        return results

    def plot_validations_plate_map(self, validations, target=None,
                                   ax=None):
        """Plot a map of the plate with passing/failing wells in green/red."""

        if not PLATEO_AVAILABLE:
            raise ImportError('This function requires Plateo installed')
        plate = Plate96('Validations map')
        constructs_colors = {}
        for well in plate.iter_wells():
            if well.name not in validations:
                well.data.bg_color = None
                well.data.pass_color = (1, 1, 1, 0.1)
            else:
                validation = validations[well.name]
                construct = validation.clone.construct_id
                if construct not in constructs_colors:
                    alpha = ((0.35 * len(constructs_colors) % 1) + .2) / 2
                    constructs_colors[construct] = (.4, .4, .76, alpha)
                well.data.bg_color = constructs_colors[construct]
                well.data.pass_color = ((.43, .92, .49) if validation.passes
                                        else (.91, 0.43, 0.43))
        bg_plotter = PlateColorsPlotter(lambda w: w.data.bg_color)
        ax, _ = bg_plotter.plot_plate(plate, figsize=(10, 5), ax=ax)
        pass_plotter = PlateColorsPlotter(lambda w: w.data.pass_color,
                                          well_radius=250)
        pass_plotter.plot_plate(plate, ax=ax)
        if target is not None:
            ax.figure.savefig(target, bbox_inches='tight')
            plt.close(ax.figure)
        else:
            return ax

    def plot_all_validations_patterns(self, validations, target=None,
                                      per_digestion_discrepancy=False):
        """Plot a Graphic report of the gel validation.

        The report displays the patterns, with green and red backgrounds
        depending on whether they passed the validation test.

        Parameters
        ----------

        target_file
          File object or file path where to save the figure. If provided, the
          function returns none, else the function returns the axes array of
          the final figure

        relative_tolerance
          Relative error tolerated on each band for the ovserved patterns to
          be considered similar to the expected patterns.

        min_band_cutoff
          Bands with a size below this value will not be considered

        max_band_cutoff
          Bands with a size above this value will not be considered

        """

        summary = self.validations_summary(validations)
        max_x = Counter([
            clone.construct_id
            for clone in self.clones.values()
        ]).most_common(1)[0][1] + 2
        digestions_by_construct = {construct: [] for construct in summary}
        for construct, validations in summary.items():
            for validation in validations:
                for digestion in validation.discrepancies:
                    if digestion not in digestions_by_construct[construct]:
                        digestions_by_construct[construct].append(digestion)
        ladder = list(validations[0].clone.digestions.values())[0].ladder
        pdf_io = BytesIO()
        with PdfPages(pdf_io) as pdf:
            for construct_id, validations in summary.items():
                digestions = digestions_by_construct[construct_id]
                band_patterns = OrderedDict()
                for i, digestion in enumerate(digestions):
                    reference_bands = self.get_digestion_bands_for_construct(
                        construct_id, digestion)
                    reference = BandsPattern(
                        bands=reference_bands,
                        ladder=ladder,
                        label="exp." if (i == 0) else None,
                        background_color="#c6dcff",
                        corner_note="Total: %d bp" % sum(reference_bands),
                        global_bands_props={"label_fontdict": {"size": 5}}
                    )
                    sorted_bands = sorted(reference.bands,
                                          key=lambda b: -b.dna_size)
                    for band_name, band in zip("abcdefghijklm", sorted_bands):
                        band.label = band_name
                    patterns = [
                        validation.to_bandwagon_bandpattern(
                            digestion=digestion,
                            per_digestion_discrepancy=per_digestion_discrepancy,
                            label='auto' if (i == 0) else None
                        )
                        for validation in validations
                    ]
                    band_patterns[digestion] = BandsPatternsSet(
                        [reference] + patterns, ladder=ladder,
                        label=" + ".join(digestion),
                        ladder_ticks=5,
                        global_patterns_props={
                            "label_fontdict": {"rotation": 60}}
                    )

                digestions = digestions_by_construct[construct_id]
                fig, axes = plt.subplots(
                    len(digestions), 1, sharex=True,
                    figsize=(1.1 * max_x, 3 * len(digestions)),
                )
                if len(digestions) == 1:
                    axes = [axes]
                axes[-1].set_xlim(xmax=max_x)

                for ax, (dig, pattern_set) in zip(axes, band_patterns.items()):
                    pattern_set.plot(ax)

                axes[-1].set_xlabel(construct_id,
                                    fontdict=dict(size=14, weight='bold'),
                                    ha='left', position=(0.0, 0.1))
                pdf.attach_note(construct_id)
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)
        pdf_data = pdf_io.getvalue()
        if target is None:
            return pdf_data
        elif target == 'base64':
            return 'data:application/pdf;base64,' + pdf_data.decode("utf-8")
        else:
            with open(target, 'wb') as f:
                f.write(pdf_data)

    def write_identification_report(self, target_file=None,
                                    relative_tolerance=0.1,
                                    min_band_cutoff=None,
                                    max_band_cutoff=None,
                                    plot_constructs_cuts=False):
        """Plot a Graphic report of the gel validation.

        The report displays the patterns, with green and red backgrounds
        depending on whether they passed the validation test.

        Parameters
        ----------

        target_file
          File object or file path where to save the figure. If provided, the
          function returns none, else the function returns the axes array of
          the final figure

        relative_tolerance
          Relative error tolerated on each band for the ovserved patterns to
          be considered similar to the expected patterns.

        min_band_cutoff
          Bands with a size below this value will not be considered

        max_band_cutoff
          Bands with a size above this value will not be considered

        """
        bands_validities = self.compute_all_bands_validities(
            relative_tolerance=relative_tolerance,
            min_band_cutoff=min_band_cutoff,
            max_band_cutoff=max_band_cutoff)
        L = len(self.constructs_records)
        max_x = max(
            len(measures)
            for measures in self.observed_bands.values()
        ) + 1
        fig, axes = plt.subplots(L, 2, figsize=(2.2 * max_x, 3 * L))
        axes_validities = zip(axes, bands_validities.items())
        for (ax1, ax2), (construct_id, validities) in axes_validities:
            reference = BandsPattern(
                self.expected_bands[construct_id], ladder=self.ladder,
                label="exp.", background_color="#c6dcff",
                corner_note="Total: %d bp" % sum(
                    self.expected_bands[construct_id]),
                global_bands_props={"label_fontdict": {"size": 5}}
            )
            sorted_bands = sorted(reference.bands, key=lambda b: -b.dna_size)
            for band_name, band in zip("abcdefghijklm", sorted_bands):
                band.label = band_name
            patterns = [

                BandsPattern(
                    self.observed_bands[construct_id][measure_name],
                    corner_note="Total: %d bp" % sum(
                        self.observed_bands[construct_id][measure_name]),
                    ladder=self.ladder,
                    label=measure_name,
                    gel_image=self.migration_images[measure_name],
                    background_color="#aaffaa" if
                                     validities[measure_name]
                                     else "#ffaaaa"
                )
                for measure_name in self.observed_bands[construct_id]
            ]
            patterns_set = BandsPatternsSet(
                [reference] + patterns, ladder=self.ladder,
                label=construct_id, ladder_ticks=5,
                global_patterns_props={"label_fontdict": {"rotation": 60}}
            )
            ax1.set_xlim(0.5, max_x + 2)
            record = self.constructs_records[construct_id]
            plot_cuts_map(record, enzymes=self.digestions[construct_id],
                          ax=ax2)
            patterns_set.plot(ax1)
        fig.subplots_adjust(hspace=0.3)
        if target_file is not None:
            fig.savefig(target_file, bbox_inches="tight")
            plt.close(fig)
        else:
            return axes

    def plot_all_constructs_cuts_maps(self, target=None, figsize=(12, 4)):
        """Plot schemas of all constructs with cuts, in a multipage PDF.

        Parameters
        ----------
        target
          Either None (at which case the function returns raw data of the PDF)
          or 'base64' (the function returns base64-encoded data, e.g. for
          web transfer), or a file path to be written to
        figsize
          The size in inches of each figure (=page of the pdf).

        """
        constructs_digestions = OrderedDict([
            (construct_id, set())
            for construct_id in self.constructs_records
        ])
        for clone in self.clones.values():
            for digestion in clone.digestions:
                constructs_digestions[clone.construct_id].add(digestion)
        return plot_all_constructs_cuts_maps([
            (self.constructs_records[construct_id], digestion)
            for construct_id, digestions in constructs_digestions.items()
            for digestion in sorted(digestions)
        ], target=target, figsize=figsize)

    def from_files(records_path, constructs_map_path, aati_zip_path,
                   digestions_map_path=None, digestion=None):
        constructs_records = {
            f._name_no_extension: load_genbank(
                f.open('r'), linear=False, name=f._name_no_extension)
            for f in flametree.file_tree(records_path)._all_files
        }
        constructs_records = OrderedDict(sorted(constructs_records.items()))

        constructs_plate = plate_from_platemap_spreadsheet(
            constructs_map_path, data_field='construct', headers=True)
        constructs_map = OrderedDict([
            (well.name, well.data.construct)
            for well in constructs_plate.iter_wells(direction='row')
            if 'construct' in well.data
            and str(well.data.construct) != 'nan'
        ])

        if digestion:
            digestion = tuple(digestion)
            digestions_map = OrderedDict([
                (wellname, digestion)
                for wellname in constructs_map
            ])
        else:
            digestions_plate = plate_from_platemap_spreadsheet(
                digestions_map_path, data_field='digestion', headers=True)
            digestions_map = OrderedDict([
                (well.name, tuple(well.data.digestion.split(', ')))
                for well in digestions_plate.iter_wells(direction='row')
                if 'digestion' in well.data
                and str(well.data.digestion) != 'nan'
            ])

        observations = BandsObservation.from_aati_fa_archive(aati_zip_path)
        clones = Clone.from_bands_observations(observations, constructs_map,
                                               digestions_map)
        clones_observations = ClonesObservations(clones, constructs_records)
        return clones_observations
