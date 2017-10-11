from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as plt
from Bio import Restriction

from dna_features_viewer import (GraphicRecord, GraphicFeature,
                                 BiopythonTranslator)
from bandwagon import BandsPattern, BandsPatternsSet
from bandwagon.ladders import ladder_from_aati_fa_calibration_table

try:
    from plateo.parsers import plate_from_aati_fragment_analyzer_zip
    from plateo.exporters.plate_to_matplotlib_plots import PlateColorsPlotter
    from plateo.containers import Plate96
    PLATEO_AVAILABLE = True
except:
    PLATEO_AVAILABLE = False

from .tools import band_patterns_discrepancy
from .bands_predictions import predict_digestion_bands


class BandsObservation:
    """One observation of a bands pattern.

    Parameters
    ----------
    name
      Name of the observation (used for the top label in plots)

    bands
      A BandPattern object

    ladder
      A BandPattern object representing a ladder.

    migration_image
      Optional RGB array (HxWx3) representing the gel "image" (which will be
      displayed on the side of the plot).

    """

    def __init__(self, name, bands, ladder, migration_image=None):
        """Initialize."""
        self.name = name
        self.bands = sorted(bands)
        self.ladder = ladder
        self.migration_image = migration_image

    @staticmethod
    def from_aati_fa_archive(archive_path, min_rfu_size_ratio=0.3):
        """Return a dictionnary of all band observations in AATI output files.

        Parameters
        ----------
        archive_path
          A path to a ZIP file containing all files output by the AATI fragment
          analyzer.

        min_rfu_size_ratio
          Cut-off ratio to filter-out bands whose intensity is below some
          threshold. The higher the value, the more bands will be filtered out.

        Returns
        -------
        A dictionary ``{'A1': BandsObservation(), 'A2': ...}`` containing the
        measured pattern information for a whole 96-well microplate.

        """
        if not PLATEO_AVAILABLE:
            raise ImportError("Plateo must be installed to parse AATI zips.")

        plate = plate_from_aati_fragment_analyzer_zip(archive_path)
        ladder_data = plate.data.ladder
        ladder = ladder_from_aati_fa_calibration_table(dataframe=ladder_data)
        # get the set of constructs sorted by order of columnwise appearance:

        def band_is_strong_enough(band):
            """Return True iff the band's intensity is above the set level."""
            return (1.0 * band["RFU"] / band["Size (bp)"] > min_rfu_size_ratio)

        return {
            well.name: BandsObservation(
                name=well.name,
                ladder=ladder,
                bands=[
                    band["Size (bp)"]
                    for band in well.data.bands.values()
                    if band_is_strong_enough(band)
                ],
                migration_image=well.data.migration_image
            )
            for well in plate.iter_wells()
        }

    def patterns_discrepancy(self, other_bands, relative_tolerance=0.1,
                             min_band_cutoff=None, max_band_cutoff=None):
        """Return the maximal discrepancy between two band patterns.

        The discrepancy is defined as the largest distance between a band in
        one pattern and the closest band in the other pattern.

        Parameters
        ----------
        other_bands
          A list of bands (integers) to be compared with the current bands

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
        ladder_min, ladder_max = self.ladder.dna_size_span
        if min_band_cutoff is None:
            min_band_cutoff = ladder_min
        if max_band_cutoff is None:
            max_band_cutoff = ladder_max
        return band_patterns_discrepancy(
            other_bands, self.bands, ladder=self.ladder,
            relative_tolerance=relative_tolerance,
            zone=[min_band_cutoff, max_band_cutoff],
            reference_and_gel=True
        )

    def to_bandwagon_bandpattern(self, background_color=None, label='auto'):
        """Return a pattern version for the plotting library Bandwagon.

        If label is left to 'auto', it will be the pattern's name.
        """
        if label == 'auto':
            label = self.name
        return BandsPattern(
            self.bands,
            corner_note="Total: %d bp." % sum(self.bands),
            ladder=self.ladder,
            label=label,
            gel_image=self.migration_image,
            background_color=background_color
        )


class CloneValidation:
    """Report from comparing a clone's pband patterns with expectations.

    Instances of this class are produced by validating a ``Clone`` object.

    Parameters
    ----------
    clone
      A Clone object.

    expected
      A dictionnary ``{digestion: expected_pattern}`` where ``digestion`` is
      of the form ``('EcoRI', 'BamHI')`` and ``expected_pattern`` is a
      BandsPattern

    discrepancies
      A dictionnary ``{digestion: float}`` where the float indicates the
      observed discrepancy between observed and predicted for this digestion.

    """

    def __init__(self, clone, expected, discrepancies):
        """Initialize."""
        self.clone = clone
        self.expected = expected
        self.discrepancies = discrepancies

    @property
    def passes(self):
        """True iff the clone is valid.

        (This means that all patterns are similar enough to expectations).
        """
        return self.max_discrepancy < 1.0

    @property
    def max_discrepancy(self):
        """Return the largest discrepancy in all patterns"""
        return max(self.discrepancies.values())

    def color(self, discrepancy='max'):
        """Return a color depending on the discrepancy level (in 0-1)."""
        if discrepancy == "max":
            discrepancy = self.max_discrepancy
        factor = min(1.0, discrepancy ** 10)
        return ((2 + factor) / 3.0, (3 - factor) / 3.0, 2 / 3.0)

    def to_bandwagon_bandpattern(self, digestion, label='auto',
                                 per_digestion_discrepancy=True):
        """Return a (plottable) bandwagon pattern for the selected digestion.

        The pattern has a color indicating the match similarity, as well as
        a corner note indicating the gap to similarity.

        Parameters
        ----------
        digestion
          A set of enzymes such as ``('EcoRI',)`` or ````

        label
          A label for the pattern, for instance the name of a microplate well.

        per_digestion_discrepancy
          If true, each pattern for each clone will have an indication of
          how close it is from the expected pattern. If False, and the
          validation involves several digestions, the biggest discrepancy
          among all patterns is shown.

        """
        if digestion not in self.clone.digestions:
            return None
        if per_digestion_discrepancy:
            discrepancy = self.discrepancies[digestion]
        else:
            discrepancy = self.max_discrepancy
        obs = self.clone.digestions[digestion]
        pattern = obs.to_bandwagon_bandpattern(
            background_color=self.color(discrepancy=discrepancy), label=label)
        pattern.corner_note += " Gap: %d" % (100 * discrepancy)
        return pattern


class Clone:
    """Gather all informations necessary to validate a clone.

    Parameters
    ----------
    name
      Name of the clone. Could be for instance a microplate well.

    digestions
      A dictionnary ``{digestion: CloneObservation}`` where ``digestion`` is
      of the form ``('EcoRI', 'BamHI')``.

    construct_id
      ID of the construct to be validated. This is used to group clones
      by construct in the validation reports.
    """

    def __init__(self, name, digestions, construct_id=None):
        """Initialize."""
        self.name = name
        self.digestions = digestions
        self.construct_id = construct_id

    def validate_bands(self, bands_by_digestion, relative_tolerance=0.1,
                       min_band_cutoff=None, max_band_cutoff=None):
        """Return a validation results (comparison of observed and expected).

        The result is a CloneValidation object.

        Parameters
        ----------
        bands_by_digestion
          A dictionnary ``{digestion: [bands]}`` where ``digestion`` is
          of the form ``('EcoRI', 'BamHI')``, and bands is a list of band sizes

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
        for digestion, observation in self.digestions.items():
            discrepancies = {
                digestion: observation.patterns_discrepancy(
                    other_bands=bands_by_digestion[digestion],
                    relative_tolerance=relative_tolerance,
                    min_band_cutoff=min_band_cutoff,
                    max_band_cutoff=max_band_cutoff
                )
                for digestion, observation in self.digestions.items()
            }
        return CloneValidation(self, bands_by_digestion,
                               discrepancies=discrepancies)


class ClonesObservations:
    """All useful informations for a collection of clones to be validated.

    Parameters
    ----------
    clones
      Either a list of Clones (each with a unique name) or a dictionnary
      ``{name: Clone}``.

    constructs_dict
      A dictionnary ``{construct_name: biopython_record}`` indicating the
      sequence of the different constructs. For each construct, set the
      attribute ``construct.linear = False`` if the construct is circular.

    """

    def __init__(self, clones, constructs_dict):
        """Initialize"""
        if isinstance(clones, (list, tuple)):
            clones = {clone.name: clone for clone in clones}
        self.clones = clones
        self.constructs_dict = constructs_dict
        self.constructs_digestions = {ct: {} for ct in constructs_dict}

    def get_digestion_bands_for_construct(self, construct_id, digestion):
        """Return the bands resulting from the digestion of the construct.

        This function enables some memoization (no digestion is computed
        twice).

        Parameters
        ----------

        construct_id
          ID of the construct as it appears in this object'
          ``constructs_dict``.

        digestion
          For instance ``('BamHI', 'XbaI')``

        """
        construct_digestions = self.constructs_digestions[construct_id]
        if digestion not in construct_digestions:
            construct_record = self.constructs_dict[construct_id]
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
                for construct_id in self.constructs_dict
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

    def plot_validations_plate_map(self, validations, target_file=None,
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
        if target_file is not None:
            ax.figure.savefig(target_file, bbox_inches='tight')
            plt.close(ax.figure)
        else:
            return ax

    def plot_validations_summary(self, validations,
                                 target_file=None,
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
        digestions_by_construct = {construct: [] for construct in summary}
        for construct, validations in summary.items():
            for validation in validations:
                for digestion in validation.discrepancies:
                    if digestion not in digestions_by_construct[construct]:
                        digestions_by_construct[construct].append(digestion)
        lines_per_construct = [
            len(digestions)
            for c, digestions in digestions_by_construct.items()
        ]
        total_lines = sum(lines_per_construct)
        max_x = max(len(clones) for construct, clones in summary.items()) + 1

        fig, axes = plt.subplots(total_lines, 2, figsize=(2.2 * max_x,
                                                          3 * total_lines))
        axes_lines = (line for line in axes)
        ladder = list(validations[0].clone.digestions.values())[0].ladder
        for construct_id, validations in summary.items():
            digestions = digestions_by_construct[construct_id]
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
                patterns_set = BandsPatternsSet(
                    [reference] + patterns, ladder=ladder,
                    label="\n".join([construct_id, " + ".join(digestion)]),
                    ladder_ticks=5,
                    global_patterns_props={"label_fontdict": {"rotation": 60}}
                )

                record = self.constructs_dict[construct_id]

                ax_left, ax_right = next(axes_lines)
                ax_left.set_xlim(0.5, max_x + 2)
                patterns_set.plot(ax_left)
                ax_left.set_xlim(0.5, max_x + 2)
                self._plot_digestion(record, enzymes=digestion, ax=ax_right)

        fig.subplots_adjust(hspace=0.3)
        if target_file is not None:
            fig.savefig(target_file, bbox_inches="tight")
            plt.close(fig)
        else:
            return axes

    def write_identification_report(self, target_file=None,
                                    relative_tolerance=0.1,
                                    min_band_cutoff=None,
                                    max_band_cutoff=None):
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
        bands_validatities = self.compute_all_bands_validities(
            relative_tolerance=relative_tolerance,
            min_band_cutoff=min_band_cutoff,
            max_band_cutoff=max_band_cutoff)
        L = len(self.constructs_records)
        max_x = max(
            len(measures)
            for measures in self.observed_bands.values()
        ) + 1
        fig, axes = plt.subplots(L, 2, figsize=(2.2 * max_x, 3 * L))
        axes_validities = zip(axes, bands_validatities.items())
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
            self._plot_digestion(record, enzymes=self.digestions[construct_id],
                                 ax=ax2)
            patterns_set.plot(ax1)
        fig.subplots_adjust(hspace=0.3)
        if target_file is not None:
            fig.savefig(target_file, bbox_inches="tight")
            plt.close(fig)
        else:
            return axes

    def _create_digestion_graphic_record(self, record, enzymes):
        """Create a DnaFeaturesViewer graphic record showing the cuts."""
        batch = Restriction.RestrictionBatch(enzymes)
        cuts_dict = batch.search(record.seq)
        all_cuts = sorted(
            set([0, len(record)] +
                [c for cc in cuts_dict.values() for c in cc])
        )
        bands = list(zip(all_cuts, all_cuts[1:]))
        if not record.__dict__.get('linear', True):
            start, end = bands.pop()
            band0 = [-(end - start), bands[0][1]]
            if bands == []:
                bands = [band0]
            else:
                bands[0] = band0
        sorted_bands = sorted(
            bands, key=lambda start_end: start_end[0] - start_end[1])
        gr_cuts = GraphicRecord(len(record), [
            GraphicFeature(_start, _end,  strand=1, label=name,
                           color="#ede15c", thickness=5)
            for name, (_start, _end) in zip("abcdefghijklmn", sorted_bands)
        ])
        gr_cuts.split_overflowing_features_circularly()
        return gr_cuts, all_cuts

    def _plot_digestion(self, record, enzymes, ax):
        """Plot the record and its cuts one above the other."""

        gr_cuts, all_cuts = self._create_digestion_graphic_record(
            record, enzymes)
        gr_cuts.plot(ax, fontsize=7, with_ruler=False)
        for cut in all_cuts:
            ax.axvline(cut, ls=":", color="k", lw=0.5)

        def features_prop(f):
            return dict(label=f.qualifiers.get("source", [False])[0],
                        color="#93a8ea",
                        thickness=5)

        def features_filter(f):
            return features_prop(f)["label"]

        translator = BiopythonTranslator([features_filter], features_prop)
        gr_record = translator.translate_record(record)
        gr_record.plot(ax, fontsize=4, level_offset=7)
