"""Classes for parsing and analyzis of DNA migration data."""
from collections import OrderedDict, Counter
import matplotlib.pyplot as plt
from io import BytesIO
import pandas
from matplotlib.backends.backend_pdf import PdfPages

from bandwagon import BandsPattern, BandsPatternsSet
import flametree

try:
    from plateo.exporters.plate_to_matplotlib_plots import PlateColorsPlotter
    from plateo.parsers import plate_from_platemap_spreadsheet
    from plateo.containers import Plate96

    PLATEO_AVAILABLE = True
except ImportError:
    PLATEO_AVAILABLE = False

try:
    import saboteurs

    SABOTEURS_AVAILABLE = True
except ImportError:
    SABOTEURS_AVAILABLE = False

from ..bands_predictions import predict_digestion_bands
from ..tools import load_record, all_subsets, record_is_linear
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

    def __init__(self, clones, constructs_records, partial_cutters=()):
        """Initialize"""
        if isinstance(clones, (list, tuple)):
            clones = OrderedDict([(clone.name, clone) for clone in clones])
        self.clones = clones
        clones_constructs = []
        for clone in list(clones.values()):
            if clone.construct_id not in clones_constructs:
                clones_constructs.append(clone.construct_id)
        self.constructs_records = OrderedDict(
            sorted(
                [
                    (construct_name, record)
                    for construct_name, record in constructs_records.items()
                ],
                key=lambda cst: +100000
                if cst[0] not in clones_constructs
                else clones_constructs.index(cst[0]),
            )
        )
        self.constructs_digestions = {ct: {} for ct in constructs_records}
        self.partial_cutters = partial_cutters

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
                linear=record_is_linear(construct_record, default=True),
                enzymes=[
                    enzyme
                    for enzyme in digestion
                    if enzyme not in self.partial_cutters
                ],
                partial_cutters=self.partial_cutters,
            )
        return construct_digestions[digestion]

    def get_clone_digestion_bands(self, clone, construct_id=None):
        """Return ``{digestion: bands}`` for all digestions of this clone."""
        if construct_id is None:
            construct_id = clone.construct_id
        return {
            digestion: self.get_digestion_bands_for_construct(
                construct_id=clone.construct_id, digestion=digestion
            )
            for digestion in clone.digestions
        }

    def validate_all_clones(
        self,
        relative_tolerance=0.05,
        min_band_cutoff=None,
        max_band_cutoff=None,
    ):
        """Return ``{clone: CloneValidation}`` for all clones."""
        return OrderedDict(
            [
                (
                    clone_name,
                    clone.validate_bands(
                        bands_by_digestion=self.get_clone_digestion_bands(
                            clone
                        ),
                        relative_tolerance=relative_tolerance,
                        min_band_cutoff=min_band_cutoff,
                        max_band_cutoff=max_band_cutoff,
                    ),
                )
                for clone_name, clone in self.clones.items()
            ]
        )

    def identify_all_clones(
        self,
        relative_tolerance=0.05,
        min_band_cutoff=None,
        max_band_cutoff=None,
    ):
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
        return OrderedDict(
            [
                (
                    clone_name,
                    {
                        construct_id: clone.validate_bands(
                            bands_by_digestion=self.get_clone_digestion_bands(
                                clone, construct_id=construct_id
                            ),
                            relative_tolerance=relative_tolerance,
                            min_band_cutoff=min_band_cutoff,
                            max_band_cutoff=max_band_cutoff,
                        )
                        for construct_id in self.constructs_records
                    },
                )
                for clone_name, clone in self.clones.items()
            ]
        )

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
                results[cst_id] = sorted(
                    validations, key=lambda v: v.max_discrepancy
                )
        return results

    def plot_validations_plate_map(self, validations, target=None, ax=None):
        """Plot a map of the plate with passing/failing wells in green/red."""

        if not PLATEO_AVAILABLE:
            raise ImportError("This function requires Plateo installed")
        plate = Plate96("Validations map")
        constructs_colors = {}
        for well in plate.iter_wells():
            if well.name not in validations:
                well.data.bg_color = None
                well.data.pass_color = (1, 1, 1, 0.1)
            else:
                validation = validations[well.name]
                construct = validation.clone.construct_id
                if construct not in constructs_colors:
                    alpha = ((0.35 * len(constructs_colors) % 1) + 0.2) / 2
                    constructs_colors[construct] = (0.4, 0.4, 0.76, alpha)
                well.data.bg_color = constructs_colors[construct]
                well.data.pass_color = (
                    (0.43, 0.92, 0.49)
                    if validation.passes
                    else (0.91, 0.43, 0.43)
                )
        bg_plotter = PlateColorsPlotter(lambda w: w.data.bg_color)
        ax, _ = bg_plotter.plot_plate(plate, figsize=(10, 5), ax=ax)
        pass_plotter = PlateColorsPlotter(
            lambda w: w.data.pass_color, well_radius=250
        )
        pass_plotter.plot_plate(plate, ax=ax)
        if target is not None:
            ax.figure.savefig(target, bbox_inches="tight")
            plt.close(ax.figure)
        else:
            return ax

    def validations_summary_table(self, validations, target=None):
        validations_summary = self.validations_summary(validations)
        records = []
        for construct, data in validations_summary.items():
            best_clone = None
            valid_clones = [clone for clone in data if clone.passes]
            best_clone_name = None
            if len(valid_clones):
                best_clone = min(valid_clones, key=lambda c: c.max_discrepancy)
                best_clone_name = best_clone.clone.name
            other_valid_clones = ", ".join(
                [
                    clone.clone.name
                    for clone in valid_clones
                    if clone.clone.name != best_clone_name
                ]
            )
            records.append(
                dict(
                    construct=construct,
                    n_clones=len(data),
                    valid_clones=len(valid_clones),
                    best_clone=best_clone_name,
                    other_valid_clones=other_valid_clones,
                )
            )
        dataframe = pandas.DataFrame.from_records(
            records,
            columns=[
                "construct",
                "n_clones",
                "valid_clones",
                "best_clone",
                "other_valid_clones",
            ],
        )
        if target is not None:
            dataframe.to_csv(target, index=False)
        return dataframe

    def plot_all_validations_patterns(
        self, validations, target=None, per_digestion_discrepancy=False
    ):
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
        max_x = (
            Counter(
                [clone.construct_id for clone in self.clones.values()]
            ).most_common(1)[0][1]
            + 2
        )
        digestions_by_construct = {construct: [] for construct in summary}
        for construct, validations in summary.items():
            for validation in validations:
                for digestion in validation.discrepancies:
                    if digestion not in digestions_by_construct[construct]:
                        digestions_by_construct[construct].append(digestion)

        # ladder = list(validations[0].clone.digestions.values())[0].ladder
        pdf_io = BytesIO()
        with PdfPages(pdf_io) as pdf:
            for construct_id, validations in summary.items():
                ladder = list(validations[0].clone.digestions.values())[
                    0
                ].ladder
                digestions = digestions_by_construct[construct_id]
                band_patterns = OrderedDict()
                seen_clones = set()
                for i, digestion in enumerate(digestions):
                    reference_bands = self.get_digestion_bands_for_construct(
                        construct_id, digestion
                    )
                    reference = BandsPattern(
                        bands=reference_bands,
                        ladder=ladder,
                        label="exp." if (i == 0) else None,
                        background_color="#c6dcff",
                        corner_note="Total: %d bp" % sum(reference_bands),
                        global_bands_props={"label_fontdict": {"size": 5}},
                    )
                    sorted_bands = sorted(
                        reference.bands, key=lambda b: -b.dna_size
                    )
                    for band_name, band in zip("abcdefghijklm", sorted_bands):
                        band.label = band_name

                    patterns = []
                    for validation in validations:
                        clone_name = validation.clone.name
                        label = None if clone_name in seen_clones else "auto"
                        pattern = validation.to_bandwagon_bandpattern(
                            digestion=digestion,
                            label=label,
                            per_digestion_discrepancy=per_digestion_discrepancy,
                        )
                        if pattern is not None:
                            seen_clones.add(clone_name)
                        patterns.append(pattern)
                    band_patterns[digestion] = BandsPatternsSet(
                        [reference] + patterns,
                        ladder=ladder,
                        label=" + ".join(digestion),
                        ladder_ticks=5,
                        global_patterns_props={
                            "label_fontdict": {"rotation": 60}
                        },
                    )

                digestions = digestions_by_construct[construct_id]
                fig, axes = plt.subplots(
                    len(digestions),
                    1,
                    sharex=True,
                    figsize=(1.1 * max_x, 3 * len(digestions)),
                )
                if len(digestions) == 1:
                    axes = [axes]
                axes[-1].set_xlim(right=max_x)

                for ax, (dig, pattern_set) in zip(axes, band_patterns.items()):
                    pattern_set.plot(ax)

                axes[-1].set_xlabel(
                    construct_id,
                    fontdict=dict(size=14, weight="bold"),
                    ha="left",
                    position=(0.0, 0.1),
                )
                pdf.attach_note(construct_id)

                # Temporary fix because of Matplotlib bug #12634
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)
        pdf_data = pdf_io.getvalue()
        if target is None:
            return pdf_data
        elif target == "base64":
            return "data:application/pdf;base64," + pdf_data.decode("utf-8")
        else:
            with open(target, "wb") as f:
                f.write(pdf_data)

    @staticmethod
    def identify_bad_parts(
        validations,
        constructs_parts,
        constructs_records=None,
        report_target=None,
        extra_failures=None,
    ):
        """Identifies parts associated with failure in the validations.

        Uses the Saboteurs library:
        https://github.com/Edinburgh-Genome-Foundry/saboteurs

        Parameters
        ----------
        validations
          validations results

        constructs_parts
          Either a dict {construct_id: [list, of, part, names]} or a function
          (biopython_record => [list, of, part, names])

        report_target
          Can be a path to file or file-like object where to write a PDF
          report. Can also be "@memory", at which case the raw binary PDF
          data is returned

        Returns
        -------

        analysis, pdf_data
          Where ``analysis`` is the result of sabotage analysis (see the
          saboteurs library), and pdf_data is None unless report_target is
          set to "@memory" (see above).
        """
        if not SABOTEURS_AVAILABLE:
            raise ImportError(
                "You must install the saboteurs library to"
                "be able to identify bad parts. Try:\n\n"
                "(sudo) pip install saboteurs"
            )
        if hasattr(constructs_parts, "__call__"):
            constructs_parts = {
                record_id: constructs_parts(record)
                for record_id, record in constructs_records.items()
            }
        extra_failures = extra_failures or {}
        constructs_stats = OrderedDict()
        for clone_name, validation in validations.items():
            construct_id = validation.clone.construct_id
            if construct_id not in constructs_stats:
                constructs_stats[construct_id] = dict(
                    exp_id=construct_id,
                    attempts=0,
                    failures=0,
                    members=constructs_parts[construct_id],
                )
            stats = constructs_stats[construct_id]
            stats["attempts"] += 1
            stats["failures"] += not validation.passes
        for construct_id, weight in extra_failures.items():
            constructs_stats[construct_id] = dict(
                exp_id=construct_id,
                attempts=weight,
                failures=weight,
                members=constructs_parts[construct_id],
            )
        # import json
        # print (json.dumps(constructs_stats, indent=2))
        analysis = saboteurs.find_statistical_saboteurs(constructs_stats)
        report_data = None
        if report_target is not None:
            report_data = saboteurs.statistics_report(analysis, report_target)
        return analysis, report_data

    def write_identification_report(
        self,
        target_file=None,
        relative_tolerance=0.05,
        min_band_cutoff=None,
        max_band_cutoff=None,
    ):
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
            max_band_cutoff=max_band_cutoff,
        )
        L = len(self.constructs_records)
        max_x = (
            max(len(measures) for measures in self.observed_bands.values()) + 1
        )
        fig, axes = plt.subplots(L, 1, figsize=(2.2 * max_x, 3 * L))
        axes_validities = zip(axes, bands_validities.items())
        for ax, (construct_id, validities) in axes_validities:
            reference = BandsPattern(
                self.expected_bands[construct_id],
                ladder=self.ladder,
                label="exp.",
                background_color="#c6dcff",
                corner_note="Total: %d bp"
                % sum(self.expected_bands[construct_id]),
                global_bands_props={"label_fontdict": {"size": 5}},
            )
            sorted_bands = sorted(reference.bands, key=lambda b: -b.dna_size)
            for band_name, band in zip("abcdefghijklm", sorted_bands):
                band.label = band_name
            patterns = [
                BandsPattern(
                    self.observed_bands[construct_id][measure_name],
                    corner_note="Total: %d bp"
                    % sum(self.observed_bands[construct_id][measure_name]),
                    ladder=self.ladder,
                    label=measure_name,
                    gel_image=self.migration_images[measure_name],
                    background_color="#aaffaa"
                    if validities[measure_name]
                    else "#ffaaaa",
                )
                for measure_name in self.observed_bands[construct_id]
            ]
            patterns_set = BandsPatternsSet(
                [reference] + patterns,
                ladder=self.ladder,
                label=construct_id,
                ladder_ticks=5,
                global_patterns_props={"label_fontdict": {"rotation": 60}},
            )
            ax.set_xlim(0.5, max_x + 2)
            patterns_set.plot(ax)
        fig.subplots_adjust(hspace=0.3)
        if target_file is not None:
            fig.savefig(target_file, bbox_inches="tight")
            plt.close(fig)
        else:
            return axes

    @staticmethod
    def from_files(
        records_path,
        constructs_map_path,
        aati_zip_path,
        digestions_map_path=None,
        digestion=None,
        direction="column",
        ignore_bands_under=None,
        topology="circular"
    ):
        constructs_records = {
            f._name_no_extension: load_record(
                f.open("r"), topology=topology, id=f._name_no_extension,
                file_format='genbank'
            )
            for f in flametree.file_tree(records_path)._all_files
        }
        constructs_records = OrderedDict(sorted(constructs_records.items()))

        constructs_plate = plate_from_platemap_spreadsheet(
            constructs_map_path, data_field="construct", headers=True
        )
        constructs_map = OrderedDict(
            [
                (well.name, well.data.construct)
                for well in constructs_plate.iter_wells(direction=direction)
                if "construct" in well.data
                and str(well.data.construct) != "nan"
            ]
        )

        if digestion:
            digestion = tuple(digestion)
            digestions_map = OrderedDict(
                [(wellname, digestion) for wellname in constructs_map]
            )
        else:
            digestions_plate = plate_from_platemap_spreadsheet(
                digestions_map_path, data_field="digestion", headers=True
            )
            digestions_map = OrderedDict(
                [
                    (well.name, tuple(well.data.digestion.split(", ")))
                    for well in digestions_plate.iter_wells(
                        direction=direction
                    )
                    if "digestion" in well.data
                    and str(well.data.digestion) != "nan"
                ]
            )

        observations = BandsObservation.from_aati_fa_archive(
            aati_zip_path,
            direction=direction,
            ignore_bands_under=ignore_bands_under,
        )
        clones = Clone.from_bands_observations(
            observations, constructs_map, digestions_map
        )
        clones_observations = ClonesObservations(clones, constructs_records)
        return clones_observations

    def partial_digests_analysis(self, relative_tolerance=0.05):
        """Compute good clones under different partial digest assumptions.

        Returns a dictionnary ``{partial: {'valid_clones': 60, 'label': 'x'}}``
        where for a given scenario ``partial`` is a tuple of all enzymes
        considered to have partia activity in this scenario (it defines the
        scenario) ``valid_clones`` is the number of good clones under this
        assumption, and ``label`` is a string representation of all enzymes
        involved, with partial activity enzymes in parenthesis.

        This result can be fed to ``ClonesObservations``'s
        ``.plot_partial_digests_analysis`` method for plotting
        """
        all_enzymes = set(
            enzyme
            for clone in self.clones.values()
            for digestion in clone.digestions
            for enzyme in digestion
        )
        all_enzymes
        results = {}
        for good_cutters in all_subsets(all_enzymes):
            partial_cutters = tuple(
                [e for e in all_enzymes if e not in good_cutters]
            )
            clones_observations = ClonesObservations(
                self.clones,
                self.constructs_records,
                partial_cutters=partial_cutters,
            )
            validations = clones_observations.validate_all_clones(
                relative_tolerance=relative_tolerance
            )
            label = " + ".join(
                sorted(good_cutters)
                + sorted(["(%s)" % c for c in partial_cutters])
            )
            valid_clones = sum([v.passes for name, v in validations.items()])

            results[partial_cutters] = {
                "label": label,
                "validations": validations,
                "valid_clones": valid_clones,
            }
        return results

    @staticmethod
    def plot_partial_digests_analysis(analysis_results, ax=None):
        """Plot partial digests analysis results.

        Parameters
        ----------
        analysis_results
          results from ``ClonesObservations.partial_digest_analysis``

        ax
          A Matplotlib ax. If none, one is created and returned at the end.

        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(5, 0.5 * len(analysis_results)))
        ax.axis("off")
        values, labels = zip(
            *sorted(
                [
                    (data["valid_clones"], data["label"])
                    for partial, data in analysis_results.items()
                ]
            )
        )
        ax.barh(range(len(values)), values)
        for i, (value, label) in enumerate(zip(values, labels)):
            ax.text(
                0,
                i,
                label + " ",
                ha="right",
                va="center",
                fontweight=("normal" if "(" in label else "bold"),
            )
            ax.text(
                value,
                i,
                str(value) + " ",
                ha="right",
                va="center",
                fontdict={"color": "white", "weight": "bold"},
            )
        ax.set_ylim(-0.5, len(values) - 0.5)
        ax.set_title("Number of good clones, by scenario")
        return ax

    @staticmethod
    def merge_observations(observations_list, prefixes=None):
        prefixes = prefixes or (len(observations_list) * [""])

        clones = OrderedDict(
            [
                (pref + k, v)
                for pref, obs in zip(prefixes, observations_list)
                for k, v in obs.clones.items()
            ]
        )
        constructs_records = {
            name: record
            for obs in observations_list
            for name, record in obs.constructs_records.items()
        }
        partial_cutters = tuple(
            enzyme
            for obs in observations_list
            for enzyme in obs.partial_cutters
        )

        return ClonesObservations(clones, constructs_records, partial_cutters)
