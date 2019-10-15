"""Implements the SeparatingDigestions- and IdealDigestionProbem classes."""


from collections import OrderedDict
from copy import deepcopy

import numpy as np

from ..tools import (
    updated_dict,
    sequence_to_biopython_record,
    record_is_linear,
    set_record_topology,
)
from .SetCoverProblem import SetCoverProblem
from ..bands_predictions import (
    predict_sequence_digestions,
    compute_sequence_digestions_migrations,
)

try:
    import bandwagon
    import matplotlib.pyplot as plt

    PLOTS_AVAILABLE = True
except ImportError:
    PLOTS_AVAILABLE = False


class DigestionProblem(SetCoverProblem):
    """General class for solving digestions problems.

    By digestion problem we mean *enzymatic* digestion problems. For other
    kinds of digestion problems, see a proctologist.

    Other digestion problems subclass this problem for particular cases such
    as construct validation and identification.

    Parameters
    ----------
    enzymes
      List of the names of the enzymes to consider, e.g. ``['EcoRI', 'XbaI']``.

    ladder
      A Ladder object representing the ladder used for migrations.

    sequences
      An (ordered) dictionary of the form {sequence_name: sequence} where the
      sequence is an ATGC string

    linear
      True for linear sequences, false for circular sequences

    max_enzymes_per_digestion
      Maximal number of enzymes that can go in a single digestion.
      Experimentally the best is 1, but you can try 2, or 3 in desperate
      situations. Not sure if more enzymes will work.

    relative_migration_precision
      Variance of the bands measured during the migration, given as a
      proportion of the total migration span (difference between the migration
      of the ladder's smallest and largest bands).

    """

    def __init__(
        self,
        enzymes,
        ladder,
        sequences,
        topology="auto",
        default_topology="linear",
        max_enzymes_per_digestion=1,
        relative_migration_precision=0.1,
        progress_logger=None,
    ):
        """Initialize."""
        if hasattr(sequences, "items"):
            if hasattr(list(sequences.values())[0], "seq"):
                self.records = {
                    name: deepcopy(record)
                    for name, record in sequences.items()
                }
                self.sequences = {
                    name: str(record.seq) for name, record in sequences.items()
                }
            else:
                self.sequences = sequences
                self.records = {
                    name: sequence_to_biopython_record(seq)
                    for name, seq in sequences.items()
                }
        else:
            self.records = OrderedDict(
                [(r.id, deepcopy(r)) for r in sequences]
            )
            self.sequences = OrderedDict(
                [(r.id, str(r.seq)) for r in sequences]
            )
        for record in self.records.values():
            if topology == "auto":
                set_record_topology(
                    record, default_topology, pass_if_already_set=True
                )
            else:
                set_record_topology(record, topology)
        self.relative_migration_precision = relative_migration_precision
        self.max_enzymes_per_digestion = max_enzymes_per_digestion

        self.sequences_names = list(self.sequences.keys())
        self.progress_logger = progress_logger
        self.enzymes = enzymes

        def part_is_linear(name):
            if topology == "auto":
                default = default_topology == "linear"
                return record_is_linear(self.records[name], default=default)
            else:
                return topology == "linear"

        self.sequences_digestions = {
            name: predict_sequence_digestions(
                sequence=sequence,
                enzymes=self.enzymes,
                linear=record_is_linear(self.records[name]),
                max_enzymes_per_digestion=self.max_enzymes_per_digestion,
            )
            for name, sequence in self.sequences.items()
        }
        digestions = self.sequences_digestions[self.sequences_names[0]]
        self.digestions = list(digestions.keys())
        self.change_ladder(ladder)

    def change_ladder(self, ladder):
        self.ladder = ladder
        mini, maxi = self.ladder.migration_distances_span
        self.migration_min, self.migration_max = mini, maxi
        self.migration_span = maxi - mini
        compute_sequence_digestions_migrations(
            self.sequences_digestions, self.ladder
        )
        SetCoverProblem.__init__(
            self,
            elements=self._compute_elements(),
            parameters=self.digestions,
            progress_logger=self.progress_logger,
        )

    @staticmethod
    def _default_heuristic(named_subset, selected):
        """Select for coverage first, digestion complexity second.

        When two digestions cover the same number of elements, the digestion
        with the less enzymes is preffered.

        """
        enzymes, subset = named_subset
        return len(subset) - 0.5 * len(enzymes)

    def bands_to_migration_pattern(self, bands_sizes):
        """Return the distance migrations from several bands sizes."""
        return self.ladder.dna_size_to_migration(np.array(bands_sizes))

    def plot_digestions(
        self,
        digestions,
        axes=None,
        bands_props=None,
        patterns_props=None,
        patternset_props=None,
        target_file=None,
        close_figure=False,
    ):
        """Plot the patterns for each sequence, for each digestion in the list.

        Requires Bandwagon.

        Parameters
        ----------
        digestions
          A list of digestions, e.g. [('EcoRI',), ('BamHI', 'XbaI'), ...]

        axes
          Axes on which to plot the plot. There should be as many axes in the
          list as there are digestions provided. Ideally the axes would be
          vertically stacked. If not provided, new axes are created and
          returned in the end

        bands_props, patterns_props, patternset_props
          Graphical properties (colors, labels, etc.) of band patterns, see
          code and BandWagon for more details.

        target_file
          The name of the (PNG, SVG, JPEG...) file in which to write the plot.

        Returns
        -------

        axes
          The axes of the generated figure (if a target file is written to,
          the figure is closed and None is returned instead)

        """
        if not PLOTS_AVAILABLE:
            raise ImportError("Plotting requires Bandwagon installed")

        bands_props = updated_dict({"size": 10, "rotation": 50}, bands_props)
        patternset_props = updated_dict({"ladder_ticks": 5}, patternset_props)
        if isinstance(self.ladder.bands, (list, tuple)):
            ladder = self.ladder
        else:
            ladder = bandwagon.custom_ladder(None, self.ladder.bands)
        ndig, nseq = len(digestions), len(self.sequences)
        if axes is None:
            fig, axes = plt.subplots(ndig, 1, figsize=(0.5 * nseq, 3 * ndig))
        if ndig == 1:
            axes = [axes]
        for ax, digestion in zip(axes, digestions):
            bandwagon.BandsPatternsSet(
                patterns=[
                    ladder if ax == axes[0] else ladder.modified(label=None)
                ]
                + [
                    bandwagon.BandsPattern(
                        bands=self.sequences_digestions[seq_name][digestion][
                            "bands"
                        ],
                        label=seq_name if (ax == axes[0]) else None,
                        ladder=ladder,
                        circularity=self.records[seq_name].annotations[
                            "topology"
                        ],
                        global_bands_props=bands_props,
                    )
                    for seq_name in self.sequences
                ],
                ladder=ladder,
                label=", ".join(
                    [" + ".join([e for e in d]) for d in [digestion]]
                ),
                global_patterns_props=patterns_props,
                **patternset_props
            ).plot(ax)
        if target_file is not None:
            axes[0].figure.savefig(target_file, bbox_inches="tight")
        if close_figure:
            plt.close(axes[0].figure)
        return axes

    def select_digestions(
        self,
        minimal_score=None,
        max_digestions=None,
        search="greedy",
        min_score_precision=0.001,
    ):
        """Select one or more digestions which collectively solve the problem.

        This lets you either find the highest-score solution of less than N
        digestions, or the best free-sized set of digestions solving the
        problem with the highest minimal score.

        Returns
        -------
        min_score, [digestion1, digestion2, ...]
          Where min_score is the lowest score that an element (sequence or
          sequence pair) gets from the coverage, and the list indicates the
          digestions in the order in which they have been selected. This list
          will be "None" for unsolvable problems.


        Parameters
        -----------
        minimal_score
          If provided and max_digestions is not provided, the method will
          return a free-sized set of parameters covering all elements.
          If max_digestions is provided, the subset will not exceed
          max_digestions in size. This can result in unsolvable problems.

        max_digestions
          When provided instead of minimal_score, the method will find the
          highest minimal_score which allows solutions in max_digestions or
          less.

        search
          Either 'greedy' for fast approximate solving or 'full' for full
          solving

        min_score_precision
          When max_digestions is provided and not minimal_score, the returned
          solution will have the optimal minimal_score, plus or minus this
          tolerance.

        """
        return self._select_parameters(
            threshold=minimal_score,
            covering_algorithm=search,
            max_set_size=max_digestions,
            heuristic="default",
            threshold_tolerance=min_score_precision,
            bisection=True,
        )
