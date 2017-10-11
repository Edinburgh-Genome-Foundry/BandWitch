"""Implements the SeparatingDigestions- and IdealDigestionProbem classes."""

import itertools
from collections import OrderedDict

import numpy as np
from .tools import (max_min_distance, updated_dict,
                    digestions_list_to_string)
from .SetCoverProblem import SetCoverProblem
from .bands_predictions import predict_sequence_digestions

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

    def __init__(self, enzymes, ladder, sequences, linear=False,
                 max_enzymes_per_digestion=1,
                 relative_migration_precision=0.1,
                 progress_logger=None):
        """Initialize."""
        if isinstance(sequences, (list, tuple)):
            sequences = OrderedDict([(r.id, str(r.seq)) for r in sequences])
        self.sequences = sequences
        self.ladder = ladder
        self.linear = linear
        self.relative_migration_precision = relative_migration_precision
        self.max_enzymes_per_digestion = max_enzymes_per_digestion

        self.sequences_names = list(sequences.keys())
        mini, maxi = self.ladder.migration_distances_span
        self.migration_min, self.migration_max = mini, maxi
        self.migration_span = maxi - mini

        self.enzymes = enzymes

        self.sequences_digestions = {
            name: predict_sequence_digestions(
                sequence=sequence, enzymes=self.enzymes, linear=self.linear,
                max_enzymes_per_digestion=self.max_enzymes_per_digestion,
                bands_to_migration=self.bands_to_migration_pattern)
            for name, sequence in self.sequences.items()
        }

        digestions = self.sequences_digestions[self.sequences_names[0]]
        self.digestions = parameters = list(digestions.keys())
        elements = self._compute_elements()
        SetCoverProblem.__init__(self, elements=elements, parameters=parameters,
                                 progress_logger=progress_logger)

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

    def plot_digestions(self, digestions, axes=None, bands_props=None,
                        patterns_props=None, patternset_props=None,
                        target_file=None):
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
                ] + [
                    bandwagon.BandsPattern(
                        bands=self.sequences_digestions[seq_name]
                                                       [digestion]['bands'],
                        label=seq_name if (ax == axes[0]) else None,
                        ladder=ladder,
                        global_bands_props=bands_props
                    )
                    for seq_name in self.sequences
                ],
                ladder=ladder,
                label=digestions_list_to_string([digestion]),
                global_patterns_props=patterns_props,
                **patternset_props
            ).plot(ax)
        if target_file is not None:
            axes[0].figure.savefig(target_file, bbox_inches='tight')
            plt.close(axes[0].figure)
        return axes

    def select_digestions(self, minimal_score=None, max_digestions=None,
                          search='greedy', min_score_precision=0.001):
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
            heuristic='default',
            threshold_tolerance=min_score_precision,
            bisection=True)


class SeparatingDigestionsProblem(DigestionProblem):
    """Problem: find best digestion(s) to identify constructs.

    Provided a set of constructs (possibly from a combinatorial assembly),
    find a list of digestions such that all pair of constructs have different
    migration patterns for at least one digestion of the list.

    """

    """General class for solving digestions problems.

    By digestion problem we mean *enzymatic* digestion problems. For other
    kinds of digestion problems, consult your proctologist.

    Other digestion problems subclass this problem and implement a computation
    of coverages which depends on the problem.

    Parameters
    ----------

    sequences
      An (ordered) dictionary of the form {sequence_name: sequence} where the
      sequence is an ATGC string

    enzymes
      List of the names of the enzymes to consider, e.g. ``['EcoRI', 'XbaI']``.

    ladder
      A Ladder object representing the ladder used for migrations.

    linear
      True for linear sequences, false for circular sequences

    max_enzymes_per_digestion
      Maximal number of enzymes that can go in a single digestion.
      Experimentally the best is 1, but you can try 2, or 3 in desperate
      situations. Not sure if more enzymes will work.

    relative_error
      Variance of the bands measured during the migration, given as a
      proportion of the total migration span (difference between the migration
      of the ladder's smallest and largest bands).
    """

    def __init__(self, enzymes, ladder, sequences=None, categories=None,
                 linear=False, max_enzymes_per_digestion=1,
                 min_discrepancy='auto',
                 relative_migration_precision=0.1):
        """Initialize."""
        if isinstance(sequences, (list, tuple)):
            sequences = OrderedDict([(r.id, str(r.seq)) for r in sequences])

        if categories is None:
            categories = OrderedDict([
                (seq_name, {seq_name: sequence})
                for seq_name, sequence in sequences.items()
            ])
        else:
            sequences = OrderedDict([
                (seq_name, sequence)
                for category, seqs in categories.items()
                for seq_name, sequence in seqs.items()
            ])
        self.categories = categories

        DigestionProblem.__init__(
            self, sequences=sequences, enzymes=enzymes, ladder=ladder,
            linear=linear, max_enzymes_per_digestion=max_enzymes_per_digestion,
            relative_migration_precision=relative_migration_precision)

    def _compute_elements(self):
        category_pairs = itertools.combinations(self.categories.values(), 2)
        return set(
            (seq1, seq2)
            for category1, category2 in category_pairs
            for seq1 in category1
            for seq2 in category2
        )

    def _parameter_element_score(self, digestion, sequences_pair):
        """See max_patterns_difference."""
        sequence1, sequence2 = sequences_pair
        digestion1, digestion2 = [self.sequences_digestions[s][digestion]
                                  for s in (sequence1, sequence2)]
        if ((digestion1['same_as'] == digestion2['same_as'])
                and digestion1['same_as'] in self.scores[sequences_pair]):
            return self.scores[sequences_pair][digestion1['same_as']]
        #     () and
        # same_as = tuple(sorted(digestion1['same_as'], digestion2['same_as']))

        mini, maxi = self.migration_min, self.migration_max
        margin = self.relative_migration_precision * self.migration_span / 2.0
        distance = max_min_distance(digestion1['migration'],
                                    digestion2['migration'],
                                    zone=[mini + margin, maxi - margin])
        return 1.0 * distance / self.migration_span

    @staticmethod
    def _score_to_color(score, maxi=0.1):
        """Transform a similarity score to a green/red color.

        Parameters
        ----------
        score
          Value between 0 (perfect similarity, green) and 1 (red)

        maxi
          Value of the score above which everything appears completely red.
          Below this value the color goes progressively from red to green in 0.

        """
        return (max(0, min(1, score / maxi)),
                min(1, max(0, 1 - score / maxi)), 0, .5)

    def plot_distances_map(self, digestions, ax=None, target_file=None):
        """Plot how well the digestions separate each construct pair.

        Parameters
        ----------
        digestions
          A list of digestions, eg ``[('EcoRV'), ('XbaI', 'MfeI')]``.

        ax
          A matplotlib ax on which to plot, if none is provided, one is created
          and returned at the end.

        target_file
          The name of the (PNG, SVG, JPEG...) file in which to write the plot.

        Returns
        -------

        axes
          The axes of the generated figure (if a target file is written to,
          the figure is closed and None is returned instead)

        """

        if not PLOTS_AVAILABLE:
            raise ImportError("Plots require Matplotlib/Bandwagon installed.")
        grid = np.zeros(2 * (len(self.sequences),))
        for i, s1 in enumerate(self.sequences):
            for j, s2 in enumerate(self.sequences):
                if i >= j:
                    grid[i, j] = np.nan
                else:
                    scores = [
                        self._parameter_element_score(digestion, (s1, s2))
                        for digestion in digestions
                    ]
                    grid[i, j] = max(scores)
        if ax is None:
            _, ax = plt.subplots(1, figsize=2 * (0.8 * len(grid),))
        ax.imshow(grid[:, ::-1], interpolation='nearest', cmap='OrRd_r',
                  vmin=0, vmax=0.2)
        for i in range(len(grid)):
            for j in range(len(grid)):
                if i > j:
                    ax.text(
                        len(self.sequences) - i - 1, j,
                        "%d%%" % (100 * grid[j, i]),
                        fontdict=dict(color='black', weight='bold', size=14),
                        horizontalalignment='center',
                        verticalalignment='center'
                    )

        ax.set_yticks(range(len(grid)))
        ax.set_yticklabels(list(self.sequences)[:-1], size=14,
                           fontdict={'weight': 'bold'})
        ax.set_xticks(range(len(grid)))
        ax.xaxis.set_ticks_position('top')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_xticklabels([" " + s for s in list(self.sequences)[1:][::-1]],
                           rotation=90, size=14, fontdict={'weight': 'bold'})
        ax.set_xlim(-0.5, len(self.sequences) - 1.5)
        ax.set_ylim(len(self.sequences) - 1.5, -0.5)

        if target_file is not None:
            ax.figure.savefig(target_file, bbox_inches='tight')
            plt.close(ax.figure)
        return ax


class IdealDigestionsProblem(DigestionProblem):
    """Find ideal digestion(s) to validate constructs.

    Other digestion problems subclass this problem and implement a computation
    of coverages which depends on the problem.

    Parameters
    ----------
    sequences
      An (ordered) dictionary of the form {sequence_name: sequence} where the
      sequence is an ATGC string

    enzymes
      List of the names of the enzymes to consider, e.g. ``['EcoRI', 'XbaI']``.

    ladder
      A Ladder object representing the ladder used for migrations.

    linear
      True for linear sequences, false for circular sequences

    max_enzymes_per_digestion
      Maximal number of enzymes that can go in a single digestion.
      Experimentally the best is 1, but you can try 2, or 3 in desperate
      situations. Not sure if more enzymes will work.

    relative_migration_error
      Variance of the bands measured during the migration, given as a
      proportion of the total migration span (difference between the migration
      of the ladder's smallest and largest bands).

    """

    def __init__(self, enzymes, ladder, sequences, min_bands=3, max_bands=7,
                 border_tolerance=0.1, linear=False,
                 max_enzymes_per_digestion=1,
                 relative_migration_precision=0.1):
        """Initialize."""
        self.min_bands = min_bands
        self.max_bands = max_bands
        self.border_tolerance = border_tolerance

        DigestionProblem.__init__(
            self, sequences=sequences, enzymes=enzymes, ladder=ladder,
            linear=linear, max_enzymes_per_digestion=max_enzymes_per_digestion,
            relative_migration_precision=relative_migration_precision)

    def _parameter_element_score(self, digestion, sequence):
        """Compute the sequence's ``.migration_score`` for each digestion."""
        digestion = self.sequences_digestions[sequence][digestion]
        if digestion['same_as'] in self.scores[sequence]:
            return self.scores[sequence][digestion['same_as']]
        migration = digestion['migration']
        return self.migration_score(migration)

    def migration_score(self, band_migrations):
        """Score the well-numbering and well-separation of all bands.

        If some bands are too high or too low, or the number of bands is out
        of bounds, return 0. Else, return the minimal distance between two
        consecutive bands.
        """
        if not self.min_bands <= len(band_migrations) <= self.max_bands:
            return 0

        t = self.border_tolerance
        mini, maxi = (1 + t) * self.migration_min, (1 - t) * self.migration_max
        if not (mini <= min(band_migrations) <= max(band_migrations) < maxi):
            return 0

        min_gap = np.diff(sorted(band_migrations)).min()
        return 1.0 * min_gap / self.migration_span

    def _compute_elements(self):
        """Return the list of digestions to serve as elements."""
        return set(self.sequences.keys())
