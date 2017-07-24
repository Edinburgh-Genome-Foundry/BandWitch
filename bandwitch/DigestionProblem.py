import itertools
from collections import OrderedDict

import numpy as np

from .tools import (max_min_distance, updated_dict,
                    digestions_list_to_string)
from .SetCoverProblem import SetCoverProblem
from .bands_predictions import predict_sequence_digestions
from .methylation import find_enzymes_with_no_methylation

try:
    import bandwagon
    import matplotlib
    import matplotlib.pyplot as plt
    PLOTS_AVAILABLE = True
except ImportError:
    PLOTS_AVAILABLE = False


class DigestionProblem(SetCoverProblem):
    """General class for solving digestions problems.

    By digestion problem we mean *enzymatic* digestion problems. For other
    kinds of digestion problems, consult your proctologist.

    Other digestion problems subclass this problem and implement a computation
    of coverages which depends on the problem.

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

    methylation
      Types of methylation to consider. Set to () for no such considerations
    """

    def __init__(self, enzymes, ladder, sequences, linear=True,
                 max_enzymes_per_digestion=1,
                 relative_migration_precision=0.1,
                 methylation=('dam', 'dcm'), progress_logger=None):
        self.sequences = sequences
        self.ladder = ladder
        self.linear = linear
        self.relative_migration_precision = relative_migration_precision
        self.max_enzymes_per_digestion = max_enzymes_per_digestion

        self.sequences_names = list(sequences.keys())
        mini, maxi = self.ladder.migration_distances_span
        self.migration_min, self.migration_max = mini, maxi
        self.migration_span = maxi - mini

        if len(methylation) > 0:
            enzymes = find_enzymes_with_no_methylation(
                enzymes,
                sequences=sequences.values(),
                methylation=methylation,
                linear=linear
            )
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
        elements = self.compute_elements()
        SetCoverProblem.__init__(self, elements=elements, parameters=parameters,
                                 progress_logger=progress_logger)

    @staticmethod
    def default_heuristic(named_subset, selected):
        enzymes, subset = named_subset
        return len(subset) - 0.5 * len(enzymes)

    def bands_to_migration_pattern(self, bands_sizes):
        """Return the distance migrations from several bands sizes"""
        return self.ladder.dna_size_to_migration(np.array(bands_sizes))

    def plot_digestions(self, digestions, axes=None, bands_props=None,
                        patterns_props=None, patternset_props=None):
        """Plot the patterns for each sequence, for each digestion in the list.

        Requires Bandwagon.

        Parameters
        ----------

        digestions, axes=None, bands_props=None,
                            patterns_props=None, patternset_props=None
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
        return axes


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

    methylation
      Types of methylation to consider.
    """

    def __init__(self, enzymes, ladder, sequences=None, categories=None,
                 linear=True, max_enzymes_per_digestion=1,
                 min_discrepancy='auto', methylation=('dam', 'dcm'),
                 relative_migration_precision=0.1):

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
            linear=True, max_enzymes_per_digestion=max_enzymes_per_digestion,
            relative_migration_precision=relative_migration_precision,
            methylation=methylation)

    def compute_elements(self):
        category_pairs = itertools.combinations(self.categories.values(), 2)
        return set(
            (seq1, seq2)
            for category1, category2 in category_pairs
            for seq1 in category1
            for seq2 in category2
        )

    def parameter_element_score(self, digestion, sequences_pair):
        """See max_patterns_difference"""
        sequence1, sequence2 = sequences_pair
        mini, maxi = self.migration_min, self.migration_max
        span = self.migration_span
        zone = (mini + self.relative_migration_precision * span / 2.0,
                maxi - self.relative_migration_precision * span / 2.0)
        distance = max_min_distance(
            self.sequences_digestions[sequence1][digestion]['migration'],
            self.sequences_digestions[sequence2][digestion]['migration'],
            zone=zone)
        return 1.0 * distance / span

    @staticmethod
    def score_to_color(score, maxi=0.1):
        """Transform a similarity score to a green/red color

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

    def plot_distances_map(self, digestions, ax=None):
        """Make a plot of how well the digestions separate each construct pair

        Parameters
        ----------

        digestions
          A list of digestions, eg ``[('EcoRV'), ('XbaI', 'MfeI')]``.

        ax
          A matplotlib ax on which to plot, if none is provided, one is created
          and returned at the end.

        """

        if not PLOTS_AVAILABLE:
            raise ImportError("Plots require Matplotlib/Bandwagon installed.")
        grid = np.zeros(2 * (len(self.sequences),))
        for i, seq1 in enumerate(self.sequences):
            for j, seq2 in enumerate(self.sequences):
                if i >= j:
                    grid[i, j] = np.nan
                else:
                    scores = [
                        self.max_patterns_difference(seq1, seq2, digestion) /
                        (self.migration_max - self.migration_min)
                        for digestion in digestions
                    ]
                    grid[i, j] = max(scores)
        if ax is None:
            _, ax = plt.subplots(1, figsize=2 * (0.8 * len(grid),))
        ax.imshow(grid[:, ::-1], interpolation='nearest', cmap='OrRd_r',
                  vmin=0, vmax=0.5)
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
        return ax


class IdealDigestionsProblem(DigestionProblem):

    def parameter_element_score(self, digestion, sequence):
        """Return True iff the pattern has the right band number and they are
           well separated."""
        migration = self.sequences_digestions[sequence][digestion]['migration']
        return self.migration_score(migration)

    def migration_score(self, migration):
        if not 1 < len(migration) < 8:
            return 0
        score_n_bands = 0.5 / (1 + (len(migration) - 4)**2)
        min_gap = np.diff(sorted(migration)).min()
        score_separation = min_gap / self.migration_span
        return score_n_bands + score_separation

    def compute_elements(self):
        return set(self.sequences.keys())
