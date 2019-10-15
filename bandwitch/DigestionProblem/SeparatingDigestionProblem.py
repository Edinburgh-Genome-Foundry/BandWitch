from collections import OrderedDict
import itertools

import numpy as np

from .DigestionProblem import DigestionProblem
from ..tools import max_min_distance

try:
    import matplotlib.pyplot as plt

    PLOTS_AVAILABLE = True
except ImportError:
    PLOTS_AVAILABLE = False


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

    def __init__(
        self,
        enzymes,
        ladder,
        sequences=None,
        categories=None,
        topology="auto",
        default_topology="linear",
        max_enzymes_per_digestion=1,
        min_discrepancy="auto",
        relative_migration_precision=0.1,
    ):
        """Initialize."""
        if isinstance(sequences, (list, tuple)):
            sequences = OrderedDict([(r.id, r) for r in sequences])

        if categories is None:
            categories = OrderedDict(
                [
                    (seq_name, {seq_name: sequence})
                    for seq_name, sequence in sequences.items()
                ]
            )
        else:
            sequences = OrderedDict(
                [
                    (seq_name, sequence)
                    for category, seqs in categories.items()
                    for seq_name, sequence in seqs.items()
                ]
            )
        self.categories = categories

        DigestionProblem.__init__(
            self,
            sequences=sequences,
            enzymes=enzymes,
            ladder=ladder,
            topology=topology,
            default_topology=default_topology,
            max_enzymes_per_digestion=max_enzymes_per_digestion,
            relative_migration_precision=relative_migration_precision,
        )

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
        digestion1, digestion2 = [
            self.sequences_digestions[s][digestion]
            for s in (sequence1, sequence2)
        ]
        # If similar pair already computed, return the previous result.
        if digestion1["same_as"] == digestion2["same_as"]:
            if digestion1["same_as"] in self.scores[sequences_pair]:
                return self.scores[sequences_pair][digestion1["same_as"]]

        mini, maxi = self.migration_min, self.migration_max
        margin = self.relative_migration_precision * self.migration_span / 2.0
        distance = max_min_distance(
            digestion1["migration"],
            digestion2["migration"],
            zone=[mini + margin, maxi - margin],
        )
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
        return (
            max(0, min(1, score / maxi)),
            min(1, max(0, 1 - score / maxi)),
            0,
            0.5,
        )

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
        ax.imshow(
            grid[:, ::-1],
            interpolation="nearest",
            cmap="OrRd_r",
            vmin=0,
            vmax=self.relative_migration_precision,
        )
        for i in range(len(grid)):
            for j in range(len(grid)):
                if i > j:
                    ax.text(
                        len(self.sequences) - i - 1,
                        j,
                        "%d%%" % (100 * grid[j, i]),
                        fontdict=dict(color="black", weight="bold", size=14),
                        horizontalalignment="center",
                        verticalalignment="center",
                    )

        ax.set_yticks(range(len(grid)))
        ax.set_yticklabels(
            list(self.sequences)[:-1], size=14, fontdict={"weight": "bold"}
        )
        ax.set_xticks(range(len(grid)))
        ax.xaxis.set_ticks_position("top")
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.set_xticklabels(
            [" " + s for s in list(self.sequences)[1:][::-1]],
            rotation=90,
            size=14,
            fontdict={"weight": "bold"},
        )
        ax.set_xlim(-0.5, len(self.sequences) - 1.5)
        ax.set_ylim(len(self.sequences) - 1.5, -0.5)

        if target_file is not None:
            ax.figure.savefig(target_file, bbox_inches="tight")
            plt.close(ax.figure)
        return ax
