import itertools
from collections import OrderedDict

import numpy as np

from .tools import (greedy_minimal_set_cover, predict_sequence_digestions,
                    predict_digestion_bands, max_min_distance, updated_dict,
                    digestions_list_to_string)
from .methylation import find_enzymes_with_no_methylation

try:
    import bandwagon
    import matplotlib.pyplot as plt
    PLOTS_AVAILABLE = True
except ImportError:
    PLOTS_AVAILABLE = False


class DigestionProblem:
    """General class for digestions problems.

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

    filter_methylated
      If True, all enzymes which may have star activity in the list of
      sequences are filtered out of the enzymes list.

    methylation
      Types of methylation to consider when ``filter_methylated`` is True.
    """

    def __init__(self, sequences, enzymes, ladder, linear=True,
                 max_enzymes_per_digestion=1, relative_error=0.1,
                 filter_methylated=True, methylation=('dam', 'dcm')):
        self.ladder = ladder
        self.sequences = sequences
        self.sequences_names = list(sequences.keys())

        if filter_methylated:
            enzymes = find_enzymes_with_no_methylation(
                enzymes,
                sequences=sequences.values(),
                methylation=methylation,
                linear=linear
            )
        self.enzymes = enzymes
        self.linear = linear
        self.max_enzymes_per_digestion = max_enzymes_per_digestion
        self.relative_error = relative_error
        mini, maxi = self.ladder.migration_distances_span
        self.migration_min, self.migration_max = mini, maxi
        self.detectable_migration_difference = relative_error * (maxi - mini)
        self.compute_coverage()

    def bands_to_migration_pattern(self, bands_sizes):
        """Return the distance migrations from several bands sizes"""
        return self.ladder.dna_size_to_migration(np.array(bands_sizes))

    def migration_patterns_look_similar(self, m1, m2):
        """Return True iff all bands migration in m1 are close from one
        band migration in m2."""
        return max_min_distance(m1, m2) < self.detectable_migration_difference

    def select_digestions(self, solver="default"):
        """Solve the digestion problem and return a list of digestion(s).

        You can provide a """
        return (self.default_solver if solver == "default" else solver)(self)

    @staticmethod
    def default_solver(problem):
        """Solver for minimum covering set problem, with penalized heuristic.

        The digestion selected at each step is the one maximizing the size
        of its coverage, minus a penalty of 1 per enzyme in the digestion,
        or 0.5 for enzymes that are already in previous digestions.
        """
        def penalized_heuristic(digestion, coverage, selected_digestions):
            selected_enzymes = [e for d in selected_digestions for e in d]
            n_covered = len(coverage[digestion])
            n_enzymes = len(digestion)
            n_common = len([e for e in digestion if e in selected_enzymes])
            return (n_covered - n_enzymes + 0.5 * n_common,
                    ' '.join(digestion))
        return greedy_minimal_set_cover(problem.coverages, problem.full_set,
                                        heuristic=penalized_heuristic)

    def compute_sequences_digestions(self):
        self.sequences_digestions = {
            name: predict_sequence_digestions(
                sequence=sequence, enzymes=self.enzymes, linear=self.linear,
                max_enzymes_per_digestion=self.max_enzymes_per_digestion,
                bands_to_migration=self.bands_to_migration_pattern)
            for name, sequence in self.sequences.items()
        }
        self.digestions = list(
            self.sequences_digestions[self.sequences_names[0]].keys())

    def digestions_set_is_covering(self, digestions):
        coverage = set().union(*[self.coverages[d] for d in digestions])
        return (coverage == self.full_set)

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

    def compute_coverage(self):
        """Compute coverage (among construct pairs) of each posisble digestion.
        """
        self.compute_sequences_digestions()
        self.full_set = set(itertools.combinations(self.sequences_names, 2))
        self.coverages = {
            digestion: self.compute_separation_coverage(digestion)
            for digestion in self.digestions
        }

    def max_patterns_difference(self, seq1, seq2, digestion):
        zone = (
            self.migration_min + 0.5 * self.detectable_migration_difference,
            self.migration_max - 0.5 * self.detectable_migration_difference
        )
        distance = max_min_distance(
            self.sequences_digestions[seq1][digestion]['migration'],
            self.sequences_digestions[seq2][digestion]['migration'],
            zone=zone)
        return 1.0 * distance

    def compute_separation_coverage(self, digestion):
        """Return all pairs of construct that are well separated
           by this digestion."""
        return [
            (seq1, seq2) for seq1, seq2 in self.full_set
            if (self.max_patterns_difference(seq1, seq2, digestion) >=
                self.detectable_migration_difference)
        ]




    @staticmethod
    def identify_constructs(observed_patterns, possible_constructs, ladder,
                            linear=False):
        """Give difference scores between observed patterns.

        The scores are from 0 to 1, 0 mening perfect match

        Parameters
        ----------

        observed_patterns
          Of the form ``{ clone_name : { digestion: [band_sizes] } }`` with all
          the observed digestion patterns for all clones.

        possible_constructs
          A dict of the form ``{ construct_name : construct_sequence_str}``
          with all candidate constructs for the observed patterns.

        ladder
          A BandWitch ladder.

        linear
          True for linear construct, else False for circular constructs.

        """
        lmini, lmaxi = ladder.migration_distances_span

        def normalized_pattern(bands):
            pattern = ladder.dna_size_to_migration(bands)
            return (1.0 * pattern - lmini) / (lmaxi - lmini)
        results = OrderedDict()
        precomputed_patterns = {
            construct: {}
            for construct in possible_constructs
        }
        for name, bands_dict in observed_patterns.items():
            results[name] = OrderedDict()
            # cst = construct
            for cst_name, cst_sequence in possible_constructs.items():
                score = -np.inf
                for digestion, bands in bands_dict.items():
                    pattern = normalized_pattern(bands)
                    if len(bands) == 0:
                        score = 1
                        continue
                    if digestion not in precomputed_patterns[cst_name]:
                        bb = predict_digestion_bands(cst_sequence,
                                                     digestion, linear=linear)
                        b_pattern = normalized_pattern(bb)
                        precomputed_patterns[cst_name][digestion] = b_pattern
                    exp_pattern = precomputed_patterns[cst_name][digestion]
                    score = max(score, max_min_distance(pattern, exp_pattern,
                                                        zone=[0.1, 0.9]))
                results[name][cst_name] = score
        return results

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

    @classmethod
    def plot_constructs_identification(cls, identification_results,
                                       observed_patterns,
                                       possible_constructs,
                                       ladder,
                                       linear=False, axes=None):
        """Plot the result of the ``.identify_construct`` method.

        Parameters
        ----------

        identification_results
          A dict ``{clone_name: {construct_name: score}}`` such as returned by
          ``SeparatingDigestionsProblem.identification_results(...)``

        observed_patterns
          Of the form ``{ clone_name : { digestion: [band_sizes] } }`` with all
          the observed digestion patterns for all clones. Same as the one
          provided to ``.identify_construct``.

        possible_constructs
          A dict of the form ``{ construct_name : construct_sequence_str}``
          with all candidate constructs for the observed patterns.
          Same as the one provided to ``.identify_construct``.

        ladder
          A BandWitch ladder.

        linear
          True for linear construct, else False for circular constructs.

        axes
          Matplotlib axes on which to plot. If non provided, new axes are
          created and returned.
        """

        nlines = len([
            digestion
            for (name, d) in observed_patterns.items()
            for digestion in d
        ])
        nconstructs = len(possible_constructs)
        fig, axes = plt.subplots(nlines, 1, sharex=True,
                                 figsize=(nconstructs, 2 * nlines))
        axes_iter = (ax for ax in axes)
        mini, maxi = min(ladder.dna_sizes), max(ladder.dna_sizes)

        for (name, scores) in identification_results.items():
            for digestion, observed_pattern in observed_patterns[name].items():
                ax = next(axes_iter)
                patterns_sets = bandwagon.BandsPatternsSet(
                    patterns=[
                        bandwagon.BandsPattern(
                            observed_pattern, ladder=ladder,
                            background_color='#c0d5f7',
                            label=("Observed" if ax == axes[0] else None))
                    ] + [
                        bandwagon.BandsPattern(
                            [
                                bandwagon.Band(
                                    b, ladder=ladder,
                                    band_color='#000000' if mini < b < maxi
                                               else '#00000044'
                                )
                                for b in predict_digestion_bands(
                                    possible_constructs[construct_name],
                                    digestion, linear=linear
                                )
                            ],
                            ladder=ladder,
                            label=construct_name if ax == axes[0] else None,
                            background_color=cls.color_function(score))
                        for construct_name, score in scores.items()

                    ],
                    label="\n".join([name, " + ".join(digestion)]),
                    ladder=ladder
                )
                patterns_sets.plot(ax)
        return axes

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
    min_bands = 3
    max_bands = 8

    def migration_pattern_is_ideal(self, migration):
        """Return True iff the pattern has the right band number and they are
           well separated."""
        if not (self.min_bands <= len(migration) <= self.max_bands):
            return False
        min_diff = np.diff(sorted(migration)).min()
        return min_diff > self.detectable_migration_difference

    def compute_coverage(self):
        self.compute_sequences_digestions()
        self.full_set = set(self.sequences_names)
        self.coverages = {
            digestion: [
                name
                for name, digestions_dict in self.sequences_digestions.items()
                if self.migration_pattern_is_ideal(
                    digestions_dict[digestion]["migration"])
            ]
            for digestion in self.digestions
        }

    def select_digestion_for_each_sequence(self, digestions):
        """Select the best digestion for each sequence."""
        for digestion in sorted(digestions, key=lambda d: len(d)):
            for sequence_name in self.digestions_coverage[digestion]:
                if sequence_name not in self.sequences_digestions:
                    self.sequences_digestions[sequence_name] = digestion

        return self.sequences_digestions
