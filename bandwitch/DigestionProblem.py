import itertools
from copy import copy

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
        return greedy_minimal_set_cover(problem.coverages, problem.full_set)

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

    def compute_separation_coverage(self, digestion):
        """Return all pairs of construct that are well separated
           by this digestion."""
        zone = (
            self.migration_min + 0.5 * self.detectable_migration_difference,
            self.migration_max - 0.5 * self.detectable_migration_difference
        )
        return [
            (n1, n2) for n1, n2 in self.full_set
            if self.detectable_migration_difference < max_min_distance(
                self.sequences_digestions[n1][digestion]['migration'],
                self.sequences_digestions[n2][digestion]['migration'],
                zone=zone)
        ]

    def identify_construct(self, digestions_bands):
        """Identify a construct from its band profile.

        Parameters
        ----------

        digestion_bands
           {(enz1, enz2): [b1, b2, b3], (enz3,): [b4, b5, b6]}

        possible_constructs
          {name: "ATGC..."}

        linear
          Whether the constructs tested are linear
        """

        possible_sequences = copy(self.sequences)
        for digestion, bands in digestions_bands.items():
            migration = self.bands_to_migration_pattern(bands)
            for name, sequence in list(possible_sequences.items()):
                bands = predict_digestion_bands(sequence, digestion,
                                                linear=self.linear)
                migration2 = self.bands_to_migration_pattern(bands)

                if not self.migration_patterns_look_similar(migration,
                                                            migration2):
                    possible_sequences.pop(name)
        return sorted(possible_sequences.keys())

    def plot_digests(self, digestions, axes=None, bands_props=None,
                     patterns_props=None, patternset_props=None):
        """Plot a series of digests."""
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
            fig, axes = plt.subplots(ndig, 1, figsize=(0.5 * nseq, 2 * ndig))
        if ndig == 1:
            axes = [axes]
        for ax, digestion in zip(axes, digestions):
            bandwagon.BandsPatternsSet(
                patterns=[ladder] + [
                    bandwagon.BandsPattern(
                        bands=self.sequences_digestions[seq_name]
                                                       [digestion]['bands'],
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


class IdealDigestionsProblem(DigestionProblem):
    min_bands = 3
    max_bands = 8

    def migration_pattern_is_ideal(self, migration):
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
        for digestion in sorted(digestions, key=lambda d: len(d)):
            for sequence_name in self.digestions_coverage[digestion]:
                if sequence_name not in self.sequences_digestions:
                    self.sequences_digestions[sequence_name] = digestion

        return self.sequences_digestions
