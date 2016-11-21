from .tools import (greedy_minimal_set_cover, predict_sequence_digestions,
                    predict_digestion_bands, max_min_distance)
import itertools
from collections import OrderedDict
from copy import copy

import numpy as np

class DigestionProblem:

    def __init__(self, sequences, enzymes, ladder, linear=True,
                 max_enzymes_per_digestion=1, relative_error=0.1):
        self.ladder = ladder
        self.sequences = sequences
        self.sequences_names = list(sequences.keys())
        self.enzymes = enzymes
        self.linear = linear
        self.max_enzymes_per_digestion = max_enzymes_per_digestion
        self.relative_error = relative_error
        mini, maxi = self.ladder.migration_distances_span()
        self.migration_min, self.migration_max = mini, maxi
        self.detectable_migration_difference = relative_error * (maxi - mini)

        self.compute_coverage()

    def bands_to_migration_pattern(self, bands):
        return self.ladder.band_size_to_migration(np.array(bands))

    def migration_patterns_look_similar(self, m1, m2):
        return max_min_distance(m1, m2) < self.detectable_migration_difference

    def select_digestions(self, solver="default"):
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

    def compute_coverage(self):
        self.compute_sequences_digestions()
        self.full_set = set(itertools.combinations(self.sequences_names, 2))
        self.coverages = {
            digestion: self.compute_separation_coverage(digestion)
            for digestion in self.digestions
        }

    def compute_separation_coverage(self, digestion):
        mask_size = self.detectable_migration_difference
        mask_bins = 5  # <= Hey ! a magic number !
        bin_size = mask_size / mask_bins
        n_bins = int(1.5 * self.migration_max / bin_size) + 1
        n_seqs = len(self.sequences)
        patterns = np.zeros((n_seqs, n_bins), dtype="uint8")
        pattern_masks = np.zeros((n_seqs, n_bins), dtype="uint8")

        kernel = np.ones(2*mask_bins+1, dtype="uint8")
        indices = {}
        items = self.sequences_digestions.items()
        for i, (name, digestions_dict) in enumerate(items):
            indices[name] = i
            migration = digestions_dict[digestion]["migration"]
            bins_indices = np.round(migration / bin_size).astype(int)
            patterns[i, bins_indices] = 1
            band_zones = np.convolve(patterns[i], kernel, mode="same")
            pattern_masks[i] = np.logical_not(band_zones)
        patterns_differ = np.zeros((n_seqs, n_seqs))
        for i, pattern in enumerate(patterns):
            bands_outside_masks = np.logical_and(pattern_masks, pattern)
            patterns_differ[i] = np.any(bands_outside_masks, axis=1)

        return [
            (n1, n2) for n1, n2 in self.full_set
            if (patterns_differ[indices[n1], indices[n2]] or
                patterns_differ[indices[n2], indices[n1]])
        ]

    def identify_construct(self, digestions_bands):
        """
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
            for sequence_name in digestions_coverage[digestion]:
                if sequence_name not in sequences_digestions:
                    sequences_digestions[sequence_name] = digestion

        return sequences_digestions
