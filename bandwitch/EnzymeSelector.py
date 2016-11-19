from .tools import (greedy_minimal_set_cover, compute_sequence_digestions,
                    compute_digestion_bands, max_min_distance)
import itertools
from collections import OrderedDict
from copy import copy

import numpy as np



class EnzymeSelector:

    set_cover_finder = greedy_minimal_set_cover

    def __init__(self, ladder, relative_error=0.1):
        self.ladder = ladder
        self.relative_error = relative_error
        mini, maxi = self.ladder.migration_distances_span()
        self.migration_min, self.migration_max = mini, maxi
        self.detectable_migration_difference = relative_error * (maxi - mini)

    def bands_to_migration_pattern(self, bands):
        return self.ladder.band_size_to_migration(np.array(bands))

    def migration_patterns_look_similar(self, m1, m2):
        return max_min_distance(m1, m2) < self.detectable_migration_difference

    def digestion_separates_sequences_pair(self, enzyme, sequences_pair,
                                           sequences_digestions_dict):
        return not self.migration_patterns_look_similar(*[
            sequences_digestions_dict[sequence][enzyme]["migrations"]
            for sequence in sequences_pair
        ])

    def find_separating_digestions(self, sequences, enzymes, linear=True,
                                   max_enzymes_per_digestion=2):
        if not isinstance(sequences, OrderedDict):
            sequences = OrderedDict([
                (name, sequence)
                for name, sequence in sorted(sequences.items())
            ])
        sequences_digestions_dict = OrderedDict(
            (name, compute_sequence_digestions(
                sequence=sequence,
                enzymes=enzymes,
                linear=linear,
                max_enzymes_per_digestion=max_enzymes_per_digestion
            ))
            for name, sequence in sequences.items()
        )
        for digestions_dict in sequences_digestions_dict.values():
            for digestion in digestions_dict.values():
                migration = self.bands_to_migration_pattern(digestion["bands"])
                digestion["migrations"] = migration

        sequences_names = list(sequences.keys())
        digestions = list(sequences_digestions_dict[sequences_names[0]].keys())
        digestions_coverage = {d: [] for d in digestions}

        all_sequences_pairs = set(itertools.combinations(sequences_names, 2))
        for sequence_pair in all_sequences_pairs:
            for digestion in digestions:
                if self.digestion_separates_sequences_pair(
                        digestion, sequence_pair, sequences_digestions_dict):
                    digestions_coverage[digestion].append(sequence_pair)

        def priority_fun(enzymes):
            return -len(enzymes)
        minimal_set_cover = self.set_cover_finder(
            digestions_coverage, full_set=all_sequences_pairs,
            priority_fun=priority_fun
        )
        return minimal_set_cover, sequences_digestions_dict


    def identify_construct(self, digestions_bands, possible_constructs,
                           linear=True):
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
        possible_constructs = copy(possible_constructs)
        for digestion, bands in digestions_bands.items():
            migration = self.bands_to_migration_pattern(bands)
            for name, sequence in list(possible_constructs.items()):
                bands = compute_digestion_bands(sequence, digestion,
                                                linear=linear)
                migration2 = self.bands_to_migration_pattern(bands)

                if not self.migration_patterns_look_similar(migration,
                                                            migration2):
                    possible_constructs.pop(name)
        return sorted(possible_constructs.keys())

    def migration_pattern_is_ideal(self, migration, min_bands=2, max_bands=6):
        if (len(migration) < min_bands) or (len(migration) > max_bands):
            return False
        min_diff = np.diff(sorted(migration)).min()
        return min_diff > self.detectable_migration_difference

    def find_enzymes_for_ideal_patterns(self, enzymes, sequences_dict,
                                        linear=True,
                                        max_enzymes_per_digestion=2,
                                        **ideal_params):

        sequences_digestions_dict = OrderedDict(
            (name, compute_sequence_digestions(
                sequence=sequence,
                enzymes=enzymes,
                linear=linear,
                max_enzymes_per_digestion=max_enzymes_per_digestion
            ))
            for name, sequence in sequences_dict.items()
        )

        sequences_names = list(sequences_dict.keys())
        digestions = list(sequences_digestions_dict[sequences_names[0]].keys())
        digestions_coverage = {d: [] for d in digestions}
        for name, digestions_dict in sequences_digestions_dict.items():
            for digestion, digest in digestions_dict.items():
                migration = self.bands_to_migration_pattern(digest["bands"])
                if self.migration_pattern_is_ideal(migration, **ideal_params):
                    digestions_coverage[digestion].append(name)
        def priority_fun(enzymes):
            return -len(enzymes)
        minimal_set_cover = self.set_cover_finder(
            digestions_coverage, full_set=set(sequences_names),
            priority_fun=priority_fun
        )

        sequences_digestions = {}
        for digestion in minimal_set_cover:
            for sequence_name in digestions_coverage[digestion]:
                if sequence_name not in sequences_digestions:
                    sequences_digestions[sequence_name] = digestion

        return sequences_digestions
