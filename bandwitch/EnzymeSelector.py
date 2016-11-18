from .tools import greedy_minimal_set_cover, compute_sequence_digestions
import itertools
from collections import OrderedDict


class EnzymeSelector:

    set_cover_finder = greedy_minimal_set_cover

    def __init__(self, relative_precision=0.2, ladder_min=100,
                 ladder_max=4000):
        self.relative_precision = relative_precision
        self.ladder_min = ladder_min
        self.ladder_max = ladder_max

    def band_is_in_ladder(self, band):
        return self.ladder_min < band < self.ladder_max

    def bands_look_similar(self, b1, b2):
        if self.band_is_in_ladder(b1) != self.band_is_in_ladder(b2):
            return False
        difference = float(abs(b1 - b2)) / min(b1, b2)
        return difference < self.relative_precision

    def merge_undistinguishable_bands_in_pattern(self, bands):
        if bands == []:
            return []
        bands = sorted(bands)
        new_bands = [bands[0]]
        for band in bands[1:]:
            if self.bands_look_similar(band, new_bands[-1]):
                new_bands[-1] = 0.5 * (band + new_bands[-1])
            else:
                new_bands.append(band)
        return new_bands

    def bands_pattern_to_digest_observation(self, bands):
        bands = [b for b in bands if self.band_is_in_ladder(b)]
        return self.merge_undistinguishable_bands_in_pattern(bands)

    def bands_patterns_look_similar(self, bands1, bands2,
                                    process_bands=True):
        if process_bands:
            bands1, bands2 = [
                self.bands_pattern_to_digest_observation(bands)
                for bands in (bands1, bands2)
            ]
        return (len(bands1) == len(bands2)) and all(
            self.bands_look_similar(b1, b2)
            for b1, b2 in zip(sorted(bands1), sorted(bands2))
        )

    def digestion_separates_sequences_pair(self, enzyme, sequences_pair,
                                           sequences_digestions_dict):
        p1, p2 = [
            sequences_digestions_dict[sequence][enzyme]["observed_bands"]
            for sequence in sequences_pair
        ]
        return not self.bands_patterns_look_similar(p1, p2)

    def find_separating_digestions(self, sequences_dict, enzymes,
                                   max_enzymes_per_digestion=2,
                                   linear=True):
        if not isinstance(sequences_dict, OrderedDict):
            sequences_dict = OrderedDict([
                (name, sequence)
                for name, sequence in sorted(sequences_dict.items())
            ])
        sequences_digestions_dict = OrderedDict(
            (name, compute_sequence_digestions(
                sequence=sequence,
                enzymes=enzymes,
                linear=linear,
                max_enzymes_per_digestion=max_enzymes_per_digestion
            ))
            for name, sequence in sequences_dict.items()
        )
        for digestions_dict in sequences_digestions_dict.values():
            for digestion in digestions_dict.values():
                bands = self.bands_pattern_to_digest_observation(digestion[
                                                                 "bands"])
                digestion["observed_bands"] = bands

        sequences_names = list(sequences_dict.keys())
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
