from Bio import Restriction
from Bio.Seq import Seq
import numpy as np

class NoSolutionError(Exception):
    pass


def compute_digestion_bands(sequence, enzymes, linear=True):
    """Return the band sizes [75, 2101, ...] resulting from
    enzymatic digestion of the sequence by all enzymes at once.

    Parameters
    ----------

    sequence
      Sequence to be digested. Either a string ("ATGC...") or
      a BioPython `Seq`

    enzymes
      list of all enzymes placed at the same time in the digestion
      mix e.g. `["EcoRI", "BamHI"]`

    linear
      True if the DNA fragment is linearized, False if it is circular
    """
    if not isinstance(sequence, Seq):
        sequence = Seq(sequence)

    batch = Restriction.RestrictionBatch(enzymes)
    cut_sites = batch.search(sequence, linear=linear)
    cut_sites = [0] + sorted(sum(cut_sites.values(), [])) + [len(sequence)]
    bands_sizes = [end - start for (start, end)
                   in zip(cut_sites, cut_sites[1:])]
    if not linear and len(bands_sizes) > 1:
        bands_sizes[0] += bands_sizes.pop()
    return sorted(bands_sizes)


def compute_bands_from_cuts(cuts, sequence_length, linear=True):
    cuts = [0] + sorted(list(set(cuts))) + [sequence_length]
    bands = [cut2 - cut1 for cut1, cut2 in zip(cuts, cuts[1:])]
    if not linear and len(bands) > 1:
        bands[0] += bands.pop()
    return bands


def merge_digestions(digestion1, digestion2, sequence_length, linear):
    all_cuts = sorted(list(set(digestion1["cuts"] + digestion2["cuts"])))
    return {
        "cuts": all_cuts,
        "bands": compute_bands_from_cuts(
            cuts=all_cuts,
            sequence_length=sequence_length,
            linear=linear
        )
    }


def compute_sequence_digestions(sequence, enzymes, linear=True,
                                max_enzymes_per_digestion=1):
    restriction_batch = Restriction.RestrictionBatch(enzymes)
    cuts_dict = restriction_batch.search(Seq(sequence))

    def get_cuts(enzyme_name):
        return {"cuts": cuts_dict[Restriction.__dict__[enzyme_name]]}
    digestions_dict = {(): {"cuts": [], "bands": [len(sequence)]}}
    for n_enzymes in range(max_enzymes_per_digestion):
        sub_enzymes = [enzs for enzs in digestions_dict.keys()
                       if len(enzs) == n_enzymes]
        for enzyme in enzymes:
            sub_sub_enzymes = [enzs for enzs in sub_enzymes
                               if enzyme not in enzs]
            for enzs in sub_sub_enzymes:
                digestions_dict[enzs + (enzyme,)] = merge_digestions(
                    digestion1=get_cuts(enzyme),
                    digestion2=digestions_dict[enzs],
                    sequence_length=len(sequence),
                    linear=linear
                )
    return digestions_dict


def greedy_minimal_set_cover(self, coverages, full_set=None, priority_fun=None):
    """{coverer: [coverage]}"""
    result = []
    if priority_fun is None:
        def priority_fun(e):
            return 0
    priorities = {k: priority_fun(k) for k in coverages.keys()}
    def maximum_key(kv):
        return (len(kv[1]), priorities[kv[0]])
    coverages = {k: set(v) for k, v in coverages.items()}
    if full_set is not None:
        full_coverage_set = set().union(*coverages.values())
        if full_coverage_set != full_set:
            raise NoSolutionError("Coverage not full.")
    while len(coverages) > 0:
        most_common, covered = max(coverages.items(),
                                   key=maximum_key)
        if len(covered) == 0:
            break
        result.append(most_common)
        covered = coverages.pop(most_common)
        for element, coverage in coverages.items():
            coverages[element] = coverage.difference(covered)
    return result

def digestions_list_to_string(digestions):
    return ", ".join([" + ".join([e for e in d]) for d in digestions])


def max_min_distance(vals1, vals2):
    vals1, vals2 = np.array(vals1), np.array(vals2)
    all_distances = abs(vals1.reshape((len(vals1),1)) - vals2)
    return max([max(all_distances.min(axis=i)) for i in (0,1)])
