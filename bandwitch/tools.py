"""Collection of useful methods and solvers for BandWitch."""

from Bio import Restriction
from Bio.Seq import Seq
import numpy as np

class NoSolutionError(Exception):
    """Specific error for the case where no enzyme set can cut it. Ah ah."""
    pass


def predict_digestion_bands(sequence, enzymes, linear=True):
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


def _compute_bands_from_cuts(cuts, sequence_length, linear=True):
    """Compute the size of the obtained bands from the position of cuts.


    Parameters
    ----------

    cuts
      Location of the different cuts on the plasmid.

    sequence_length
      Length of the DNA molecule

    linear
      True for a linear DNA molecule, False for a circular DNA molecule.
    """
    cuts = [0] + sorted(list(set(cuts))) + [sequence_length]
    bands = [cut2 - cut1 for cut1, cut2 in zip(cuts, cuts[1:])]
    if not linear and len(bands) > 1:
        bands[0] += bands.pop()
    return bands


def _merge_digestions(digestion1, digestion2, sequence_length, linear):
    """Merges and sorts the cuts from two different digestions"""
    all_cuts = sorted(list(set(digestion1["cuts"] + digestion2["cuts"])))
    return {
        "cuts": all_cuts,
        "bands": _compute_bands_from_cuts(
            cuts=all_cuts,
            sequence_length=sequence_length,
            linear=linear
        )
    }


def predict_sequence_digestions(sequence, enzymes, linear=True,
                                max_enzymes_per_digestion=1,
                                bands_to_migration=None):
    """Return a dict giving the bands sizes pattern for all possible digestions
    (i.e. all subsets of the provided enzymes).

    The result if of the form ``{digestion: {'cuts': [], 'bands': []}}``
    Where ``digestion`` is a tuple of enzyme names e.g. ``('EcoRI', 'XbaI')``,
    'cuts' is a list of cuts locations, 'bands' is a list of bands sizes

    Parameters
    ----------

    sequence
      The sequence to be digested

    enzymes
      List of all enzymes to be considered

    max_enzymes_per_digestion
      Maximum number of enzymes allowed in one digestion

    bands_to_migration
      Function associating a migration distance to a band size. If provided,
      each digestion will have a ``'migration'`` field (list of migration
      distances) in addition to 'cuts' and 'bands'.
    """
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
                digestion = tuple(sorted(enzs + (enzyme,)))
                if digestion not in digestions_dict:
                    digestions_dict[digestion] = _merge_digestions(
                        digestion1=get_cuts(enzyme),
                        digestion2=digestions_dict[enzs],
                        sequence_length=len(sequence),
                        linear=linear
                    )
                    if bands_to_migration is not None:
                        bands = digestions_dict[digestion]["bands"]
                        migration = bands_to_migration(bands)
                        digestions_dict[digestion]["migration"] = migration
    digestions_dict.pop(())
    return digestions_dict


def greedy_minimal_set_cover(coverages, full_set=None, heuristic=None):
    """Return a (hopefully) 'minimal' full-covering set, with a greedy method.

    The greedy method consists in selecting first the element with the biggest
    coverage, then the element with the biggest coverage among yet-uncovered
    targets, etc. until all targets are covered.

    Parameters
    ----------

    coverages
      A dictionary ``{element: [covered targets]}``

    full_set
      The full set of elements to be covered. Providing this full_set enables
      to quickly check that there is a solution to the problem, i.e. that
      the union of all coverages from all elements is the full set.

    heuristic
      Function ``(element) => score`` to select the element with the highest
      score at each iteration. By default, the heuristic is the length of the
      element's coverage of yet-uncovered targets.
    """
    current_selection = []
    if heuristic is None:
        def heuristic(element, coverage, current_selection):
            return len(coverage[element])

    coverages = {k: set(v) for k, v in coverages.items()}
    def key(e):
        """Function with regard to which the next best element is selected."""
        return heuristic(e, coverages, current_selection)
    if full_set is not None:
        full_coverage_set = set().union(*coverages.values())
        if full_coverage_set != full_set:
            raise NoSolutionError("Coverage not full.")
    while len(coverages) > 0:
        selected = max(coverages, key=key)
        covered = coverages.pop(selected)
        if len(covered) == 0:
            break
        current_selection.append(selected)
        for element, coverage in coverages.items():
            coverages[element] = coverage.difference(covered)
    return current_selection


def digestions_list_to_string(digestions):
    """Return a nicely formatted string of a list of 'digestion'
    (where a digestion is a list of enzyme names)"""
    return ", ".join([" + ".join([e for e in d]) for d in digestions])

def max_min_distance(values_1, values_2, zone=None):
    """Return the maximum distance between one value in one set and its
    nearest neighbor in the other set.

    A large max-min-distance means that the sets values_1 and values_2 are
    visually different, i.e. that one value in one of these sets has no close
    equivalent in the other set.
    """
    values_1, values_2 = np.array(values_1), np.array(values_2)
    all_distances = abs(values_1.reshape((len(values_1), 1)) - values_2)
    if zone is not None:
        mini, maxi = zone
        m1 = all_distances.min(axis=0)[(mini <= values_2) & (values_2 <= maxi)]
        m1 = m1.max() if len(m1) else 0
        m2 = all_distances.min(axis=1)[(mini <= values_1) & (values_1 <= maxi)]
        m2 = m2.max() if len(m2) else 0
        return max(m1, m2)

    return max([all_distances.min(axis=i).max() for i in (0, 1)])


def band_patterns_look_similar(bands1, bands2, ladder, relative_tolerance=0.1,
                               reference_and_gel=False):
    """Return True iff the observed bands for the two patterns will be similar.

    Similarity is determined by the fact that the min-max distance between
    the two migration patterns is below some threshold. The threshold is some
    proportion of the ladder's migration span.

    If mode ``reference_and_gel`` is set to True, then ``bands1`` is considered
    a "truth" (= expected bands pattern) and the function will automatically
    return False if the observed pattern on gel (bands2) has more bands than
    bands1.

    Parameters
    ----------

    bands1
      A list of bands sizes forming the first pattern

    bands2
      A list of bands sizes forming the second pattern

    ladder
      A Ladder object used to compute migration distances.

    relative tolerance
      Proportion of the ladder's migration span (between lowest and highest
      migration distances) that is the threshold for the min-max distance.

    reference_and_gel
      Set to True if ``bands1`` is a reference (= expected bands) and
      ``bands2`` is the observed bands.
    """
    if bands1 == bands2 == []:
        return True
    elif min(len(bands1), len(bands2)) == 0:
        return False
    if reference_and_gel and (len(bands2) > len(bands1)):
        return False

    m1, m2 = (ladder.dna_size_to_migration(np.array(b))
              for b in (bands1, bands2))
    mini, maxi = mini, maxi = ladder.migration_distance_span
    return max_min_distance(m1, m2) < relative_tolerance * (maxi - mini)


def updated_dict(dic1, dic2):
    """Return dic1 updated with dic2 if dic2 is not None.

    Example
    -------

    >>> my_dict = updated_dict({"size": 7, "color": "r"}, some_defaults_dict)
    """
    if dic2 is not None:
        dic1.update(dic2)
    return dic1
