"""Collection of useful methods and solvers for BandWitch."""

from Bio import Restriction
from .data.enzymes_infos import enzymes_infos
from Bio.Seq import Seq
from Bio import SeqIO
import numpy as np

class NoSolutionError(Exception):
    """Specific error for the case where no enzyme set can cut it. Ah ah."""
    pass

def list_common_enzymes(site_length=(6,), opt_temp=(37,), min_suppliers=1,
                        uniquify_sites=True,
                        avoid_methylations=('Dam', 'Dcm')):
    result = [
    str(e)
        for e in Restriction.AllEnzymes
        if str(e) in enzymes_infos
        and (len(e.site) in site_length)
        and ((opt_temp is None) or (e.opt_temp == 37))
        and enzymes_infos[str(e)]['suppliers'] >= min_suppliers
        and not any([enzymes_infos[str(e)][meth]
                     for meth in avoid_methylations])
    ]
    if uniquify_sites:
        sites_dict = {}
        for e in result:
            enzyme = Restriction.__dict__[e]
            sites_dict[enzyme.site] = e
        result = list(sites_dict.values())
    return sorted(result)

def load_genbank(filename, linear=True, name="unnamed"):
    """Load a genbank from a file."""
    record = SeqIO.read(filename, "genbank")
    record.linear = linear
    record.id = name
    record.name = name.replace(" ", "_")[:20]
    return record

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

    values_1, values2
      Two lists of values

    zone
      The max distance will be computed only on couples (v1, v2) where v1 is
      in the zone. Allows to remove border effects.
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


def band_patterns_discrepancy(bands1, bands2, ladder, relative_tolerance=0.1,
                              reference_and_gel=False, zone=None):
    """Return a discrepancy indicating whether the band patterns are
    perfectly matching (0) or dissimilar (>=1) or in-between (0-1).

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

    zone
    """
    if bands1 == bands2 == []:
        return 0
    elif min(len(bands1), len(bands2)) == 0:
        return 2.0
    if reference_and_gel and (len(bands2) > len(bands1)):
        return 2.0

    m1, m2 = (ladder.dna_size_to_migration(np.array(b))
              for b in (bands1, bands2))
    mini, maxi = ladder.migration_distances_span
    tolerance = relative_tolerance * (maxi - mini)
    zone = [ladder.dna_size_to_migration(b) for b in zone][::-1]
    return 1.0 * max_min_distance(m1, m2, zone=zone) / tolerance

def updated_dict(dic1, dic2):
    """Return dic1 updated with dic2 if dic2 is not None.

    Example
    -------

    >>> my_dict = updated_dict({"size": 7, "color": "r"}, some_defaults_dict)
    """
    if dic2 is not None:
        dic1.update(dic2)
    return dic1
