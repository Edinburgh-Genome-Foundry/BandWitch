import numpy as np
from ..tools import max_min_distance


def band_patterns_discrepancy(
    bands1,
    bands2,
    ladder,
    relative_tolerance=0.1,
    reference_and_gel=False,
    zone=None,
):
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

    m1, m2 = (
        ladder.dna_size_to_migration(np.array(b)) for b in (bands1, bands2)
    )
    mini, maxi = ladder.migration_distances_span
    tolerance = relative_tolerance * (maxi - mini)
    zone = [ladder.dna_size_to_migration(b) for b in zone][::-1]
    return 1.0 * max_min_distance(m1, m2, zone=zone) / tolerance
