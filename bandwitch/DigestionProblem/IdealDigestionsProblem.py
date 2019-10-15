from .DigestionProblem import DigestionProblem
import numpy as np


class IdealDigestionsProblem(DigestionProblem):
    """Find ideal digestion(s) to validate constructs.

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

    relative_migration_error
      Variance of the bands measured during the migration, given as a
      proportion of the total migration span (difference between the migration
      of the ladder's smallest and largest bands).

    """

    def __init__(
        self,
        enzymes,
        ladder,
        sequences,
        min_bands=3,
        max_bands=7,
        border_tolerance=0.1,
        topology="auto",
        default_topology="linear",
        max_enzymes_per_digestion=1,
        relative_migration_precision=0.1,
    ):
        """Initialize."""
        self.min_bands = min_bands
        self.max_bands = max_bands
        self.border_tolerance = border_tolerance

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

    def _parameter_element_score(self, digestion, sequence):
        """Compute the sequence's ``.migration_score`` for each digestion."""
        digestion = self.sequences_digestions[sequence][digestion]
        if digestion["same_as"] in self.scores[sequence]:
            return self.scores[sequence][digestion["same_as"]]
        migration = digestion["migration"]
        return self.migration_score(migration)

    def migration_score(self, band_migrations):
        """Score the well-numbering and well-separation of all bands.

        If some bands are too high or too low, or the number of bands is out
        of bounds, return 0. Else, return the minimal distance between two
        consecutive bands.
        """
        if not self.min_bands <= len(band_migrations) <= self.max_bands:
            return 0

        t = self.border_tolerance
        mini, maxi = (1 + t) * self.migration_min, (1 - t) * self.migration_max
        if not (mini <= min(band_migrations) <= max(band_migrations) < maxi):
            return 0

        min_gap = np.diff(sorted(band_migrations)).min()
        return 1.0 * min_gap / self.migration_span

    def _compute_elements(self):
        """Return the list of digestions to serve as elements."""
        return set(self.sequences.keys())
