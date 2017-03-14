"""Minimal implementation of a Ladder, which is used to predict band migration
patterns from band DNA sizes."""

import numpy as np
from scipy.interpolate import CubicSpline

class Ladder:
    """Class to represent gel ladders. These ladders serve as a scale for
    plotting any other gel simulation.

    Parameters
    ----------

    bands
      A dictionnary of the form {dna_size: migration distance}
    """

    def __init__(self, bands):
        self.bands = bands
        self.dna_sizes, self.migration_distances = (
            np.array(e)
            for e in zip(*sorted(bands.items()))
        )
        self._dna_size_to_migration_interpolator = CubicSpline(
            self.dna_sizes, self.migration_distances, bc_type='natural'
        )

    def dna_size_to_migration(self, dna_sizes):
        """Return the migration distances for the given dna sizes"""
        return self._dna_size_to_migration_interpolator(dna_sizes)

    @property
    def migration_distances_span(self):
        return [
            self.dna_size_to_migration(band)
            for band in (self.dna_sizes.max(), self.dna_sizes.min())
        ]

LADDER_100_to_4k = Ladder(bands={
    100: 205,
    200: 186,
    300: 171,
    400: 158,
    500: 149,
    650: 139,
    850: 128,
    1000: 121,
    1650: 100,
    2000: 90,
    3000: 73,
    4000: 65
})
