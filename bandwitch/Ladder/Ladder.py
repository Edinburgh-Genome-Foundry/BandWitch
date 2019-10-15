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

    def __init__(self, bands, name=None, infos=None):
        self.bands = bands
        self.dna_sizes, self.migration_distances = (
            np.array(e)
            for e in zip(*sorted(bands.items()))
        )
        self._dna_size_to_migration_interpolator = CubicSpline(
            self.dna_sizes, self.migration_distances, bc_type='natural'
        )
        self.infos = infos

    def dna_size_to_migration(self, dna_sizes):
        """Return the migration distances for the given dna sizes"""
        return self._dna_size_to_migration_interpolator(dna_sizes)



    @property
    def migration_distances_span(self):
        return [
            self.dna_size_to_migration(band)
            for band in (self.dna_sizes.max(), self.dna_sizes.min())
        ]

    def __repr__(self):
        """Represent."""
        return "Ladder(%d-%d)" % (min(self.dna_sizes), max(self.dna_sizes))

    def __str__(self):
        """Represent."""
        return "Ladder(%d-%d)" % (min(self.dna_sizes), max(self.dna_sizes))
