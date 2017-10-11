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

    def __init__(self, bands, infos=None):
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

LADDERS = {

    '100_to_10k': Ladder(bands={
        100: 318,
        200: 287,
        300: 261,
        400: 239,
        500: 220,
        600: 203,
        700: 190,
        800: 178,
        900: 167,
        1000: 157,
        1200: 138,
        1500: 116,
        2000: 91,
        3000: 64,
        4000: 46,
        5000: 38,
        6000: 33,
        8000: 27,
        10000: 21
    }, infos='Also known as the 2-Log DNA ladder'),

    '100_to_4k': Ladder(bands={
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
    }, infos='Some typical ladder.'),

    '35_to_5k': Ladder(bands={
        5000: 234,
        3000: 407,
        2000: 528,
        1500: 607,
        1200: 666,
        1000: 733,
        900: 778,
        800: 848,
        700: 929,
        600: 1081,
        500: 1244,
        400: 1488,
        300: 1734,
        200: 2034,
        100: 2378,
        35: 2581
    }, infos="AATI Fragment analyzer ladder, from calibration file."),

    '75_to_15k': Ladder(bands={
        15000: 28,
        10000: 46,
        8000: 55,
        6000: 69,
        5000: 75,
        4000: 86,
        3500: 91,
        3000: 97,
        2500: 104,
        2000: 111,
        1500: 124,
        1000: 158,
        750: 189,
        500: 229,
        200: 276,
        75:	310
    }, infos="AATI Fragment analyzer ladder, values determined from image")
}
