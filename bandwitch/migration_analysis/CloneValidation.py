
class CloneValidation:
    """Report from comparing a clone's pband patterns with expectations.

    Instances of this class are produced by validating a ``Clone`` object.

    Parameters
    ----------
    clone
      A Clone object.

    expected
      A dictionnary ``{digestion: expected_pattern}`` where ``digestion`` is
      of the form ``('EcoRI', 'BamHI')`` and ``expected_pattern`` is a
      BandsPattern

    discrepancies
      A dictionnary ``{digestion: float}`` where the float indicates the
      observed discrepancy between observed and predicted for this digestion.

    """

    def __init__(self, clone, expected, discrepancies):
        """Initialize."""
        self.clone = clone
        self.expected = expected
        self.discrepancies = discrepancies

    @property
    def passes(self):
        """True iff the clone is valid.

        (This means that all patterns are similar enough to expectations).
        """
        return self.max_discrepancy < 1.0

    @property
    def max_discrepancy(self):
        """Return the largest discrepancy in all patterns"""
        return max(self.discrepancies.values())

    def color(self, discrepancy='max'):
        """Return a color depending on the discrepancy level (in 0-1)."""
        if discrepancy == "max":
            discrepancy = self.max_discrepancy
        factor = min(1.0, discrepancy ** 10)
        return ((2 + factor) / 3.0, (3 - factor) / 3.0, 2 / 3.0)

    def to_bandwagon_bandpattern(self, digestion, label='auto',
                                 per_digestion_discrepancy=True):
        """Return a (plottable) bandwagon pattern for the selected digestion.

        The pattern has a color indicating the match similarity, as well as
        a corner note indicating the gap to similarity.

        Parameters
        ----------
        digestion
          A set of enzymes such as ``('EcoRI',)`` or ````

        label
          A label for the pattern, for instance the name of a microplate well.

        per_digestion_discrepancy
          If true, each pattern for each clone will have an indication of
          how close it is from the expected pattern. If False, and the
          validation involves several digestions, the biggest discrepancy
          among all patterns is shown.

        """
        if digestion not in self.clone.digestions:
            return None
        if per_digestion_discrepancy:
            discrepancy = self.discrepancies[digestion]
        else:
            discrepancy = self.max_discrepancy
        obs = self.clone.digestions[digestion]
        pattern = obs.to_bandwagon_bandpattern(
            background_color=self.color(discrepancy=discrepancy), label=label)
        pattern.corner_note += " Gap: %d" % (100 * discrepancy)
        return pattern
