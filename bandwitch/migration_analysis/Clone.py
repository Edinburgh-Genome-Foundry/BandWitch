from collections import OrderedDict
from .CloneValidation import CloneValidation


class Clone:
    """Gather all informations necessary to validate a clone.

    Parameters
    ----------
    name
      Name of the clone. Could be for instance a microplate well.

    digestions
      A dictionnary ``{digestion: CloneObservation}`` where ``digestion`` is
      of the form ``('EcoRI', 'BamHI')``.

    construct_id
      ID of the construct to be validated. This is used to group clones
      by construct in the validation reports.
    """

    def __init__(self, name, digestions, construct_id=None):
        """Initialize."""
        self.name = name
        self.digestions = digestions
        self.construct_id = construct_id

    def validate_bands(self, bands_by_digestion, relative_tolerance=0.1,
                       min_band_cutoff=None, max_band_cutoff=None):
        """Return a validation results (comparison of observed and expected).

        The result is a CloneValidation object.

        Parameters
        ----------
        bands_by_digestion
          A dictionnary ``{digestion: [bands]}`` where ``digestion`` is
          of the form ``('EcoRI', 'BamHI')``, and bands is a list of band sizes

        relative_tolerance
          Tolerance, as a ratio of the full ladder span. If =0.1, then the
          discrepancy will have a value of 1 when a band's nearest
          correspondent in the other pattern is more that 10% of the ladder
          span apart.

        min_band_cutoff
          Discrepancies involving at least one band below this minimal band
          size will be ignored. By default, it will be set to the smallest
          band size in the ladder.

        max_band_cutoff
          Discrepancies involving at least one band above this minimal band
          size will be ignored. By default, it will be set to the smallest
          band size in the ladder.

        """
        for digestion, observation in self.digestions.items():
            discrepancies = {
                digestion: observation.patterns_discrepancy(
                    other_bands=bands_by_digestion[digestion],
                    relative_tolerance=relative_tolerance,
                    min_band_cutoff=min_band_cutoff,
                    max_band_cutoff=max_band_cutoff
                )
                for digestion, observation in self.digestions.items()
            }
        return CloneValidation(self, bands_by_digestion,
                               discrepancies=discrepancies)

    @staticmethod
    def from_bands_observations(observations, constructs_map, digestions_map,
                                clones_map=None):
        if isinstance(observations, (dict, OrderedDict)):
            observations = list(observations.values())
        if clones_map is None:
            clones_map = {name: name for name in constructs_map}
        clones = OrderedDict()

        for obs in observations:
            if obs.name not in clones_map:
                continue
            clone_name = clones_map[obs.name]
            if clone_name not in clones:
                clones[clone_name] = Clone(
                    name=clone_name,
                    construct_id=constructs_map.get(clone_name, 'no_clone'),
                    digestions={}
                )
            digestion = digestions_map[obs.name]
            clones[clone_name].digestions[digestion] = obs

        return clones
