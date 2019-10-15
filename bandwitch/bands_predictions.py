"""Module for digestion pattern prediction."""

from Bio import Restriction
from Bio.Seq import Seq
from collections import OrderedDict
import numpy as np

from .tools import all_subsets


def _compute_bands_from_cuts(cuts, sequence_length, linear=True):
    """Compute the size of the obtained bands from the position of cuts.

    Returns a list of band sizes.

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


def find_cuts(sequence, enzymes, linear=True):
    batch = Restriction.RestrictionBatch(enzymes)
    return [
        cut
        for cuts in batch.search(sequence, linear=linear).values()
        for cut in cuts
    ]


def predict_digestion_bands(
    sequence, enzymes, linear=True, partial_cutters=()
):
    """Return the band sizes from digestion by all enzymes at once.

    Returns a list of bands sizes sorted from smallest to largest

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

    partial_cutters
      List of enzymes that are tired or inhibited would randomly miss
      some sites, thus creating extra bands in the final digest.

    """
    if not isinstance(sequence, Seq):
        sequence = Seq(sequence)
    cuts = find_cuts(sequence, enzymes, linear=linear)
    if len(partial_cutters) == 0:
        bands = _compute_bands_from_cuts(cuts, len(sequence), linear=linear)
    else:
        propartial_cuts = find_cuts(sequence, partial_cutters, linear=linear)
        bands = set(
            band
            for partial_cuts in all_subsets(propartial_cuts)
            for band in _compute_bands_from_cuts(
                cuts + list(partial_cuts), len(sequence), linear
            )
        )
    return sorted(bands)


def compute_sequence_digestions_migrations(sequences_digestions, ladder):
    """Add a 'migration' field to the data in sequences_digestions.

    ``sequences_digestions`` is a dict ``{seq: digestions_data}``
    where ``digestions_data`` is a result of ``predict_sequence_digestions``,
    ``ladder`` is a Ladder object.
    """
    all_bands = set()
    for seq, digestions_data in sequences_digestions.items():
        for digestion, data in digestions_data.items():
            all_bands.update(data["bands"])
    all_bands = sorted(all_bands)
    all_migrations = ladder.dna_size_to_migration(list(all_bands))
    band_to_migration = {
        band: migration for band, migration in zip(all_bands, all_migrations)
    }
    for seq, digestions_data in sequences_digestions.items():
        for digestion, data in digestions_data.items():
            data["migration"] = np.array(
                [band_to_migration[b] for b in data["bands"]]
            )


def predict_sequence_digestions(
    sequence, enzymes, linear=True, max_enzymes_per_digestion=1
):
    """Return a dict giving bands sizes pattern for all possible digestions.

    The digestions, double-digestions, etc. are listed and for each the
    sequence band sizes are computed.

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

    def _merge_digestions(digestion1, digestion2, sequence_length, linear):
        """Merge and sort the cuts from two different digestions."""
        all_cuts = sorted(list(set(digestion1["cuts"] + digestion2["cuts"])))
        return {
            "cuts": all_cuts,
            "bands": _compute_bands_from_cuts(
                cuts=all_cuts, sequence_length=sequence_length, linear=linear
            ),
        }

    empty_digestion = ((), {"cuts": [], "bands": [len(sequence)]})
    digestions_dict = OrderedDict([empty_digestion])
    for n_enzymes in range(max_enzymes_per_digestion):
        sub_enzymes = [
            enzs for enzs in digestions_dict.keys() if len(enzs) == n_enzymes
        ]
        for enzyme in enzymes:
            sub_sub_enzymes = [
                enzs for enzs in sub_enzymes if enzyme not in enzs
            ]
            for enzs in sub_sub_enzymes:
                digestion = tuple(sorted(enzs + (enzyme,)))
                if digestion not in digestions_dict:
                    no_enzyme_band = len(get_cuts(enzyme)["cuts"]) == 0
                    no_enzs_band = len(digestions_dict[enzs]["cuts"]) == 0
                    one_no_bands = no_enzs_band or no_enzyme_band
                    if ((enzyme,) in digestions_dict) and one_no_bands:
                        if no_enzyme_band:
                            digestions_dict[digestion] = digestions_dict[enzs]
                            digestions_dict[digestion]["same_as"] = enzs
                        elif no_enzs_band:
                            dig = (enzyme,)
                            digestions_dict[digestion] = digestions_dict[dig]
                            digestions_dict[digestion]["same_as"] = dig
                    else:
                        digestions_dict[digestion] = _merge_digestions(
                            digestion1=get_cuts(enzyme),
                            digestion2=digestions_dict[enzs],
                            sequence_length=len(sequence),
                            linear=linear,
                        )
                        digestions_dict[digestion]["same_as"] = digestion
    digestions_dict.pop(())

    # Reordering the dictionnary makes the computation of scores faster
    # using 'same_as' to avoid recomputing scores involving similar patterns
    digestions_dict = OrderedDict(
        sorted(
            digestions_dict.items(),
            key=lambda item: (len(item[0]), len(item[1]["cuts"])),
        )
    )
    return digestions_dict
