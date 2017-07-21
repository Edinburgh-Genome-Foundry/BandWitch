from Bio import Restriction
from Bio.Seq import Seq

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

    # TODO: support methylation ?
    # if len(methylation) > 0:
    #     inds = get_methylated_indices(sequence, methylation=methylation,
    #                                   linear=linear)

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
