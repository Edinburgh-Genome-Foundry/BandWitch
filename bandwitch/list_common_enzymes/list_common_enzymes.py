import os
import pandas
from Bio import Restriction

this_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)))
csv_path = os.path.join(this_dir, "data", "enzymes_infos.csv")
dataframe = pandas.read_csv(csv_path, index_col="enzyme", sep=";")
enzymes_infos = dataframe.T.to_dict()


def list_common_enzymes(
    site_length=(6,),
    opt_temp=(37,),
    min_suppliers=1,
    uniquify_sites=True,
    avoided_methylations=("Dam", "Dcm"),
    accepted_methylation_effects=("-", "cut"),
):
    """List enzymes with many suppliers, no star activity, etc.

    site_length
      Tuple of the accepted lengths of enzymes sites, e.g. (6, 4) for 6bp and
      4bp.

    opt_temp
      Tuple of the accepted optimal temperatures. You would better keep it
      to a single value though, to ensure that all enzymes returned can
      co-digest.

    min_suppliers
      Minimal number of suppliers. Set to 3 or 4 to get enzymes that should
      be largely commercially available.

    uniquify_sites
      If True and two enzymes in the list have the same recognition site,
      only one will be kept. this speeds up computations in bandwitch in
      practice.

    avoided_methylations
      Tuple of the methylation types to avoid, of Dam, Dcm, CpG... see the
      next parameter for the tolerance level

    accepted_methylation_effects
      Accepted effects of the methylations specified by "avoided_methylations".
      This data comes from the REBASE online database.
      The effects are among "-" (no overlap, no effect), "cut" (overlap with
      a methylation site, but the enzyme still cuts), "variable" (whatever
      that means), "impaired", "some blocked", "some impaired", and "blocked".

    """
    result = [
        str(e)
        for e in Restriction.AllEnzymes
        if str(e) in enzymes_infos
        and (len(e.site) in site_length)
        and ((opt_temp is None) or (e.opt_temp == 37))
        and enzymes_infos[str(e)]["suppliers"] >= min_suppliers
        and not any(
            [
                enzymes_infos[str(e)][meth] not in accepted_methylation_effects
                for meth in avoided_methylations
            ]
        )
    ]
    if uniquify_sites:
        sites_dict = {}
        for e in result:
            enzyme = Restriction.__dict__[e]
            sites_dict[enzyme.site] = e
        result = list(sites_dict.values())
    return sorted(result)
