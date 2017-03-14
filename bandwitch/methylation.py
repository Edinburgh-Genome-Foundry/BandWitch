"""Methods to deal with the methylation (=protection) of some sites."""

import re
from Bio import Restriction
from Bio.Seq import Seq
import numpy as np
from tqdm import tqdm

def get_methylated_indices(sequence, methylation=('dam', 'dcm'), linear=True):
    """Return the list of indices of the centers of methylation sites.

    sequence
      A string of the DNA sequence to analyze

    methylation
      The methylation type. Either ``'dam'`` or ``'dcm'`` or both:
      ``(dam, dcm)``

    linear
      True for linear DNA molecules, False for circular.
    """
    if isinstance(methylation, (list, tuple)):
        return sorted([
            cut
            for meth in methylation
            for cut in get_methylated_indices(sequence, meth, linear)
        ])
    site = {"dam": "GATC", "dcm": "CC[AT]GG"}[methylation]
    if not linear:
        used_sequence = sequence + sequence[:len(site)]
    else:
        used_sequence = sequence
    return [
        (match.start() + len(site)/2) % len(sequence)
        for match in re.finditer(site, used_sequence)
    ]

def sites_are_spread_apart(sites1, sites2, min_distance, linear=True,
                           sequence_size=None):
    """Return whether all sites locations in ``sites1`` are at least at
    ``min_distance`` from all sites in ``sites2``.

    Parameters
    ----------

    (sites1, sites2)
      Lists of integers representing locations in a DNA sequence.

    min_distance
      True for linear DNA molecules, False for circular.

    linear
      True for linear DNA molecules, False for circular.

    sequence_size
      Size of the DNA sequence in which these sites are. Only required for
      circular sequences (i.e. linear = False)
    """
    diffs = abs(np.array(sites1) - np.array([sites2]).T)
    if linear:
        return diffs.min() > min_distance
    else:
        return min(diffs.min(), (sequence_size - diffs).min()) > min_distance

def find_enzymes_with_no_methylation(enzymes, sequences, min_distance=5,
                                     methylation=('dam', 'dcm'),
                                     linear=True, progress_bar=False):
    """Find which enzymes have no methylation problem in provided sequences.

    At the moment a bit naive and inexact, we leave a big min_distance between
    the dam sites (GATC) and enzyme sites, as a conservative measure to be
    sure.

    Parameters
    ----------

    enzymes
      List of enzyme names. They have to be in the Biopython.Restriction list.

    sequences
      List of ATGC strings of sequences in which we want to avoid methylation
      sites.

    min_distance
      minimal Distance between a methylation site and an enzyme restriction
      site for the enzyme to be considered "methylation-problem-free".

    methylation
      The methylation type. Either ``'dam'`` or ``'dcm'`` or both:
      ``(dam, dcm)``

    linear
      True for linear DNA molecules, False for circular.

    progress_bar
      True for a progress bar as the sequences are evaluated.

    """
    enzymes = set(enzymes)
    if progress_bar:
        sequences = tqdm(sequences)
    for sequence in sequences:
        if len(enzymes) == 0:
            continue
        methylated_sites = get_methylated_indices(
            sequence, methylation, linear=linear)
        if methylated_sites == []:
            continue
        batch = Restriction.RestrictionBatch(enzymes)
        cuts_dict = batch.search(Seq(sequence))
        for enzyme, cuts in cuts_dict.items():
            if cuts == []:
                continue
            if not sites_are_spread_apart(methylated_sites, cuts, min_distance,
                                          linear=linear,
                                          sequence_size=len(sequence)):
                if enzyme in enzymes:
                    enzymes.remove(enzyme)
                else:
                    enzymes.remove(str(enzyme))
    return enzymes
