"""Collection of useful methods and solvers for BandWitch."""
import itertools
import os

import numpy as np
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import DNAAlphabet
from snapgene_reader import snapgene_file_to_seqrecord


def set_record_topology(record, topology, pass_if_already_set=False):
    """Set record.annotations['topology'], optionally passing if already set.
    """
    record_topology = record.annotations.get("topology", None)
    do_nothing = pass_if_already_set and (record_topology is not None)
    if not do_nothing:
        record.annotations["topology"] = topology


def record_is_linear(record, default=True):
    """Return true if record.annotations['topology'] == 'linear'"""
    if "topology" not in record.annotations:
        return default
    else:
        return record.annotations["topology"] == "linear"


def load_record(
    record_file,
    topology="auto",
    default_topology="linear",
    id="auto",
    upperize=True,
    max_name_length=20,
    file_format=None,
):
    """Read a record (from many different input formats).

    Parameters
    ----------

    record_file
      A genbank file, a fasta file, a snapgene file, or a filelike object
      (at which case the format, genbank or fasta, must be given with
      ``file_format``)

    topology
      Either circular or linear or auto. If auto, then will attempt to read
      record.annotations['topology], and default to ``default_topology``.

    default_topology
      Default topology to use when topology is set to "auto" and a record
      has no designated topology.

    id
      Will be used for the record ID and name. If auto, the record id will
      be unchanged unless it is ".", " ", etc. at which case it will be
      replaced by the file name.

    upperize
      If true, the sequence will get upperized.

    max_name_length
      The name of the record will be truncated if too long to avoid Biopython
      exceptions being raised.

    file_format
      Indicates the file format for the parser, when record_file is a filelike
      object.

    """
    if file_format is not None:
        record = SeqIO.read(record_file, file_format)
    elif record_file.lower().endswith(("gb", "gbk")):
        record = SeqIO.read(record_file, "genbank")
    elif record_file.lower().endswith(("fa", "fasta")):
        record = SeqIO.read(record_file, "fasta")
    elif record_file.lower().endswith(".dna"):
        record = snapgene_file_to_seqrecord(record_file)
    else:
        raise ValueError("Unknown format for file: %s" % record_file)
    if upperize:
        record = record.upper()
    if topology == "auto":
        set_record_topology(record, default_topology, pass_if_already_set=True)
    else:
        set_record_topology(record, topology)
    if id == "auto":
        id = record.id
        if id in [None, "", "<unknown id>", ".", " "]:
            id = os.path.splitext(os.path.basename(record_file))[0]
            record.name = id.replace(" ", "_")[:max_name_length]
        record.id = id
    elif id is not None:
        record.id = id
        record.name = id.replace(" ", "_")[:max_name_length]

    return record


def sequence_to_biopython_record(
    sequence, id="<unknown id>", name="<unknown name>", features=()
):
    """Return a SeqRecord of the sequence, ready to be Genbanked."""
    return SeqRecord(
        Seq(sequence, alphabet=DNAAlphabet()),
        id=id,
        name=name,
        features=list(features),
    )


def updated_dict(dic1, dic2):
    """Return dic1 updated with dic2 if dic2 is not None.

    Example
    -------

    >>> my_dict = updated_dict({"size": 7, "color": "r"}, some_defaults_dict)
    """
    if dic2 is not None:
        dic1.update(dic2)
    return dic1


def all_subsets(mylist):
    return itertools.chain.from_iterable(
        itertools.combinations(mylist, k) for k in range(len(mylist) + 1)
    )


def max_min_distance(values_1, values_2, zone=None):
    """Return the maximum distance between one value in one set and its
    nearest neighbor in the other set.

    A large max-min-distance means that the sets values_1 and values_2 are
    visually different, i.e. that one value in one of these sets has no close
    equivalent in the other set.

    values_1, values2
      Two lists of values

    zone
      The max distance will be computed only on couples (v1, v2) where v1 is
      in the zone. Allows to remove border effects.
    """
    values_1, values_2 = np.array(values_1), np.array(values_2)
    all_distances = abs(values_1.reshape((len(values_1), 1)) - values_2)
    if zone is not None:
        mini, maxi = zone
        m1 = all_distances.min(axis=0)[(mini <= values_2) & (values_2 <= maxi)]
        m1 = m1.max() if len(m1) else 0
        m2 = all_distances.min(axis=1)[(mini <= values_1) & (values_1 <= maxi)]
        m2 = m2.max() if len(m2) else 0
        return max(m1, m2)

    return max([all_distances.min(axis=i).max() for i in (0, 1)])
