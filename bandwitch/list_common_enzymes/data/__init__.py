"""Loads REBASE enzymes infos (mehtylation sensitivity and providers)"""

__all__ = ['enzymes_infos']

import os.path as osp

csv_path = osp.join(osp.dirname(osp.realpath(__file__)), "enzymes_infos.csv")

with open(csv_path, "r") as f:
    _lines = f.read().split("\n")
    _fields = _lines[0].split(";")
    _replacements = dict([("N/A", False), ("+", True), ("-", True)] +
                         [(str(i), i) for i in range(50)])
    enzymes_infos = {
        _line.split(";")[0]: dict(zip(_fields, [
            _replacements.get(e, e) for e in _line.split(";")]))
        for _line in _lines[1:]
    }
