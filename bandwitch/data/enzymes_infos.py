import os

this_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)))
csv_path = os.path.join(this_dir, "enzymes_data", "enzymes_infos.csv")

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
