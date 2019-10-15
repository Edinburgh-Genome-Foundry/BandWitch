import os
from .Ladder import Ladder

this_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)))
ladders_path = os.path.join(this_dir, "data")
infos_file = "ladders_infos.csv"
with open(os.path.join(ladders_path, infos_file), "r") as f:
    infos = dict(line.split(";") for line in f.read().split("\n")[1:])


class LadderDict(dict):
    pass


LADDERS = LadderDict()
for filename in os.listdir(ladders_path):
    if filename == infos_file:
        continue
    ladder_name = filename.split(".")[0]
    with open(os.path.join(ladders_path, filename), "r") as f:
        ladder = Ladder(
            infos=infos[ladder_name],
            name=ladder_name,
            bands=dict(
                tuple(int(e) for e in line.split(";"))
                for line in f.read().split("\n")[1:]
            ),
        )
    LADDERS[ladder_name] = LADDERS.__dict__["ladder_" + ladder_name] = ladder
