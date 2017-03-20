

from bandwitch import (SeparatingDigestionsProblem, plot_separating_digests,
                       LADDERS, digestions_list_to_string)
from dnachisel.biotools import random_dna_sequence
from gelsimulator import GelSimulator

enzymes = ["EcoRI", "BamHI", "XhoI", "EcoRV", "SpeI", "XbaI", "NotI",
           "SacI", "SmaI", "HindIII", "PstI"]
sequences = {"C%02d" % (i + 1): random_dna_sequence(4000) for i in range(8)}

problem = SeparatingDigestionsProblem(sequences, enzymes, linear=False,
                                      ladder=LADDERS['100_to_4k'],
                                      max_enzymes_per_digestion=2,
                                      relative_error=0.15)
digestions = problem.select_digestions()
print (digestions_list_to_string(digestions))
axes = plot_separating_digests(GelSimulator(LADDERS['100_to_4k']), digestions,
                               problem.sequences_digestions)
axes[0].figure.savefig("simple_combinatorial.png")
