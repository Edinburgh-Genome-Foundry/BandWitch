

from bandwitch import (EnzymeSelector, plot_separating_digests,
                       LADDER_100_to_4k, digestions_list_to_string)
from dnaweaver.biotools import random_dna_sequence
from gelsimulator import GelSimulator

enzymes = ["EcoRI", "BamHI", "XhoI", "EcoRV", "SpeI", "XbaI", "NotI",
           "SacI", "SmaI", "HindIII", "PstI"]
sequences = {"C%02d" % (i+1): random_dna_sequence(4000) for i in range(50)}

selector = EnzymeSelector(LADDER_100_to_4k, relative_error=0.15)
digestions, sequences_digestions_dict = selector.find_separating_digestions(
    sequences, enzymes, linear=False, max_enzymes_per_digestion=3)
print (digestions_list_to_string(digestions))
axes = plot_separating_digests(GelSimulator(LADDER_100_to_4k), digestions,
                               sequences_digestions_dict)
axes[0].figure.savefig("simple_combinatorial.png")
