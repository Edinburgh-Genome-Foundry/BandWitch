from gelsimulator import GelSimulator, LADDER_100_to_4k
from dnaweaver.biotools import random_dna_sequence
from bandwitch import EnzymeSelector, plot_separating_digests

enzymes = ["EcoRI", "BamHI", "XhoI", "EcoRV", "SpeI", "XbaI", "NotI"]
sequences = {"C%02d" % i: random_dna_sequence(4000) for i in range(20)}

selector = EnzymeSelector()
digestions, sequences_digestions_dict = selector.find_separating_digestions(
    sequences, enzymes, linear=False, max_enzymes_per_digestion=2)
axes = plot_separating_digests(GelSimulator(LADDER_100_to_4k), digestions,
                               sequences_digestions_dict)

axes[0].figure.savefig("simple_combinatorial.png")
