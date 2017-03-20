from bandwitch import IdealDigestionsProblem, LADDERS
import os
from Bio import SeqIO
from collections import OrderedDict

# DEFINE THE CRITERIA FOR ACCEPTABLE BAND PATTERNS

class MyIdealDigestionsProblem(IdealDigestionsProblem):
    def migration_pattern_is_ideal(self, migration):
        """Are there 2-3 bands between 30% and 70% of the migration span?"""
        min_migration = 0.7 * self.migration_min + 0.3 * self.migration_max
        max_migration = 0.3 * self.migration_min + 0.7 * self.migration_max
        bands_in_central_zone = [band for band in migration
                                 if min_migration <= band <= max_migration]
        return 2 <= len(bands_in_central_zone) <= 3

# DEFINE THE SEQUENCES AND THE ENZYME SET

enzymes = ["EcoRI", "BamHI", "XhoI", "EcoRV", "SpeI", "XbaI", "NotI",
           "SacI", "SmaI", "HindIII", "PstI"]

sequences = OrderedDict([
    (f, str(SeqIO.read(os.path.join('example_data', f), 'genbank').seq))
    for f in sorted(os.listdir('example_data'))
])

# DEFINE AND SOLVE THE PROBLEM

problem = MyIdealDigestionsProblem(sequences, enzymes, linear=False,
                                  ladder=LADDERS['100_to_4k'],
                                  max_enzymes_per_digestion=2)
selected_digestions = problem.select_digestions()

# GENERATE A FIGURE OF THE BAND PATTERNS

axes = problem.plot_digestions(
    selected_digestions,
    patterns_props={'label_fontdict': {'rotation': 35}}
)
axes[0].figure.savefig("ideal_digestions.png", bbox_inches="tight")
