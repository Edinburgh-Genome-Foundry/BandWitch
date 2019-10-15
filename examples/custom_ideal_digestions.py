from bandwitch import IdealDigestionsProblem, LADDERS, load_record
import os

# DEFINE THE CRITERIA FOR ACCEPTABLE BAND PATTERNS


class MyIdealDigestionsProblem(IdealDigestionsProblem):
    def migration_score(self, migration):
        """Are there 2-3 bands between 30% and 70% of the migration span?"""
        min_migration = 0.7 * self.migration_min + 0.3 * self.migration_max
        max_migration = 0.3 * self.migration_min + 0.7 * self.migration_max
        bands_in_central_zone = [
            band
            for band in migration
            if min_migration <= band <= max_migration
        ]
        return 2 <= len(bands_in_central_zone) <= 3


# DEFINE THE SEQUENCES AND THE ENZYME SET

enzymes = [
    "EcoRI",
    "BamHI",
    "XhoI",
    "EcoRV",
    "SpeI",
    "XbaI",
    "NotI",
    "SacI",
    "SmaI",
    "HindIII",
    "PstI",
    "NheI",
    "AfeI",
]

records_folder = os.path.join("example_data", "digestion_selection_data")
sequences = [
    load_record(os.path.join(records_folder, f), id=f)
    for f in sorted(os.listdir(records_folder))
]

# DEFINE AND SOLVE THE PROBLEM

problem = MyIdealDigestionsProblem(
    enzymes=enzymes,
    ladder=LADDERS["100_to_4k"],
    sequences=sequences,
    max_enzymes_per_digestion=2,
)
score, selected_digestions = problem.select_digestions(max_digestions=2)
print("Score:", score)

# GENERATE A FIGURE OF THE BAND PATTERNS

axes = problem.plot_digestions(
    digestions=selected_digestions,
    patterns_props={"label_fontdict": {"rotation": 35}},
)
axes[0].figure.savefig("custom_ideal_digestions.png", bbox_inches="tight")
