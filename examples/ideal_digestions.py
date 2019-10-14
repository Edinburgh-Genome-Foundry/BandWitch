from bandwitch import IdealDigestionsProblem, LADDERS, load_record
import os

# DEFINE THE SEQUENCES AND THE ENZYME SET

enzymes = ["EcoRI", "BamHI", "XhoI", "EcoRV", "SpeI", "XbaI",
           "NotI", "SacI", "SmaI", "HindIII", "PstI"]

records_folder = os.path.join('example_data', 'digestion_selection_data')
sequences = [
    load_record(os.path.join(records_folder, f), id=f)
    for f in sorted(os.listdir(records_folder))
]

# DEFINE AND SOLVE THE PROBLEM
problem = IdealDigestionsProblem(enzymes=enzymes,
                                 ladder=LADDERS['100_to_4k'],
                                 sequences=sequences,
                                 max_enzymes_per_digestion=2)

# Here we ask for one digestion
score, selected_digestions = problem.select_digestions(max_digestions=1)

# GENERATE A FIGURE OF THE BAND PATTERNS

problem.plot_digestions(
  digestions=selected_digestions,
  patterns_props={'label_fontdict': {'rotation': 35}},
  target_file="ideal_digestions.png"
)
