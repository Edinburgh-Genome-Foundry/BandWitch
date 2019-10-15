
import os
from bandwitch import (SeparatingDigestionsProblem, list_common_enzymes,
                       LADDERS, load_record)


# DEFINE THE SEQUENCES AND THE ENZYME SET

# all common enzymes with a site of 6bp and at least 3 known providers
enzymes = list_common_enzymes(site_length=(6, ), min_suppliers=3)

records_folder = os.path.join('example_data', 'digestion_selection_data')
sequences = [
    load_record(os.path.join(records_folder, f), topology="circular", id=f)
    for f in sorted(os.listdir(records_folder))
]

# DEFINE AND SOLVE THE PROBLEM

problem = SeparatingDigestionsProblem(enzymes=enzymes,
                                      ladder=LADDERS['100_to_4k'],
                                      sequences=sequences,
                                      max_enzymes_per_digestion=1)
score, selected_digestions = problem.select_digestions(max_digestions=2)
print("Score: %.03f" % score)

# GENERATE A FIGURE OF THE BAND PATTERNS

axes = problem.plot_digestions(
    selected_digestions,
    patterns_props={'label_fontdict': {'rotation': 35}},
    target_file="separating_digestions.png"
)
problem.plot_distances_map(digestions=selected_digestions,
                           target_file="separating_digestions_distances.png")
