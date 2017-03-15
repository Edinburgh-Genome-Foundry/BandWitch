
import os
from bandwitch import SeparatingDigestionsProblem, LADDER_100_to_4k
from Bio import SeqIO, Restriction
from collections import OrderedDict


# DEFINE THE SEQUENCES AND THE ENZYME SET

# all common enzymes with a site of 6bp and at least 3 known providers
enzymes = sorted([
    str(enzyme) for enzyme in Restriction.CommOnly
    if (enzyme.size == 6) and (len(enzyme.supplier_list()) >= 3)
])

sequences = OrderedDict([
    (f, str(SeqIO.read(os.path.join('example_data', f), 'genbank').seq))
    for f in sorted(os.listdir('example_data'))
])

# DEFINE AND SOLVE THE PROBLEM

problem = SeparatingDigestionsProblem(sequences, enzymes, linear=False,
                                      ladder=LADDER_100_to_4k,
                                      max_enzymes_per_digestion=2,
                                      relative_error=0.05)
selected_digestions = problem.select_digestions()

# GENERATE A FIGURE OF THE BAND PATTERNS

axes = problem.plot_digestions(
    selected_digestions,
    patterns_props={'label_fontdict': {'rotation': 35}}
)
axes[0].figure.savefig("separating_digestions.png", bbox_inches="tight")

ax = problem.plot_distances_map(digestions=selected_digestions)
ax.figure.savefig("distances.png", bbox_inches="tight")
