"""
Basic tests to check that the core functionalities are at least running.
"""

import os
from collections import OrderedDict
from bandwitch import (SeparatingDigestionsProblem, IdealDigestionsProblem,
                       LADDER_100_to_4k)
from Bio import SeqIO, Restriction
import pytest


@pytest.fixture
def sequences():
    data_path = os.path.join('tests', 'test_data')
    return OrderedDict([
        (fname, str(SeqIO.read(os.path.join(data_path, fname), 'genbank').seq))
        for fname in sorted(os.listdir(data_path))
    ])

def test_separating_digestions(tmpdir, sequences):
    enzymes = sorted([
        str(enzyme) for enzyme in Restriction.CommOnly
        if (enzyme.size == 6) and (len(enzyme.supplier_list()) >= 3)
    ])
    problem = SeparatingDigestionsProblem(sequences, enzymes, linear=False,
                                          ladder=LADDER_100_to_4k,
                                          max_enzymes_per_digestion=2,
                                          relative_error=0.05)
    selected_digestions = problem.select_digestions()
    assert selected_digestions == [('AvaI',)]

    axes = problem.plot_digestions(
        selected_digestions,
        patterns_props={'label_fontdict': {'rotation': 35}}
    )
    axes[0].figure.savefig(
        os.path.join(str(tmpdir), "separating_digestions.png"),
        bbox_inches="tight"
    )
    ax = problem.plot_distances_map(digestions=selected_digestions)
    ax.figure.savefig(os.path.join(str(tmpdir), "distances.png"))


def test_ideal_digestions(sequences):

    class MyIdealDigestionsProblem(IdealDigestionsProblem):
        def migration_pattern_is_ideal(self, migration):
            """2-3 bands between 30% and 70% of the migration span?"""
            min_migr = 0.7 * self.migration_min + 0.3 * self.migration_max
            max_migr = 0.3 * self.migration_min + 0.7 * self.migration_max
            bands_in_central_zone = [band for band in migration
                                     if min_migr <= band <= max_migr]
            return 2 <= len(bands_in_central_zone) <= 3

    # DEFINE THE SEQUENCES AND THE ENZYME SET

    enzymes = ["EcoRI", "BamHI", "XhoI", "EcoRV", "SpeI", "XbaI", "NotI",
               "SacI", "SmaI", "HindIII", "PstI"]

    # DEFINE AND SOLVE THE PROBLEM

    problem = MyIdealDigestionsProblem(sequences, enzymes, linear=False,
                                       ladder=LADDER_100_to_4k,
                                       max_enzymes_per_digestion=2)
    selected_digestions = problem.select_digestions()
    assert selected_digestions == [('EcoRI', 'SmaI')]
