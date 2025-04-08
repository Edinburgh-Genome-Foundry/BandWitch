"""
Basic tests to check that the core functionalities are at least running.
"""

import matplotlib

matplotlib.use("Agg")

import os
from collections import OrderedDict
from bandwitch import (
    SeparatingDigestionsProblem,
    IdealDigestionsProblem,
    LADDERS,
    list_common_enzymes,
)
from Bio import SeqIO, Restriction
import pytest


@pytest.fixture
def sequences():
    data_path = os.path.join("tests", "test_data", "test_sequences")
    return OrderedDict(
        [
            (
                fname,
                str(SeqIO.read(os.path.join(data_path, fname), "genbank").seq),
            )
            for fname in sorted(os.listdir(data_path))
        ]
    )


def test_separating_digestions(tmpdir, sequences):
    enzymes = sorted(
        [
            str(enzyme)
            for enzyme in Restriction.CommOnly
            if (enzyme.size == 6) and (len(enzyme.supplier_list()) >= 3)
        ]
    )
    problem = SeparatingDigestionsProblem(
        sequences=sequences,
        enzymes=enzymes,
        topology="circular",
        ladder=LADDERS["100_to_4k"],
        max_enzymes_per_digestion=1,
        relative_migration_precision=0.05,
    )
    score, selected_digestions = problem.select_digestions(max_digestions=1)
    assert selected_digestions == [("AvaI",)]

    axes = problem.plot_digestions(
        selected_digestions,
        patterns_props={"label_fontdict": {"rotation": 35}},
    )
    axes[0].figure.savefig(
        os.path.join(str(tmpdir), "separating_digestions.png"),
        bbox_inches="tight",
    )
    ax = problem.plot_distances_map(digestions=selected_digestions)
    ax.figure.savefig(os.path.join(str(tmpdir), "distances.png"))


def test_ideal_digestions(sequences):
    class MyIdealDigestionsProblem(IdealDigestionsProblem):
        def migration_pattern_is_ideal(self, migration):
            """2-3 bands between 30% and 70% of the migration span?"""
            min_migr = 0.7 * self.migration_min + 0.3 * self.migration_max
            max_migr = 0.3 * self.migration_min + 0.7 * self.migration_max
            bands_in_central_zone = [
                band for band in migration if min_migr <= band <= max_migr
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
    ]

    # DEFINE AND SOLVE THE PROBLEM

    problem = MyIdealDigestionsProblem(
        sequences=sequences,
        enzymes=enzymes,
        ladder=LADDERS["100_to_4k"],
        topology="circular",
        max_enzymes_per_digestion=2,
    )
    score, selected_digestions = problem.select_digestions(max_digestions=1)
    assert score > 0.09
    assert selected_digestions == [("BamHI", "SacI")]


def test_list_common_enzymes():
    l = list_common_enzymes(min_suppliers=4, avoided_methylations=("Dam", "Dcm"))
    assert len(l) == 38
    l = list_common_enzymes(min_suppliers=4, avoided_methylations=())
    assert len(l) == 45
