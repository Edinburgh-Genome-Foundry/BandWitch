"""
Complex validation test.
"""

import matplotlib

matplotlib.use("Agg")

import os
from bandwitch import ClonesObservations


def test_complex_validation(tmpdir):
    data_dir = os.path.join("tests", "test_data", "complex_validation_data")
    records_path = os.path.join(data_dir, "constructs_sequences.zip")
    constructs_map_path = os.path.join(data_dir, "constructs_map.xls")
    digestions_map_path = os.path.join(data_dir, "digestions_map.xls")
    aati_zip_path = os.path.join(data_dir, "digestion_results.zip")
    clones = ClonesObservations.from_files(
        records_path=records_path,
        constructs_map_path=constructs_map_path,
        aati_zip_path=aati_zip_path,
        digestions_map_path=digestions_map_path,
    )
    validations = clones.validate_all_clones(relative_tolerance=0.03)
    clones.plot_validations_plate_map(validations)
    clones.plot_all_validations_patterns(validations)
    partial_digest_analysis = clones.partial_digests_analysis()
    clones.plot_partial_digests_analysis(partial_digest_analysis)
    summary_table = clones.validations_summary_table(
        validations, target=os.path.join(str(tmpdir), "table.csv")
    )


# Chosen digestions, for reference:
# C00	EcoRI,PvuI
# C01	EcoRI,PvuI
# C02	EcoRI,PvuI
# C03	EcoRI,PvuI
# C04	EcoRI,PvuI
# C05	EcoRI,PvuI
# C06	AseI,EcoRI
# C07	AseI,EcoRI
# C08	AseI,EcoRI
