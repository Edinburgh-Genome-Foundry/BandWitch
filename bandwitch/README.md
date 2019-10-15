# Code organization

This document walks you trough the DNA Chisel code. Please request changes if anything is unclear.

#### Ladder/

- **Ladder.py** implements the Ladder class, used throughout the library to represent band ladders.A ladder's main method is ``dna_size_to_migration``, which predicts the coordinates of a band of a given size.
- **preset_ladders.py** provides a set of commonly used ladders obtained by parsing the spreadsheets in ``data/``.

#### DigestionProblem/

This submodule implements the classes and methods for finding enzyme sets adapted to a set of sequences.
- **SetCoverProblem** (implemented in SetCoverProblem.py) implements a class for solving generic Minimal Subset Cover problems. This class inherited by all other classes in this module to find one or more digestions collectively covering all the given records.
- **DigestionProblem** (implemented in DigestionProblem.py) inherits from *SetCoverProblem* and is the base class for all enzyme selection problems. It implements the initialization (where the digestion of every sequence by every enzyme mixes are computed), the plotting of the sequences by selected enzymes, and the scaffold for enzyme selection. All children classes simply need to define a custom ``_parameter_element_score`` method.
- **IdealDigestionsProblem** (implemented in IdealDigestionsProblem.py) is a subclass of DigestionProblem to select enzymes that will give "acceptable" patterns for every sequence provided. Its ``_parameter_element_score`` method computes a ``migration_score()`` depending on the number of bands, the spacing between bands, etc.
- **SeparatingDigestionsProblem** (implemented in SeparatingDigestionsProblem.py) is a subclass of DigestionProblem to select enzymes that will give distinct patterns (not necessarily ideal) for each given sequence. It's ``_parameter_element_score`` computes the distance between the bands of 2 digestions.

#### ClonesObservations/

This submodule implements classes for representing and validating the results of experimental restriction digest experiments. In particular, it allows to import observations from the  AATI fragment analyzer, and compare these to expected results to create validation or identification reports.

The class structure in this module is complicated because life is complicated. In a typical assembly batch, there are several different constructs assembled. For each construct we can pick several clones and submit the DNA extracted from each clone to several restriction digests.

- **BandsObservations** represents the observation of a digest. It has a name (e.g. the name of the microplate well in which it was observed), bands (more precisely, a set of bands sizes), a ladder with which it was observed, and optionally a picture of the gel. It has a method to be compared with another pattern (using the *band_patterns_discrepancy()* method defined in *band_patterns_discrepancy.py*)
- **Clone** represents a clone, with a name (could be a microplate well's name), the construct associated with the clone (if known) and the (possibly multiple) restriction digest observations for this clone.
- **CloneValidation** is a structure representing the comparison between the observations of a Clone and expected patterns. It contains a Clone, the expected patterns for each digestion, and the discrepancy between observed and expected for each digestion. It is the "building block" of validation and identification reports.
- **ClonesObservations** represents all the information necessary to process a full experimental outcome involving several clones and restriction digests. At it's core, it groups several Clone instances together with the biopython records representing the expected constructs. As the highest-level class in this module, ClonesObservations has methods to be directly created from an AATI Fragment Analyzer zip file (or several zip files, as ClonesObservations can be merged togther). The class also has methods to generate CloneValidations for all clones in its set, in order to validate or identify all clones and compile the results in a PDF report.

#### list_common_enzymes/

- **list_common_enzymes.py** provides a method for getting common enzymes (filtered for star-activity using data in *enzymes_data/*)
- **enzymes_data/** provides the dictionnary ``enzymes_infos``, generated in the *__init__.py* by parsing the spreadsheet **enzymes_infos.csv**, which contains enzyme data, notably methylation sensitivity, compiled from rebase.neb.com using the script **update_enzymes_list.py**.

#### Files at the root

- **bands_predictions.py** implements methods to predict which band sizes given records and enzymes will create, in particular, some methods efficiently compute band sizes for batches of sequences and combinations of enzymes. The methods are used extensively by the *DigestionProblem* class.
- **tools.py** implements generic methods, notably for Genbank record manipulation.