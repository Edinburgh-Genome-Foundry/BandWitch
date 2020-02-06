.. raw:: html

    <p align="center">
    <img alt="BandWitch Logo" title="BandWitch Logo" src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/BandWitch/master/docs/_static/images/title.png" width="400">
    <br /><br />
    </p>

Bandwitch
=========

.. image:: https://travis-ci.org/Edinburgh-Genome-Foundry/BandWitch.svg?branch=master
   :target: https://travis-ci.org/Edinburgh-Genome-Foundry/BandWitch
   :alt: Travis CI build status
.. image:: https://coveralls.io/repos/github/Edinburgh-Genome-Foundry/BandWitch/badge.svg?branch=master
   :target: https://coveralls.io/github/Edinburgh-Genome-Foundry/BandWitch?branch=master

Bandwitch (full documentation `here <https://edinburgh-genome-foundry.github.io/BandWitch/>`_)
is a Python library for the planning and analysis of restriction
experiments in DNA assembly operations. Bandwitch implements method to select the best enzyme(s) to validate or identify DNA assemblies. It also provides report generation methods to automatically validate/identify assemblies from experimental data.

You can try BandWitch's enzyme suggestion feature in `this web demo <https://cuba.genomefoundry.org/select_digestions>`_, and the sequence validation (from AATI fragment analyzer files) in `this other demo <http://cuba.genomefoundry.org/analyze-digests>`_

Installation
-------------

You can install DnaCauldron through PIP


.. code:: shell

    sudo pip install bandwitch

On Ubuntu at least, you may need to install libblas first:

.. code::

    sudo apt-get install libblas-dev liblapack-dev

Alternatively, you can unzip the sources in a folder and type

.. code:: shell

    sudo python setup.py install


Enzyme selection with BandWitch
-------------------------------

In the following examples, we assume that we have a set of 12 constructs which we will
need to either validate (i.e. we digest these constructs and compare each pattern
with the expected pattern for that construct) or identify (i.e. we will digest an
a-priori unknown construct and use the migration patterns to un-ambiguously
identify each construct among the 12 possible candidates).

For validation purposes, the difficulty is to find a digestion that will produce
harmonious patterns for all the constructs at once: well-spaced bands, and not
too many or too few of them. For identification purposes, the difficulty is to
find a digestion giving very distant patterns for each construct in the set of
candidates.

Every time when the problem cannot be solved with a single digestion, BandWitch
can propose 2 or 3 digestions which collectively solve the problem.

**Important:** when providing BandWitch with a record, make sure to set the
topology, defined by ``record.annotations['topology'] = 'linear'|'circular'``.


Finding enzymes that "work well" for many constructs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is the code to select enzymes that will produce nice patterns for all constructs, for validation:

.. code:: python

  from bandwitch import IdealDigestionsProblem, LADDERS, load_record

  # DEFINE THE SEQUENCES AND THE ENZYME SET
  enzymes = ["EcoRI", "BamHI", "XhoI", "EcoRV", "SpeI", "XbaI",
             "NotI", "SacI", "SmaI", "HindIII", "PstI"]
  sequences = [
      load_record(genbank_file_paidname=f, topology='circular')
      for genbank_file_path in some_list_of_files)
  ]

  # SELECT THE BEST SINGLE DIGESTION WITH AT MOST ENZYMES
  problem = IdealDigestionsProblem(enzymes=enzymes,
                                   ladder=LADDERS['100_to_4k'],
                                   sequences=sequences,
                                   max_enzymes_per_digestion=2)
  score, selected_digestions = problem.select_digestions(max_digestions=1)

  # PLOT THE BAND PATTERNS PRODUCED BY THE SELECTED DIGESTION
  problem.plot_digestions(
      digestions=selected_digestions,
      patterns_props={'label_fontdict': {'rotation': 35}},
      target_file="ideal_digestions.png"
  )

Result:

.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/BandWitch/master/examples/ideal_digestions.png
   :alt: [logo]
   :align: center

Finding enzymes that will differentiate many constructs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To select enzymes that will produce **different patterns for each construct, for identification:**

.. code:: python

    from bandwitch import (SeparatingDigestionsProblem, list_common_enzymes,
                           LADDERS, load_record)


    # DEFINE SEQUENCES AND ENZYME SET (6-CUTTERS WITH >3 COMMERCIAL PROVIDERS)
    enzymes = list_common_enzymes(site_length=(6,), min_suppliers=3)
    sequences = [
        load_record(genbank_file_path, id=f)
        for genbank_file_path in some_list_of_files)
    ]

    # SELECT THE BEST DIGESTION PAIRS (AT MOST 1 ENZYME PER DIGESTION)
    problem = SeparatingDigestionsProblem(enzymes=enzymes,
                                          ladder=LADDERS['100_to_4k'],
                                          sequences=sequences,
                                          max_enzymes_per_digestion=1)
    score, selected_digestions = problem.select_digestions(max_digestions=2)

    # GENERATE A FIGURE OF THE BAND PATTERNS

    problem.plot_digestions(
        selected_digestions,
        patterns_props={'label_fontdict': {'rotation': 35}},
        target_file="separating_digestions.png"
    )

    problem.plot_distances_map(digestions=selected_digestions,
                               target_file="separating_digestions_distances.png")

Result:

.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/BandWitch/master/examples/separating_digestions.png
   :alt: [logo]
   :align: center

In the result above, each construct has a unique "fingerprint". Assuming that you
have an unlabelled DNA sample which could be any of these assemblies, then simply
digesting the sample with MspA1I and BsmI will give you 2 pattern which collectively
will correspond to a unique assembly.

Usage: Construct validation or identification from experimental data
---------------------------------------------------------------------

This part is still under construction.

Bandwitch can process output files from an automated fragment analyzer and produce
informative reports as illustrated below:

.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/BandWitch/master/docs/_static/images/bands_validation.png
   :alt: [logo]
   :align: center
   :width: 750px


License = MIT
--------------

BandWitch is an open-source software originally written at the `Edinburgh Genome Foundry <http://edinburgh-genome-foundry.github.io/home.html>`_ by `Zulko <https://github.com/Zulko>`_ and `released on Github <https://github.com/Edinburgh-Genome-Foundry/Primavera>`_ under the MIT licence (¢ Edinburg Genome Foundry). Everyone is welcome to contribute !

More biology software
-----------------------

.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/Edinburgh-Genome-Foundry.github.io/master/static/imgs/logos/egf-codon-horizontal.png
   :target: https://edinburgh-genome-foundry.github.io/

BandWitch is part of the `EGF Codons <https://edinburgh-genome-foundry.github.io/>`_ synthetic biology software suite for DNA design, manufacturing and validation.

.. _Github: https://github.com/EdinburghGenomeFoundry/BandWitch
.. _PYPI: https://pypi.python.org/pypi/bandwitch
