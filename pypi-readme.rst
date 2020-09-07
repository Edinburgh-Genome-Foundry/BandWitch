BandWitch
=========

BandWitch is a Python library for the planning and analysis of restriction
experiments in DNA assembly operations. BandWitch implements method to select
the best enzyme(s) to validate or identify DNA assemblies. It also provides
report generation methods to automatically validate/identify assemblies from
experimental data.

You can try BandWitch's enzyme suggestion feature in
`this web demo <http://cuba.genomefoundry.org/digestion-selector>`_,
and the sequence validation (from AATI fragment analyzer files) in
`this other demo <http://cuba.genomefoundry.org/analyze-digests>`_

Infos
------

**PIP installation:**

.. code:: bash

   pip install bandwitch

**Web documentation:**

`<https://edinburgh-genome-foundry.github.io/BandWitch/>`_

**Github Page:**

`<https://github.com/Edinburgh-Genome-Foundry/BandWitch>`_

**Live demo:**

Enzyme suggestion: `<http://cuba.genomefoundry.org/digestion-selector>`_

Fragment analysis validation: `<http://cuba.genomefoundry.org/analyze-digests>`_

**License:** MIT, Copyright Edinburgh Genome Foundry

Screenshots
-----------

Enzymes selected by bandwitch to obtain **clear, optimal** patterns for all tested DNA constructs:

.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/BandWitch/master/examples/ideal_digestions.png
   :alt: [logo]
   :align: center

Enzymes selected by bandwitch to obtain **significant differences** between the patterns of the tested constructs, so that a construct can be identified by its pattern:

.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/BandWitch/master/examples/separating_digestions.png
   :alt: [logo]
   :align: center

More biology software
-----------------------

.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/Edinburgh-Genome-Foundry.github.io/master/static/imgs/logos/egf-codon-horizontal.png
   :target: https://edinburgh-genome-foundry.github.io/

BandWitch is part of the `EGF Codons <https://edinburgh-genome-foundry.github.io/>`_ synthetic biology software suite for DNA design, manufacturing and validation.
