BandWitch Reference manual
==========================

DigestionProblem
--------------------
.. mermaid::

  graph TD
    SC[SetCoverProblem] --> DP[DigestionProblem]
    DP --> SDP[SeparatingDigestionProblem]
    DP --> IDP[IdealDigestionProblem]

.. automodule:: bandwitch.DigestionProblem
   :members:

Clone Observations
------------------

.. mermaid::

  graph TD;
    aati[AATI fragment analysis file] -- parser --> bobs[BandsObervation-s];
    bobs --> clone[Clone-s]
    bobs --> clone
    digestions[digestion enzymes infos] -- parser --> clone

    clone --> ClonesObservations
    clone --> ClonesObservations
    cs[constructs sequences] --> ClonesObservations
    ClonesObservations --pdf generator--> reports
    style aati stroke:none, fill: none
    style digestions stroke:none, fill: none
    style cs stroke:none, fill: none
    style reports stroke:none, fill: none

Note: the classes in this module have a complicated organization, mostly due
to the history of this module and the heterogeneity of the sources of data
necessary for clone validation. It may get better in the future.

.. automodule:: bandwitch.ClonesObservations
   :members:

Ladder
---------

.. automodule:: bandwitch.Ladder
   :members:

Bands Predictions
-----------------
.. mermaid::

  graph TD
    pcomp[_compute_digestion_bands] --> pr[predict_digestion_bands]
    pcomp --> pds[predict_sequences_digestions]

.. automodule:: bandwitch.bands_predictions
   :members:

SetCoverProblem
---------------

.. automodule:: bandwitch.SetCoverProblem
   :members:
