.. reference ::

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


SetCoverProblem
---------------
