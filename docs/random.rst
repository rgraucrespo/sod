.. _random-sampling:

Random sampling (randomsod)
===========================

``randomsod`` draws uniform random configurations at the target substitution
level and writes them as an ``ENSEMBLE``, with **no energy evaluation at all**.
It is the sampling counterpart of ``combsod``: where ``combsod`` enumerates
*every* symmetry-inequivalent configuration, ``randomsod`` draws a random sample
of them — the practical route when the target level is too large to enumerate.

Because the sample is purely geometric (every uniform draw is accepted),
energies are computed *a posteriori*: feed ``nXX/random/ENSEMBLE`` into the
normal structure-writer → DFT → ``statsod`` path, exactly as for an enumerated
``combsod`` ensemble.

Usage
-----

.. code-block:: bash

   sod_random.sh -nconfigs <N> [-symmetry on|off] [-seed clock|<int>]

Run it from the main problem directory, which must contain ``INSOD`` and
``SGO`` (always), plus ``EQMATRIX`` when ``-symmetry on`` (generate ``EQMATRIX``
with ``sod_comb.sh``).

Options
~~~~~~~

``-nconfigs <N>``
  Number of uniform **draws** (required).  This is the number of draws, not the
  number of distinct configurations — the latter is not known a priori.

``-symmetry on|off``  (default ``on``)
  With ``on``, draws are folded to symmetry representatives and the degeneracy
  column of the ``ENSEMBLE`` holds **visit counts**.  A uniform draw visits each
  symmetry orbit in proportion to its size, so visit counts are exactly the
  importance weights ``statsod`` needs for canonical Boltzmann averages — no
  separate orbit-size calculation is required.  With ``off``, every draw is
  written as its own row with degeneracy 1.

``-seed clock|<int>``  (default ``clock``)
  RNG seed.  Use a positive integer for a reproducible sample; ``clock`` seeds
  from the system clock for an independent sample each run.

Output
------

``nXX/random/ENSEMBLE`` (``XX`` = target level), in the standard ``ENSEMBLE``
format.  It is consumed by ``sod_sqs.sh`` (SQS selection), ``sod_gener.sh``
(calculation inputs), and ``statsod`` (canonical thermodynamics once energies
are available).

SQS selection from a random sample
-----------------------------------

Once ``nXX/random/ENSEMBLE`` exists, the standard SQS workflow runs entirely
from within ``nXX/random/``:

1. Place an ``INSQS`` file in ``nXX/random/`` (see :doc:`sqs_gqs` for the
   format).
2. Score the sample::

      cd nXX/random/
      sod_sqs.sh

   ``sqssod`` reads ``ENSEMBLE`` and ``INSQS`` from the current directory and
   writes ``OUTSQS`` there.  ``EQMATRIX``, ``supercell.cif``, and ``INSOD``
   are taken from SODPROJECT/ automatically.

3. Pick the best SQS from ``OUTSQS`` (rank 0) and generate calculator input
   files for it::

      sod_gener.sh -choose <index>

   Calculation subdirectories (``c01/``, etc.) are created under
   ``nXX/random/``, keeping the random-sample workflow self-contained.

Relationship to Metropolis MC
-----------------------------

``randomsod`` and ``mcsod`` are complementary:

- ``randomsod`` samples geometry uniformly with no Hamiltonian, ideal for
  building large candidate sets (e.g. for SQS/GQS selection) before any DFT.
- ``mcsod`` (see :doc:`pme_mc`) runs an energy-biased Metropolis walk with the
  PME Hamiltonian to sample the low-energy, finite-temperature ensemble.
