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

``randomsod`` supports the same general substitution model as ``combsod``:
multiple target sublattices (multisite) and multiple substituting species on a
single site (multinary).  See :ref:`random-multinary` below.

Usage
-----

.. code-block:: bash

   sod_random.sh -nconf <N> [-sym on|off] [-seed clock|<int>]

Run it from the main problem directory, which must contain ``INSOD`` and
``SGO`` (always), plus ``EQMATRIX`` when ``-sym on`` (generate ``EQMATRIX``
with ``sod_comb.sh``).

Options
~~~~~~~

``-nconf <N>``
  Number of uniform **draws** (required).  This is the number of draws, not the
  number of distinct configurations — the latter is not known a priori.

``-sym on|off``  (default ``off``)
  With ``on``, draws are folded to symmetry representatives and the degeneracy
  column of the ``ENSEMBLE`` holds **visit counts**.  A uniform draw visits each
  symmetry orbit in proportion to its size, so visit counts are exactly the
  importance weights ``statsod`` needs for canonical Boltzmann averages — no
  separate orbit-size calculation is required.  With ``off``, every draw is
  written as its own row with degeneracy 1.  Folding uses a *colored* orbit
  representative, so it is correct for multinary and multisite runs (a
  configuration is only merged with another if their sites match species by
  species, not merely as an unlabelled set).  ``-sym on`` requires ``EQMATRIX``,
  which ``sod_comb.sh`` writes with one block per target sublattice.

``-seed clock|<int>``  (default ``clock``)
  RNG seed.  Use a positive integer for a reproducible sample; ``clock`` seeds
  from the system clock for an independent sample each run.

Output
------

``nXX/random/ENSEMBLE`` in the standard ``ENSEMBLE`` format, where ``XX`` is the
substitution level.  The directory name follows ``combsod``: for a single binary
target it is ``n<nsubs>``; for multinary or multisite runs the per-species counts
are joined with underscores (e.g. ``n04_04`` for a 4/4 multinary split,
``n02_01`` for two sublattices substituted 2 and 1).  It is consumed by
``sod_sqs.sh`` (SQS selection), ``sod_gener.sh`` (calculation inputs), and
``statsod`` (canonical thermodynamics once energies are available).

.. _random-multinary:

Multinary and multisite substitutions
--------------------------------------

The number and identity of targets and substituting species are read from
``INSOD`` exactly as for ``combsod`` — ``randomsod`` needs no extra flags. Each
draw independently places, on **every** target sublattice, a uniform *colored*
assignment of the requested per-species counts:

- **Multisite** (several ``sptarget`` sublattices): every sublattice is
  substituted in the same draw, at its own requested level.
- **Multinary** (``nsubs`` line with more than one count, e.g. ``4 4``): the
  substituted sites on that sublattice are partitioned uniformly among the
  substituting species.

The resulting ``ENSEMBLE`` has one column group per (target, species) and is
byte-compatible with the ``combsod`` format, so ``statsod`` and the structure
writers consume it unchanged. With ``-sym on`` the same crystal symmetry acts on
all sublattices simultaneously, and the visit counts reproduce the exact orbit
degeneracies from ``combsod`` (in the enumerable limit) up to Monte-Carlo noise.

The colored orbit folding uses the **same canonical representative as**
``combsod`` — the ranking is nested target-major and, within each target,
support-major (sorted combined support, then each species colouring). As a
result, for any orbit that appears in both, ``randomsod -sym on`` and ``combsod``
write the *identical* configuration row: the two ``ENSEMBLE`` files match
line-for-line on the sampled orbits, differing only in the degeneracy column
(``combsod`` writes exact orbit sizes, ``randomsod`` writes visit counts).

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
