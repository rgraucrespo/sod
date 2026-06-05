Special Quasirandom Structures (SQS and GQS)
============================================

SQS and GQS are methods for identifying configurations that best represent
random disorder in a solid solution, by comparing pair correlations (Warren
parameters) to ideal random-mixing values.

Concept
-------

For a disordered solid solution at composition *x*, ideal random mixing
produces a characteristic set of pair correlation values π for each cluster
distance. The ``sqssod`` program scores and ranks all enumerated configurations
by their deviation from these ideal values — the configuration with the
smallest deviation is the Special Quasirandom Structure (SQS).

- **SQS** (``sqssod``): identifies the best-matching configuration at a fixed
  composition (0 K weighting).
- **GQS** (``gqssod``): extends SQS to finite temperature by computing
  Boltzmann-weighted thermal averages of pair correlations from calculated
  energies.

Workflow
--------

Prerequisites: run ``sod_comb.sh`` first to generate ``EQMATRIX``,
``supercell.cif``, and ``n*/ENSEMBLE``.

1. Create an ``INSQS`` file in the working directory (see below).
2. Run SQS scoring::

      sod_sqs.sh        # from SODPROJECT: scores all nXX/ compositions
      sod_sqs.sh        # from nXX/: scores that composition only

3. Optionally run GQS (requires ``ENERGIES`` and ``TEMPERATURES``)::

      sod_gqs.sh        # from SODPROJECT or nXX/

Input file: INSQS
-----------------

``INSQS`` controls which clusters are included and how configurations are
scored. Example::

   # Maximum cluster order (2-6)
   4

   # Cutoff radii (Angstroms) for orders 2..MaxOrder
   8.0 6.0 4.0

   # Weights for orders 2..MaxOrder
   1.0 1.0 1.0

   # omega eps_tol
   10  1.0E-8

Key parameters:

- **MaxOrder**: maximum cluster order (2–6).
- **Cutoff radii**: one per order from 2 to MaxOrder, in Å.
- **Weights**: relative importance of each order from 2 to MaxOrder in the
  normalized weighted ``AbsErr``; set to 0.0 to remove that order from
  ``AbsErr``.
- **Scoring**: van de Walle matched-diameter scoring, ``Q = -omega * L +
  AbsErr``, where ``AbsErr`` is the normalized weighted mean of \|Δπ\|.

Output files
------------

``sqssod`` writes ``OUTSQS`` in each composition folder — a ranked list of
configurations with their pair correlations and scores. Rank 0 is the best
SQS.

``gqssod`` writes ``OUTGQS`` — thermal averages of pair correlations at each
temperature in ``TEMPERATURES``.

See the ``README.md`` for full format descriptions of ``OUTSQS`` and
``OUTGQS``.

Selecting a SQS
---------------

- **Rank 0** is closest to ideal randomness. Prefer the configuration with
  the lowest weighted van de Walle score, ``Q``.
- **Degeneracy**: higher-degeneracy configurations may be thermodynamically
  preferred; use judgement when multiple configurations score similarly.
- **Supercell size**: larger supercells generally yield better SQS. For small
  cells, no configuration may achieve exact π = target.

Further information
-------------------

The ``README.md`` contains complete documentation including the full
``INSQS``/``OUTSQS``/``OUTGQS`` format specifications and a worked example
using ``example01/FILER1_gulp``.
