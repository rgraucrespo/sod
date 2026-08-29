Special Quasirandom Structures (SQS and GQS)
============================================

SQS and GQS are methods for identifying configurations that best represent
random disorder in a solid solution, by comparing short-range ordering
statistics to ideal random-mixing values.

Concept
-------

For a disordered solid solution, ideal random mixing produces characteristic
pair probabilities. The ``sqssod`` program scores and ranks all configurations
by comparing species-resolved pair probabilities against independent-random
targets. For a pair orbit connecting targets ``t`` and ``u``, the target for
species channel ``alpha,beta`` is ``x(t,alpha) * x(u,beta)``. This target-aware
form covers binary, multinary, multi-target, and multinary multi-target
substitution without allowing species to move onto forbidden sublattices.

- **SQS** (``sqssod``): identifies the best-matching configuration at a fixed
  composition (0 K weighting).
- **GQS** (``gqssod``): extends SQS to finite temperature by computing
  Boltzmann-weighted thermal averages of the multi-site cluster correlation
  functions (orders 1 to ``MaxOrder``) from calculated energies.

Ensemble source: enumeration or random sampling
------------------------------------------------

``sqssod`` and ``gqssod`` score whatever configurations are listed in
``ENSEMBLE`` — they do **not** enumerate or search the configurational space
themselves. The SQS is simply the best-scoring member of the ensemble you
provide. That ensemble can come from either route:

- **Full enumeration** (``sod_comb.sh`` → ``nXX/ENSEMBLE``): the
  symmetry-inequivalent configurations at the composition. This guarantees the
  true best SQS *within the supercell* is present, but is only feasible when the
  space is small enough to enumerate.
- **Uniform random sampling** (``sod_random.sh`` →
  ``nXX/random/ENSEMBLE``): a large random sample of configurations. This is
  the practical route when the full space is **too large to enumerate** —
  generate a big uniform ensemble and extract the best SQS from it. Point
  ``sqssod`` at the ``random/`` directory (which holds the sampled ``ENSEMBLE``)
  while ``EQMATRIX``, ``supercell.cif``, and ``INSOD`` remain in the parent
  folder. Random sampling needs no reference energies, so this works before any
  DFT is run; for GQS, supply ``ENERGIES`` for the sampled configurations.

With random sampling the result is the best SQS *found in the sample*, not
provably the global optimum — enlarge the sample to improve it.

.. note::

   **Future direction.** SOD does not yet perform a *directed* search for the
   SQS (e.g. simulated annealing or basin-hopping that minimises the score
   ``Q`` as a pseudo-energy to drive the configuration toward ideal
   randomness). At present the SQS can only be selected from a pre-generated
   ensemble (enumerated or randomly sampled). A directed-search mode is a
   planned enhancement.

Workflow
--------

**From a full enumeration** (``sod_comb.sh`` → ``nXX/ENSEMBLE``):

Prerequisites: ``sod_comb.sh`` must have been run to generate ``EQMATRIX``,
``supercell.cif``, and ``nXX/ENSEMBLE``.

1. Create an ``INSQS`` file in SODPROJECT/ (or in ``nXX/`` for a per-composition
   override).
2. Run SQS scoring::

      sod_sqs.sh        # from SODPROJECT/: scores all nXX/ compositions
      sod_sqs.sh        # from nXX/:        scores that composition only

3. Optionally run GQS (requires ``ENERGIES`` and ``TEMPERATURES``)::

      sod_gqs.sh        # from SODPROJECT/ or nXX/

**From a random sample** (``sod_random.sh`` → ``nXX/random/ENSEMBLE``):

Prerequisites: ``sod_comb.sh`` must have been run (for ``EQMATRIX`` and
``supercell.cif``); ``sod_random.sh`` must have produced ``nXX/random/ENSEMBLE``.

1. Create ``INSQS`` in ``nXX/random/`` (alongside ``ENSEMBLE``).
2. Run SQS scoring from that directory::

      cd nXX/random/
      sod_sqs.sh        # reads ENSEMBLE and INSQS here; writes OUTSQS here

3. Inspect ``OUTSQS`` and pick the best configuration (rank 1 is nearest to
   ideal randomness).  Generate calculator input files for the chosen index::

      sod_gener.sh -choose <index>   # still from nXX/random/
                                     # writes cYY/ under nXX/random/

Input file: INSQS
-----------------

``INSQS`` controls the pair cutoff and scoring parameters. Example::

   # Maximum cluster order (2-6)
   4

   # Cutoff radii (Angstroms) for orders 2..MaxOrder
   8.0 6.0 4.0

   # Weights for orders 2..MaxOrder
   1.0 1.0 1.0

   # omega eps_tol
   10  1.0E-8

   # n_top_sqs: number of top configurations listed in OUTSQS (0 = rank and list all)
   10

Key parameters:

- **MaxOrder**: accepted for backward compatibility. The generalized
  ``sqssod`` scorer currently uses pair correlations only, so order 2 controls
  the SQS score and entries for orders 3+ are ignored by ``sqssod``.
- **Cutoff radii**: one per order from 2 to MaxOrder, in Å. ``sqssod`` uses the
  order-2 value as the pair cutoff.
- **Weights**: one per order from 2 to MaxOrder. ``sqssod`` uses the order-2
  weight in the pair-family error normalization.
- **Scoring**: van de Walle matched-diameter scoring, ``Q = -omega * L +
  Error``, where ``Error`` is the normalized mean absolute deviation of
  species-resolved pair probabilities from their target random probabilities.
- **n_top_sqs**: number of best-ranked configurations listed in ``OUTSQS``;
  ``0`` sorts and lists the whole ensemble. Optional line; defaults to 10 if
  absent.

Output files
------------

``sqssod`` writes two files in each composition or sample folder:

- ``OUTSQS``: compact ranked list of the ``n_top_sqs`` best configurations
  (default 10; ``0`` = all) with matched distance, total error, score, and
  target-pair family errors such as ``E11`` and ``E12``. Equal-score
  configurations rank by configuration index.
- ``SQS_CORRELATIONS``: detailed species-resolved pair-channel probabilities
  for the best-ranked configuration, including ``P_cell``, ``P_random``, and
  ``Delta`` for each pair orbit/channel.

Rank 1 is the best SQS.

``gqssod`` writes two files:

- ``OUTGQS``: thermal averages of the cluster correlation functions at each
  temperature in ``TEMPERATURES``, plus the T → ∞ (equiprobable) limit.
- ``wc_parameters.dat``: the corresponding Warren-Cowley short-range order
  parameters α\ :sub:`n` for each symmetrically distinct pair shell. These are
  reported as an SRO diagnostic; the SQS/GQS ranking itself is done on the
  correlation functions, not on α\ :sub:`n`.

See the repository `README.md <https://github.com/rgraucrespo/sod/blob/master/README.md>`_ for full format descriptions of
``OUTSQS``, ``SQS_CORRELATIONS`` and ``OUTGQS``.

Selecting a SQS
---------------

- **Rank 1** is closest to ideal randomness. Prefer the configuration with
  the lowest weighted van de Walle score, ``Q``.
- **Degeneracy**: higher-degeneracy configurations may be thermodynamically
  preferred; use judgement when multiple configurations score similarly.
- **Supercell size**: larger supercells generally yield better SQS. For small
  cells, no configuration may achieve exact π = target.

Further information
-------------------

The repository `README.md <https://github.com/rgraucrespo/sod/blob/master/README.md>`_ contains the full
``INSQS``/``OUTSQS``/``OUTGQS`` format specifications and a worked example
using ``example01/FILER1_gulp``.
