Periodic Motif Expansion (PME) and Monte Carlo (MC)
====================================================

PME builds an effective Hamiltonian for a disordered crystal by decomposing
its energy into contributions from periodic motifs.  The fitted Hamiltonian can
then be evaluated for every configuration at a target composition, or used by
the MC sampler when full enumeration is impractical.

Concept
-------

The energy of a configuration is written as a sum of site interaction energies
up to four-body motifs:

.. code-block:: none

   E = V0 + sum_i V1(i) + sum_{i<j} V2(i,j)
       + sum_{i<j<k} V3(i,j,k) + ...

The interaction terms are fitted from :term:`reference energies` at low
substitution levels (``n00`` to ``n04``) and, when available, from symmetric
high-substitution levels (``nM`` to ``n(M-4)``), where ``M`` is the number of
substitutable sites.

Three Hamiltonian variants are used:

- ``PME0``: low-side expansion only.
- ``PME1``: high-side, or hole, expansion only.
- ``PMEh``: weighted hybrid of ``PME0`` and ``PME1``.

The per-order ``epsilon`` scale factors and the hybrid ``alpha`` / ``eta``
parameters are controlled by ``SODPROJECT/pme.model`` for the target level
specified in ``INSOD``.  Optional calibration energies are read from
``nXX/ENERGIES``.

PMEh Hybrid Weighting
---------------------

``PMEh`` blends the low-side energy estimate :math:`E^{(0)}(x)` and the
high-side estimate :math:`E^{(N)}(x)` using a piecewise power-law scheme
parameterised by two values in ``pme.model``:

- **alpha** (default 2.0): sharpness of the composition-space transition.
  ``alpha > 1`` guarantees continuous first derivatives at the blend boundaries.
  ``alpha = 2`` gives a smooth sigmoidal curve.
- **eta** (default 0.0): log-scale asymmetry.  ``eta > 0`` keeps the low-end
  PME dominant for a larger composition range; ``eta < 0`` does the same for
  the high-end PME.

Let :math:`x = n/N` be the composition (substitutions / total sites) and
:math:`x_\mathrm{ref} = p/N` where :math:`p` is the PMEOrder truncation.

For :math:`x \le x_\mathrm{ref}` the low-end PME is within its explicitly
fitted motif range and receives full weight:
:math:`w_0 = 1,\; w_N = 0`.
Similarly, for :math:`x \ge 1 - x_\mathrm{ref}`:
:math:`w_0 = 0,\; w_N = 1`.

In the central region the rescaled coordinate
:math:`u = (x - x_\mathrm{ref}) / (1 - 2x_\mathrm{ref})`
runs from 0 to 1, and the weights are:

.. math::

   w_0 = \frac{e^\eta (1-u)^\alpha}{e^\eta (1-u)^\alpha + u^\alpha},
   \qquad
   w_N = 1 - w_0.

The crossover point where :math:`w_0 = w_N = 0.5` shifts with ``eta``:

.. math::

   u_\mathrm{cross} = \frac{e^{\eta/\alpha}}{1 + e^{\eta/\alpha}}.

With ``eta = 0`` the weighting is symmetric around the mid-composition.

PME Workflow
------------

PME building proceeds in two phases.

**Phase 1 — Referencing** fits the motif interaction terms V from DFT at boundary
concentrations and produces initial energy predictions for the target level:

1. Run ``sod_comb.sh`` to create ``EQMATRIX`` and the required ``nXX/ENSEMBLE`` files.

2. Run DFT for the low-side reference levels (``n00``–``n04``) and, for PMEh, the
   symmetric high-side levels.  Place energies in ``nXX/ENERGIES`` using the
   two-column format ``m  E_nm`` (configuration index and energy in eV).

3. Run ``sod_pme.sh`` from :term:`SODPROJECT` (no ``pme.model`` needed).

``sod_pme.sh`` reads the target substitution level from ``INSOD`` and writes
PME outputs under the target folder:

- ``nXX/PME0/ENERGIES`` — low-side predictions.
- ``nXX/PME1/ENERGIES`` — high-side predictions, when high-side training data
  are available.
- ``nXX/PMEh/ENERGIES`` — hybrid predictions, when both sides are available.
- ``pme.model.tmp`` in :term:`SODPROJECT` — suggested control file with a
  bisection-selected ``calib_config_list``; rename to ``pme.model`` to use on
  future runs.
- ``nXX/pme.model`` — copy of the ``pme.model`` file used for the run, when
  present (kept as a record).

**Phase 2 — Calibrating** (optional) fits the epsilon correction terms ε from a
small set of DFT calculations at the target concentration, improving accuracy:

1. Run ``sod_gener.sh -choose <indices>`` for the indices listed in
   ``calib_config_list`` from ``pme.model.tmp``.  Compute their DFT energies and
   add them to ``nXX/ENERGIES``.

2. Rename ``pme.model.tmp`` → ``pme.model``, then re-run ``sod_pme.sh``.

After Phase 2, ``pme.model.tmp`` is rewritten with the fitted epsilon values
(``n_calib=0``) and can be renamed to ``pme.model`` for subsequent runs without
re-fitting.

pme.model Format
----------------

``SODPROJECT/pme.model`` selects the Hamiltonian variant and calibration policy for
the target level specified by ``INSOD``.  Blank lines and comment lines starting
with ``#`` are ignored.  When present, SOD copies this file to ``nXX/pme.model`` as
a record of the control data used for that target.

.. code-block:: none

   # PME Hamiltonian: 0 (PME0 low-side); 1 (PME1 high-side); 2 (PMEh hybrid)
   2

   # PMEorder cap (2, 3, or 4)
   3

   # n_calib: 0 no calibration (manual epsilons); 1..9 number of configurations for calibration
   9

   # calib_config_list: configuration indices; only first n_calib used
   25 3 20 34 8 33 12 17 7

   # epsilon_low (epsilon_0 ... epsilon_PMEorder)
   0.0 1.0 1.0 1.0

   # epsilon_high (epsilon_N ... epsilon_{N-PMEorder})
   0.0 1.0 1.0 1.0

   # PMEh_alpha (>1.0 for smooth blending, default 2.0) and PMEh_eta (default 0.0; >0 favours low-end, <0 favours high-end)
   2.0  0.0

The calibrated Hamiltonian is:

.. code-block:: none

   E(c) = V0 + ε₀ + ε₁·T₁(c) + ε₂·T₂(c) + ... + εK·TK(c)

``ε₀`` is an additive energy offset in eV (default 0.0); ``ε₁``–``εK`` are
multiplicative scale factors for the per-order motif contributions (default 1.0).

``n_calib`` controls calibration:

- ``0``: no calibration; ``epsilon``, ``alpha``, and ``eta`` are read directly
  from ``pme.model``.
- ``1`` to ``9``: fit ``epsilon`` from the first ``n_calib`` entries in
  ``calib_config_list`` using energies from ``nXX/ENERGIES``.  ``alpha`` and
  ``eta`` are always taken from ``pme.model`` (they are not fitted).  When
  ``n_calib`` < PMEorder+1, the higher-order epsilons are taken from
  ``pme.model`` and only the lower-order ones are fitted.

After a successful calibration run (Phase 2), ``pme.model.tmp`` is rewritten with the
fitted epsilon values and ``n_calib=0``.  Rename it to ``pme.model`` to lock in those
parameters for subsequent runs without re-fitting.

**Calibration scoring.**  The score reported after fitting is leave-one-out
cross-validation (LOO-CV) RMSE when there are enough configurations for each
fold to remain a determined system (``n_calib - 1 > PMEorder``).  Otherwise
it falls back to in-sample RMSE.  The output label reads ``-> LOO-CV X eV``
or ``-> X eV`` accordingly.  LOO-CV is a more honest estimate of
generalisation error than the in-sample RMSE, which is always optimistic when
the number of calibration configurations is small relative to the number of
fitted parameters.

MC Workflow
-----------

``sod_mc.sh`` is the user interface for Metropolis MC sampling.  It reads
``INMC`` and loops over ``TEMPERATURES`` by calling the underlying
single-temperature executable once per temperature.

.. note::

   ``sod_mc.sh`` runs **Metropolis** sampling only.  For energy-free uniform
   random sampling, use ``sod_random.sh`` instead (see :ref:`random-sampling`).

Example ``INMC``:

.. code-block:: none

   # Symmetry reduction (0 = off, 1 = on)
   1
   # Production steps
   24000
   # Starting configuration ('random' or space-separated site indices)
   random
   # Write trace (0 = off, 1 = write MCTRACE)
   0
   # Equilibration steps
   2400
   # Restart probability
   0.01
   # Random seed (-1 = system clock, >0 = fixed integer)
   12345

Create ``TEMPERATURES`` in :term:`SODPROJECT`, one temperature in K per line,
then run:

.. code-block:: bash

   sod_mc.sh

Metropolis sampling drives the walk with the PME Hamiltonian, so it always
requires the reference energies (``n00``/``n01``/``n02``/…).

.. note::

   **Incremental swap energies.**  Each Metropolis step swaps one occupied
   site for one hole, so ``mcsod`` updates the PME energy *incrementally*:
   only the cluster terms that contain the removed or added site change, which
   reduces the per-step cost from :math:`O(L^k)`/:math:`O(H^k)` to
   :math:`O(L^{k-1})`/:math:`O(H^{k-1})` (expansion order :math:`k`,
   :math:`L` substitutions, :math:`H = N_\text{pos} - L` holes).  The full
   recompute is still used for the starting configuration, the
   equilibration→production transition, the occasional restart move, and a
   periodic resync that bounds floating-point drift; results are numerically
   equivalent to a full recompute at every step.  The speedup grows with the
   larger of :math:`L` or :math:`H` and with the expansion order, so it is
   largest for dilute (or near-full) substitutions in large cells.  The swap
   also draws its added site directly from a maintained list of holes
   (:math:`O(1)`) rather than by rejection sampling, keeping the move proposal
   cheap and robust at any cell size and filling fraction.

MC Output
---------

The Metropolis walk is driven by the Hamiltonian, so its output lives under the
``PMEx`` variant directory together with the energies:

- ``nXX/MCT_TTTK/PMEx/ENSEMBLE``, ``ENERGIES``, ``OUTMC``, and optionally
  ``MCTRACE`` for Metropolis at integer-labelled temperature ``TTT`` K.
- ``INMC`` is copied next to the sampling output (``nXX/MCT_TTTK/PMEx/INMC``).

``ENSEMBLE`` and ``ENERGIES`` share the same configuration index ``m``.  ``OUTMC`` describes the
sampler, symmetry handling, ENSEMBLE row semantics, acceptance statistics, sample
energy statistics, epsilon values, and block-average SEM estimates.

Metropolis ``ENSEMBLE`` rows are already sampled with the PME energy bias, so
their weights are visit counts.

MC Thermodynamics
-----------------

After running Metropolis at multiple temperatures, run ``sod_mcstat.sh`` from
the target-level directory, for example:

.. code-block:: bash

   cd n12
   sod_mcstat.sh

The script reads ``../TEMPERATURES`` and all ``MCT_*K`` subdirectories.  The PME
variant is taken from ``../pme.model`` when present (otherwise it defaults to
``PMEh``, matching ``mcsod``), and ``MCT_*K/PMEx/ENSEMBLE`` + ``ENERGIES`` are
read for each temperature.  It writes ``thermodynamics.dat`` in the same style
as ``statsod`` and ``gcstatsod``.

Example
-------

``examples/example15/`` demonstrates the PME/MC workflow for Si/Ge substitution
in alpha-quartz:

.. code-block:: bash

   cd examples/example15
   sod_pme.sh
   sod_mc.sh

The committed regression references compare ``n12/PMEh/ENERGIES`` (pmesod
enumeration energies) and ``n12/MCT_300K/PMEh/OUTMC`` (Metropolis MC output).
