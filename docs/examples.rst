Examples
========

The SOD release includes a set of worked examples that illustrate the main
workflow, supported calculator backends, and post-processing tools. These
examples are the best starting point for learning how to prepare inputs, run
configuration generation, and analyse the resulting ensembles.

The examples are distributed in the ``examples/`` directory of the released
package.

Example layout
--------------

The examples directory contains:

- ``README``
- ``example01`` through ``example19``

The ``example01`` family demonstrates the same substitution problem with
different calculator backends. The remaining examples cover a broader range of
disorder models, statistical-mechanical workflows, and post-processing techniques.

Recommended starting point
--------------------------

New users should begin with ``example01``. This example is especially useful
because it keeps the physical problem fixed while changing only the calculator
backend. It therefore shows clearly how the SOD workflow stays the same while
the generated input files depend on the selected ``FILER`` option.

Calculator backends illustrated in ``example01`` include:

- GULP (``FILER1_gulp``)
- LAMMPS (``FILER2_lammps``)
- VASP (``FILER11_vasp``)
- CASTEP (``FILER12_castep``)
- Quantum ESPRESSO (``FILER13_QE``)
- CIF, used here for machine-learning potentials (``FILER0_mace``)

This is the best place to learn how SOD generates calculator-specific inputs.

``FILER0_mace`` writes a plain ``configuration.cif`` per configuration and is the
reference workload for :doc:`mace`, since ASE reads CIF directly.  It carries a
commented ``mace_settings.yaml``, so the whole machine-learning-potential
workflow runs with no arguments:

.. code-block:: bash

   cd examples/example01/FILER0_mace
   sod_mace.sh          # relaxes all 71 configurations, writes n04/ENERGIES
   cd n04 && sod_stat.sh

The energies MACE writes are in the same format the ``sod_*_ener.sh`` collectors
produce, so ``sod_stat.sh`` consumes them unchanged and applies the ENSEMBLE
degeneracies as usual.  Because ``FILER0_mace`` and ``FILER11_vasp`` are the same
physical problem, this is also a convenient way to compare an MLIP against DFT
on identical configurations.

Other examples
--------------

Examples 02–19 cover a wide range of disorder models and statistical workflows:

.. list-table::
   :header-rows: 1
   :widths: 8 30 62

   * - Example
     - Physical system
     - Key features
   * - 02
     - Fe/Al in magnetite
     - Binary substitution across octahedral and tetrahedral Fe sites
   * - 03
     - Fe/Sb in rutile
     - High degeneracy; enumeration focus (FILER=-1)
   * - 04
     - Al/Fe in LaFeO₃
     - Grand-canonical statistics
   * - 05
     - Zr/Sn in pyrochlore
     - Full composition range (nsubs=0:16); NMR spectra averaging
   * - 06
     - Li/Mg + H vacancy
     - **Multi-target** substitution (charge-neutral defect pair)
   * - 07
     - Fe vacancies in maghemite
     - **Vacancy** model (``%Fe`` syntax)
   * - 08
     - MA in CsPbI₃ perovskite
     - **Molecule** substitution (``@MA`` syntax)
   * - 09
     - Mg/La + O vacancy in LaFeO₃
     - **Multi-target** on two sites simultaneously
   * - 10
     - Ti₂ZrNb alloy (BCC)
     - **Multi-nary** ternary substitution
   * - 11
     - NiCoFeCr Cantor alloy (FCC)
     - **Multi-nary** quaternary substitution on primitive cell
   * - 12
     - La₀.₇₅Sr₀.₂₅Mn₀.₂₅Fe₀.₇₅O₃
     - **Multi-target** on two sites (La and Fe); committed SQS reference ranks the best multi-target SQS
   * - 13
     - La₁₋ₓSrₓFe₁₋ᵧMnᵧO₃₋ᵤ
     - **Multi-target** on three sites (La, Fe, O vacancy)
   * - 14
     - La₁₋ₓ₋ᵧSrₓBaᵧMnᵤFe₁₋ᵤO₃
     - **Multi-target multi-nary** (ternary on La + binary on Fe)
   * - 15
     - Si/Ge in α-quartz (2×2×2, 24 Si sites)
     - **CPME/MC effective Hamiltonian** fitting and finite-temperature sampling
   * - 16
     - Ni/Mg in MgO rocksalt (8 substitutions)
     - **SQS/GQS workflow** with thermal-weighted quasirandom structure selection
   * - 17
     - Al/Fe in LaFeO₃ (3×3×3, 27 Fe sites)
     - **CPMEh** third-order hybrid CPME for target level n04 with full DFT reference energies
   * - 18
     - MAPbI₃–MAPbBr₃ equimolar solid solution (4×4×4 Pm-3m, 192 halide sites)
     - **SQS via random sampling**: 50 000 draws with ``sod_random.sh``, ranked with ``sod_sqs.sh``; ``@MA`` parent-molecule syntax (all A-sites occupied by MA)
   * - 19
     - Equimolar fcc CoCrFeNi alloy (2×2×2 conventional fcc, 32 metal sites)
     - **Multinary SQS via random sampling**: 20 000 draws with ``sod_random.sh``, ranked with ``sod_sqs.sh``; selected SQS written as CIF with ``sod_gener.sh -choose``

Many examples include reference outputs such as ``ENSEMBLE``, ``ENERGIES``,
``DATA``, ``SPECTRA``, ``OUTSQS``, ``OUTGQS``, and grand-canonical ``x???/`` folders.
These are useful for understanding expected workflows and validating your own runs.

Typical way to use an example
-----------------------------

A typical workflow is:

1. Copy an example directory to a separate working location.
2. Make sure the SOD executables and wrapper scripts are available in your
   ``PATH``.
3. Run ``sod_comb.sh`` in the example directory.
4. If calculator inputs are generated, run the external calculator in the
   corresponding configuration directories.
5. Use the appropriate post-processing wrapper to extract energies or carry out
   statistical analysis.

Examples and regression testing
-------------------------------

The examples serve as a regression test suite for the released code. The test
script ``bin/sod_run_tests.sh`` validates the build against committed reference
outputs:

- **``example02–14`` combsod tests**: Regenerate all ``n*/ENSEMBLE`` files and
  compare against committed references. Validates configuration enumeration.

- **``example01`` genersod tests**: Run combsod + genersod for each calculator
  backend (GULP, LAMMPS, VASP, CASTEP, QE). Validates input file generation.

- **``example05`` statsod tests**: Extract and average energies (canonical
  ensemble). Validates statistical-mechanical analysis.

- **``example05`` gcstatsod tests**: Perform grand-canonical analysis over a
  composition range. Validates multi-level workflows.

- **``example15`` CPME/MC tests**: Fit the CPME Hamiltonian and run the
  Metropolis reduced MC workflow. Validates CPME fitting, MC output semantics,
  and the committed ``OUTMC`` reference.

- **``example15`` mcstatsod test**: Run the Metropolis MC workflow over a
  three-temperature ladder (1000/600/300 K), then perform thermodynamic
  integration with ``sod_mcstat.sh`` and diff the committed
  ``thermodynamics.dat`` reference. Validates the Gibbs-Helmholtz TI from the
  exact T→∞ reference.

- **``example15`` TI-vs-enumeration cross-check**: A physics consistency test.
  The same CPMEh Hamiltonian is evaluated two ways — exactly, via ``statsod``
  over the full 56846-configuration enumeration, and approximately, via
  ``mcsod`` Monte Carlo sampling plus ``mcstatsod`` thermodynamic integration
  over a 300–4000 K ladder. The test asserts the free energies agree to within
  a tolerance (max\|ΔF\| ≈ 2 meV observed, 5 meV tolerance), confirming that TI
  reproduces the exact statistical mechanics up to MC sampling and
  discretization error.

- **``example16`` SQS/GQS tests**: Run ``sqssod`` and ``gqssod`` on the
  8-substitution MgO rocksalt enumeration and diff against committed ``OUTSQS``
  and ``OUTGQS`` references. Validates the quasirandom-structure scoring and
  the finite-temperature Boltzmann-weighted pair-correlation averaging.

- **``example16`` statsod test**: Run ``statsod`` on the full 8043-configuration
  ``n08`` enumeration with the real GULP energies and the committed
  ``TEMPERATURES`` (0/1000/1e6 K), and diff the exact canonical
  ``thermodynamics.dat`` (E/F/S vs T). Provides the exact-energy reference that
  the GQS Boltzmann re-analysis approximates.

- **``example17`` CPMEh test**: Fit a third-order CPMEh Hamiltonian for LaFeO₃
  (target level n04) and diff the committed ``ENERGIES`` reference. Validates
  the hybrid low/high-side blending with the ``alpha``/``eta`` power-law
  weighting scheme.

See :doc:`installation` for how to run the full test suite.

Further information
-------------------

For more detailed scientific context and example-specific notes, see the
repository `README.md <https://github.com/rgraucrespo/sod/blob/master/README.md>`_ and the files included within each
example directory.
