SOD Overview
============

SOD (standing for Site Occupancy Disorder) is a toolkit for 
modelling site-occupancy disorder in periodic solids with
the supercell ensemble method. It generates symmetry-inequivalent atomic
configurations consistent with a specified substitutional disorder pattern, and
supports the statistical-mechanical analysis of the resulting configurational
ensemble.

The code is designed to work with external atomistic calculators, allowing the
user to evaluate energies or other properties for each inequivalent
configuration and then compute ensemble averages, free energies, and related
thermodynamic quantities.

What SOD does
-------------

SOD supports a workflow in which a disordered crystalline solid is represented
by an ensemble of ordered supercell configurations. Its main tasks are:

- **Enumeration** of symmetry-inequivalent configurations under crystal symmetry
- **Degeneracy calculation** for each inequivalent configuration
- **Input generation** for external calculators (GULP, LAMMPS, VASP, CASTEP, QE)
- **Energy and property extraction** from calculated results
- **Statistical-mechanical analysis** in canonical and grand-canonical ensembles
- **Constrained Periodic Motif Expansion (CPME)** effective Hamiltonian fitting from ab initio training data
- **Monte Carlo (MC)** sampling of configuration space at finite temperature using the CPME Hamiltonian
- **Special Quasirandom Structures (SQS)** identification for optimized short-range order
- **Thermal averaging** of pair correlations (GQS)

Substitution modes
------------------

SOD supports several substitution patterns:

- **Binary**: one new species on a single site (e.g., Ni/Mg in MgO)
- **Multi-nary**: multiple species simultaneously on one site (e.g., NiCoFeCr on Ti in BCC)
- **Multi-target**: substitutions on multiple independent sites (e.g., Sr on La and Mn on Fe)
- **Multi-target multi-nary**: combined multi-nary on each of multiple sites
- **Molecules**: rigid molecular groups substituting at crystal sites (e.g., MA in perovskites)
- **Parent molecules**: rigid molecular groups that are part of the parent structure, represented by a spherical placeholder for the symmetry analysis and materialised at write time (e.g., a free-rotor cation on every A-site)
- **Vacancies**: removal of atoms (e.g., oxygen vacancies)

These modes can be combined to handle complex disorder in real materials.

Typical workflow
----------------

A standard SOD calculation proceeds as follows:

1. Prepare the input files describing the parent structure, symmetry, and
   substitutional disorder pattern (``INSOD``, ``SGO``) — :doc:`enumeration`.
2. Run ``sod_comb.sh`` to enumerate symmetry-inequivalent configurations.
3. If calculator input files are needed, ``sod_comb.sh`` invokes ``genersod``
   automatically.  To regenerate them after changing ``FILER`` or a template,
   run ``sod_gener.sh`` — :doc:`calculators`.
4. Run external calculations (GULP, LAMMPS, VASP, etc.) for each configuration,
   or evaluate them directly with MACE — :doc:`mace`.
5. Extract energies or other properties with the appropriate wrapper scripts.
6. Analyse the ensemble with ``statsod`` (canonical) or ``gcstatsod``
   (grand-canonical) — :doc:`thermodynamics`.

When full enumeration is out of reach, steps 1–2 are replaced by random
sampling (:doc:`random`) or Monte Carlo with an effective Hamiltonian
(:doc:`cpme_mc`), and a single representative structure can be selected instead
of averaging over the ensemble (:doc:`sqs_gqs`).

Main executables
----------------

The release build provides the following core executables:

- ``combsod``: enumerates inequivalent configurations and writes files such as
  ``ENSEMBLE`` and ``EQMATRIX``
- ``genersod``: generates configuration-specific input files and folder trees
- ``statsod``: performs canonical statistical-mechanical analysis, including
  Metropolis-sampled ENSEMBLE files when sampling metadata is present
- ``gcstatsod``: performs grand-canonical statistical-mechanical analysis
- ``cpmesod``: fits a periodic motif expansion Hamiltonian and evaluates energies
- ``mcsod``: explores configuration space via Metropolis Monte Carlo using the
  CPME Hamiltonian, one temperature per run
- ``mcstatsod``: thermodynamic integration over the Monte Carlo temperatures
- ``randomsod``: draws uniform random configurations without evaluating energies
- ``sqssod`` / ``gqssod``: identify Special and Generalized Quasirandom Structures
- ``invertENSEMBLE``: post-processes or transforms ``ENSEMBLE`` data
- ``peaks2spec``: converts peak data into spectra

Wrapper scripts and where to run them
-------------------------------------

The ``sod/bin/`` directory provides shell wrappers that drive the workflow.
Always use these rather than the bare executables.  Each one expects to be
called from a particular level of the project tree; those marked
**SODPROJECT/ or nXX/** detect the calling level automatically — run them from
``SODPROJECT/`` to process every substitution level at once, or from a single
``nXX/`` folder to process only that level.

.. list-table::
   :header-rows: 1
   :widths: 28 26 46

   * - Script
     - Call from
     - What it does
   * - ``sod_comb.sh``
     - ``SODPROJECT/``
     - Enumeration, and generation of calculator input files
   * - ``sod_gener.sh``
     - ``SODPROJECT/``
     - Regenerates input files from an existing ``ENSEMBLE``, skipping the
       combinatorics
   * - ``sod_*_ener.sh``, ``sod_*_cell.sh``
     - ``SODPROJECT/`` or ``nXX/``
     - Extract energies, cells and other results from calculator output
       (:doc:`calculators`)
   * - ``sod_mace.sh``
     - ``SODPROJECT/`` or ``nXX/``
     - Evaluates or relaxes the configurations with MACE (:doc:`mace`)
   * - ``sod_stat.sh``
     - ``SODPROJECT/`` or ``nXX/``
     - Canonical statistical mechanics (:doc:`thermodynamics`)
   * - ``sod_gcstat.sh``
     - ``x???/``
     - Grand-canonical statistical mechanics (:doc:`thermodynamics`)
   * - ``sod_p2s.sh``
     - folder holding ``PEAKS``
     - Broadens peak lists into spectra for ensemble averaging
   * - ``sod_random.sh``
     - ``SODPROJECT/``
     - Uniform random sampling of configuration space (:doc:`random`)
   * - ``sod_sqs.sh`` / ``sod_gqs.sh``
     - ``SODPROJECT/`` or ``nXX/``
     - SQS scoring and GQS thermal averaging (:doc:`sqs_gqs`)
   * - ``sod_cpme.sh``
     - ``SODPROJECT/``
     - Fits and evaluates the CPME Hamiltonian (:doc:`cpme_mc`)
   * - ``sod_mc.sh``
     - ``SODPROJECT/``
     - Monte Carlo sampling with the CPME Hamiltonian (:doc:`cpme_mc`)
   * - ``sod_mcstat.sh``
     - ``nXX/``
     - Thermodynamic integration over the ``MCT_*K/CPMEx/`` Monte Carlo outputs

Directory conventions
---------------------

Generated calculations follow a consistent naming convention:

- ``nXX/`` for a substitution level or composition
- ``cYY/`` for an inequivalent configuration within that composition
- ``x???/`` for grand-canonical analysis inputs

Every wrapper script relies on these names to find its inputs.  The full
layout is described in :doc:`enumeration`.

Where to go next
----------------

- :doc:`installation` — building SOD and making the scripts available
- :doc:`enumeration` — the ``INSOD`` file and the enumeration step
- :doc:`examples` — nineteen worked examples, starting with ``example01``
- :doc:`citing` — the papers to cite when you publish
