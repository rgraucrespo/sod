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
- **Periodic Motif Expansion (PME)** effective Hamiltonian fitting from ab initio training data
- **Monte Carlo (MC)** sampling of configuration space at finite temperature using the PME Hamiltonian
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
   substitutional disorder pattern (``INSOD``, ``SGO``).
2. Run ``sod_comb.sh`` to enumerate symmetry-inequivalent configurations.
3. If calculator input files are needed, ``sod_comb.sh`` invokes ``genersod``
   automatically.  To regenerate them after changing ``FILER`` or a template,
   run ``sod_gener.sh``.
4. Run external calculations (GULP, LAMMPS, VASP, etc.) for each configuration.
5. Extract energies or other properties with the appropriate wrapper scripts.
6. Analyse the ensemble with ``statsod`` (canonical) or ``gcstatsod``
   (grand-canonical), or find Special Quasirandom Structures with ``sqssod``.

Main executables
----------------

The release build provides the following core executables:

- ``combsod``: enumerates inequivalent configurations and writes files such as
  ``ENSEMBLE`` and ``EQMATRIX``
- ``genersod``: generates configuration-specific input files and folder trees
- ``statsod``: performs canonical statistical-mechanical analysis, including
  Metropolis-sampled ENSEMBLE files when sampling metadata is present
- ``gcstatsod``: performs grand-canonical statistical-mechanical analysis
- ``pmesod``: fits a periodic motif expansion Hamiltonian and evaluates energies
- ``mcsod``: explores configuration space via Metropolis or Uniform Monte Carlo
  using the PME Hamiltonian
- ``sqssod`` / ``gqssod``: identify Special and Generalized Quasirandom Structures
- ``invertENSEMBLE``: post-processes or transforms ``ENSEMBLE`` data
- ``peaks2spec``: converts peak data into spectra

Wrapper scripts
---------------

The ``sod/bin/`` directory also provides shell wrappers that simplify the
standard workflow. These include:

- ``sod_comb.sh`` for configuration generation and workflow setup
- ``sod_stat.sh`` for canonical statistics
- ``sod_gcstat.sh`` for grand-canonical statistics
- ``sod_pme.sh`` for periodic motif expansion (PME) fitting and evaluation
- ``sod_mc.sh`` for Monte Carlo sampling using the PME Hamiltonian
- ``sod_sqs.sh`` / ``sod_gqs.sh`` for SQS/GQS structure identification
- calculator-specific extraction scripts such as
  ``sod_gulp_ener.sh``, ``sod_vasp_ener.sh``, ``sod_castep_ener.sh``, and
  ``sod_qe_ener.sh``

Supported calculators
---------------------

SOD can generate inputs for several external calculators, depending on the
``FILER`` value in ``INSOD``:

+---------+----------------------+----------------------------------+
| FILER   | Output mode          | Template file                    |
+=========+======================+==================================+
| ``-1``  | No calculator files  | None                             |
+---------+----------------------+----------------------------------+
| ``0``   | CIF per config       | None                             |
+---------+----------------------+----------------------------------+
| ``1``   | GULP                 | ``template_input.gin``           |
+---------+----------------------+----------------------------------+
| ``2``   | LAMMPS               | ``template_in.lammps``           |
+---------+----------------------+----------------------------------+
| ``11``  | VASP                 | None; ``POSCAR`` is generated    |
+---------+----------------------+----------------------------------+
| ``12``  | CASTEP               | ``template_castep.cell``         |
+---------+----------------------+----------------------------------+
| ``13``  | Quantum ESPRESSO     | ``template_pw.in``               |
+---------+----------------------+----------------------------------+

Directory conventions
---------------------

Generated calculations follow a consistent naming convention:

- ``nXX/`` for a substitution level or composition
- ``cYY/`` for an inequivalent configuration within that composition
- ``x???/`` for grand-canonical analysis inputs

These conventions are used by the wrapper scripts during post-processing and
analysis.

Source tree
-----------

In the released package:

- ``sod/src/`` contains the Fortran source code
- ``sod/bin/`` contains compiled binaries and shell wrappers
- ``sod/sgo/`` contains the space-group operator library
- ``sod/examples/`` contains worked examples

Building SOD
------------

The release build is controlled by ``sod/Makefile``. To compile from the
repository root, run:

.. code-block:: bash

   make all

More information
----------------

The main SOD paper describes the underlying methodology and should be cited if
you use SOD in your research:

   Grau-Crespo, R., Hamad, S., Catlow, C. R. A., & De Leeuw, N. H. (2007).
   *Symmetry-adapted configurational modelling of fractional site occupancy in
   solids.* **Journal of Physics: Condensed Matter**, 19(25), 256201.

The main ``README.md`` in the repository provides:

- Comprehensive documentation of all input/output file formats
- Detailed descriptions of each of the 14 worked examples
- Guides for each supported calculator (GULP, LAMMPS, VASP, CASTEP, QE)
- Extended workflow examples (PME, MC, SQS, GQS, etc.)

For questions or issues, contact the SOD developers.

Where to go next
----------------

For detailed installation instructions, see :doc:`installation`. For worked examples,
see :doc:`examples`.
