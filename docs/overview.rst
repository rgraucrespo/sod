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

- enumeration of symmetry-inequivalent configurations
- calculation of configurational degeneracies
- generation of calculator-specific input files for each configuration
- post-processing of calculated energies or properties
- statistical-mechanical analysis in canonical and grand-canonical ensembles

Typical workflow
----------------

A standard SOD calculation proceeds as follows:

1. Prepare the input files describing the parent structure, symmetry, and
   substitutional disorder.
2. Run ``sod_comb.sh`` to enumerate symmetry-inequivalent configurations.
3. If required, let ``genersod`` generate calculator-specific input files and
   directory trees for each configuration.
4. Run external calculations for the generated configurations.
5. Extract energies or other quantities with the appropriate wrapper scripts.
6. Analyse the ensemble with ``statsod``, ``gcstatsod``, or SPBE-based tools,
   depending on the problem.

Main executables
----------------

The release build provides the following core executables:

- ``combsod``: enumerates inequivalent configurations and writes files such as
  ``OUTSOD`` and ``EQMATRIX``
- ``genersod``: generates configuration-specific input files and folder trees
- ``statsod``: performs canonical statistical-mechanical analysis
- ``gcstatsod``: performs grand-canonical statistical-mechanical analysis
- ``spbesod``: performs simple point-based energy extrapolation
- ``invertOUTSOD``: post-processes or transforms ``OUTSOD`` data
- ``peaks2spec``: converts peak data into spectra

Wrapper scripts
---------------

The ``sod/bin/`` directory also provides shell wrappers that simplify the
standard workflow. These include:

- ``sod_comb.sh`` for configuration generation and workflow setup
- ``sod_stat.sh`` for canonical statistics
- ``sod_gcstat.sh`` for grand-canonical statistics
- ``sod_spbe0.sh`` and ``sod_spbe1.sh`` for SPBE workflows
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
- ``spbe0/`` and ``spbe1/`` for SPBE extrapolation branches

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

Where to go next
----------------

For detailed installation instructions, see :doc:`installation`. For worked examples,
see :doc:`examples`.
