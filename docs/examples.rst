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
- ``example01`` through ``example14``

The ``example01`` family demonstrates the same substitution problem with
different calculator backends. The remaining examples cover a broader range of
disorder models and statistical-mechanical workflows.

Recommended starting point
--------------------------

New users should begin with ``example01``. This example is especially useful
because it keeps the physical problem fixed while changing only the calculator
backend. It therefore shows clearly how the SOD workflow stays the same while
the generated input files depend on the selected ``FILER`` option.

Calculator backends illustrated in ``example01`` include:

- GULP
- LAMMPS
- VASP
- CASTEP
- Quantum ESPRESSO

This is the best place to learn how SOD generates calculator-specific inputs.

Other examples
--------------

The later examples illustrate a broader range of use cases, including:

- binary substitutions
- multi-nary substitutions
- multi-species substitutions
- vacancy models
- molecule substitutions
- canonical statistics workflows
- grand-canonical statistics workflows
- SPBE extrapolation workflows

Many examples already include reference outputs such as ``OUTSOD``,
``ENERGIES``, ``DATA``, ``SPECTRA``, or grand-canonical ``x???/`` folders.
These are useful both for understanding the expected workflow and for checking
that a run has behaved as intended.

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

The examples are also used as a regression test set for the released code.

In particular:

- ``example02`` through ``example14`` are used in ``combsod`` regression tests
- ``example01`` variants are used to test ``genersod`` output generation
- selected ``nXX/`` and ``x???/`` folders are used to test ``statsod`` and
  ``gcstatsod`` workflows

Further information
-------------------

For more detailed scientific context and example-specific notes, see the main
``README.md`` file and the files included within each example directory.
