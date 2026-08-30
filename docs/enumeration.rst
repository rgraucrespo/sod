Enumerating the configurations
==============================

Enumeration is the core step of SOD: given a parent structure, its symmetry
operators, and a substitution pattern, ``combsod`` identifies every
symmetry-inequivalent configuration of the supercell and computes its
degeneracy.  Everything else in SOD — calculator input generation, statistics,
SQS/GQS scoring, CPME fitting — consumes the output of this step.

Setting up a project
--------------------

Create one working directory per project (a family of compositions within a
given supercell).  This directory is called ``SODPROJECT`` throughout the
documentation.  It must contain:

- ``INSOD`` — the parent structure and the substitution pattern
- ``SGO`` — the symmetry operators of the parent space group
- a template file for the target calculator, if ``FILER > 0``
  (see :doc:`calculators`)

Then run::

   sod_comb.sh

Symmetry data without enumeration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

   sod_comb.sh -eqmatrix-only

stops after the symmetry analysis, writing ``supercell.cif``, ``EQMATRIX`` and
``OPERATORS`` but no ``ENSEMBLE``, and exits successfully.

Use it when the supercell is too large to enumerate but the symmetry data is
wanted on its own.  ``sod_random.sh -sym on`` folds its draws onto symmetry
orbits and so requires an ``EQMATRIX``, and a plain ``sod_comb.sh`` cannot
supply one at those sizes: on example18, C(192,96) is about 3x10^17
configurations, and combsod writes the three files in under a second before
aborting in the allocator with ``problem too large for SOD``.  The data was
therefore obtainable only as the by-product of a run that failed.  The files
are identical either way -- the flag changes only where combsod stops -- and
the ``INSOD`` substitution counts are still validated, so a level that exceeds
the available sites is still rejected.

.. _insod:

The INSOD file
--------------

``INSOD`` contains all the information needed by the combinatorics program.
**The format is rigid**: keep the comment lines and the blank lines exactly as
they appear in the examples, and edit only the data lines.  The simplest way to
start is to copy the ``INSOD`` from an example and modify it.

A complete file, from ``example01/FILER1_gulp``:

.. code-block:: text

   # Title
   Ni/Mg substitutions in MgO rocksalt

   # a,b,c,alpha,beta,gamma
   4.22  4.22  4.22  90.000  90.000  90.000

   # nsp: Number of species in the parent structure
   2

   # symbol(1:nsp): Atom symbols
   Mg O

   # natsp0(1:nsp): Number of atoms for each species
   1 1

   # coords0(1:nat0,1:3): Coordinates of each atom (one line per atom)
   0.0  0.0  0.0
   0.5  0.5  0.5

   # na,nb,nc (supercell definition)
   2 2 2

   # sptarget: Species to be substituted
   1

   # nsubs: Substitution counts. One line per target site.
   4

   # newsymbol: One line per target site.
   Ni Mg

   # FILER
   1

Field by field:

``Title``
   Free-text label, carried into the generated files.

``a,b,c,alpha,beta,gamma``
   Cell parameters of the **parent** cell, in ångström and degrees.

``nsp`` / ``symbol`` / ``natsp0`` / ``coords0``
   The parent structure: number of species, their symbols, the number of atoms
   given for each species, and their fractional coordinates (one atom per
   line).  It is enough to give the asymmetric unit — the ``SGO`` operators
   generate the rest — but a full cell is also accepted.  A species symbol
   prefixed with ``@NAME`` is a :ref:`parent molecule <parent-molecules>`.

``na,nb,nc``
   Supercell multipliers along **a**, **b** and **c**.

``sptarget``
   Index (or indices) of the parent species to substitute, referring to the
   order in ``symbol``.  Up to 5 target sites may be listed, separated by
   spaces; each one gets its own line in ``nsubs`` and ``newsymbol``.

``nsubs``
   Substitution counts — see :ref:`nsubs-formats` below.

``newsymbol``
   One line per target site.  Each line lists the *k* new-species symbols
   followed by the symbol kept for the remaining sites (*k*\ +1 symbols in
   total).  The prefixes ``@`` and ``%`` select molecules and vacancies; see
   :ref:`molecules-vacancies`.

``FILER``
   Which calculator input files to generate.  See :doc:`calculators`.

.. _nsubs-formats:

Substitution counts (nsubs)
---------------------------

The ``nsubs`` field controls how many atoms of each new species are placed on
the target site(s).  Several formats are supported:

**Fixed count** (e.g. ``4``)
   Enumerate configurations with exactly 4 substitutions.  The endpoints are
   allowed: ``0`` gives the single unsubstituted parent configuration, and a
   count equal to the number of target sites gives the single fully
   substituted one.

**Range** (e.g. ``1:8``)
   Loop over all integer values from 1 to 8, creating one ``nXX/`` folder per
   concentration.  The range may start at ``0`` (e.g. ``0:16`` spans parent to
   full substitution).  Only valid with a single target site and a single new
   species.

**Multi-nary** (e.g. ``1 2`` or ``2 2 2``)
   Place the given numbers of each new species simultaneously on one target
   site.  The file format supports up to 5 new species per target, but
   **full enumeration with** ``combsod`` **is limited to 3 new species per
   target** (quaternary disorder): its combinatorial ranking algorithm does not
   generalise further, and full enumeration is infeasible for quinary and
   senary compositions in any case.  For 4 or more new species use
   :doc:`random` sampling with :doc:`sqs_gqs` scoring instead — see
   ``example19``.  ``randomsod``, ``sqssod`` and ``genersod`` all support the
   full 5-species range.

**Multi-target** (two lines, e.g. ``2`` then ``1``)
   One line per target site, in the order listed in ``sptarget``.  This example
   places 2 substitutions on the first target and 1 on the second, enumerating
   all joint configurations simultaneously under the full crystal symmetry.

**Multi-target multi-nary** (e.g. line 1 ``1 1``, line 2 ``1``)
   Each target site independently carries several new species.  Each line holds
   the space-separated counts for the new species on that site.

.. _molecules-vacancies:

Molecules (@NAME) and vacancies (%NAME)
---------------------------------------

Two prefixes extend ``newsymbol`` beyond simple atomic substitution:

``@NAME`` — molecule
   SOD reads ``NAME.xyz`` from the working directory (standard XYZ format:
   number of atoms, comment line, then ``symbol x y z`` per line in ångström),
   computes its centre of mass, and places the molecule at the substituted site
   with a uniformly random orientation.  Each site gets an independent
   rotation.  All output formats expand the molecule into individual atoms.

``%Symbol`` — vacancy
   The atom at the substituted site is omitted from all output files.
   ``Symbol`` is informational only (e.g. ``%Fe`` or ``%FeB``).

Molecules and vacancies may appear in any slot of any target site — including
the remaining-species slot, and in multi-nary and multi-target substitutions —
and up to 10 distinct molecule types may be used at once.  The field accepts up
to 5 characters (prefix plus a 4-character name), e.g. ``@MA``, ``@FA``,
``@CO2``, ``%FeB``.

.. _parent-molecules:

Parent-structure molecules
~~~~~~~~~~~~~~~~~~~~~~~~~~

A molecule can also belong to the **parent structure** rather than only being
substituted in.  This is the right representation when a molecular group is
dynamically isotropic (a free rotor): the parent site is a single spherical
placeholder species that obeys the ``SGO`` site symmetry and takes part in the
enumeration, and the molecule is materialised only when ``genersod`` writes the
explicit calculation inputs.

Write the parent species in the ``symbol`` list with an ``@NAME`` prefix — the
same convention used for substituted molecules.  To make every A-site of the
parent a methylammonium molecule (``MA.xyz``):

.. code-block:: text

   # symbol(1:nsp): Atom symbols
   @MA Pb I

SOD strips the ``@``, treats the site as an ordinary spherical placeholder
species ``MA`` for the symmetry analysis and enumeration, and records
``MA.xyz`` as the molecule to materialise.  The placeholder must **not** also be
a substitution target: a site cannot be both a disorder target and a parent
molecule.

.. note::

   Parent molecules are expanded by ``genersod`` only.  The ``cpmesod`` and
   ``mcsod`` energy-model programs read the ``@NAME`` symbol as a plain
   placeholder species and do not materialise parent molecules.

The SGO file
------------

``SGO`` holds the matrix–vector representations of the symmetry operators of
the parent space group.  Each line gives one row of an operator: the first
three numbers are a row of the rotation matrix and the fourth is the
corresponding component of the translation vector.

First check whether your space group is in the library shipped with SOD, and if
so copy it in under the name ``SGO``:

.. code-block:: bash

   cp ROOTSOD/sod(version)/sgo/225.sgo ./SGO

Otherwise build the file from the International Tables of Crystallography, or
from the `Bilbao Crystallographic Server <https://www.cryst.ehu.es>`_.

Output of sod_comb.sh
---------------------

``sod_comb.sh`` prints the total number of configurations and the number of
inequivalent ones to standard output, and writes:

``ENSEMBLE``
   The inequivalent configurations, one line each: configuration index,
   degeneracy, then the substitution sites.  See :ref:`ensemble-format`.

``EQMATRIX``
   How each supercell operator transforms each atom position.

``nXX/cYY/``
   One folder per substitution level, and inside it one folder per
   inequivalent configuration, holding that configuration's calculator input
   (:doc:`calculators`).  Nothing is generated when ``FILER = -1``.

.. _ensemble-format:

The ENSEMBLE file
~~~~~~~~~~~~~~~~~

``ENSEMBLE`` is written in **format version 3**.  Its first line reads::

   Enumerated ensemble: 71 configurations; sum_degeneracies = 35960

or ``Uniform random ensemble: ...`` / ``Metropolis ensemble (300.0 K): ...``
for sampled ensembles.  It is followed by one target line per substituted
sublattice, a column-header comment, and then the configuration rows.

This is the only format SOD reads.  Version 2 files written before release 0.80
(a ``#`` comment header followed by ``N substitutions in M sites``) are no
longer accepted; regenerate the level with ``sod_comb.sh`` or
``sod_random.sh``.

.. note::

   ``ENSEMBLE`` was called ``OUTSOD`` in versions before 0.80.

Directory layout
----------------

A project ends up with the following hierarchy::

   SODPROJECT/          INSOD, SGO, template files, ENSEMBLE, EQMATRIX
     n01/               one folder per substitution level (zero-padded)
       ENSEMBLE
       ENERGIES
       c01/  c02/ ...   one folder per inequivalent configuration
     n02/  n04/ ...
     x250/              grand-canonical working folder (user-created)
       INGC

The ``nXX``/``cYY`` convention is what every downstream script relies on to
find its inputs.

Where to go next
----------------

- :doc:`calculators` — generating and running the calculation inputs
- :doc:`thermodynamics` — configurational averages from the computed energies
- :doc:`random` — sampling instead of enumerating, when the level is too large
