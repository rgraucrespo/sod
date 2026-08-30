SOD Documentation
=================

**SOD** (Site Occupancy Disorder) is a toolkit for modelling site-occupancy
disorder in periodic solids using the supercell ensemble method. It enumerates
symmetry-inequivalent atomic configurations, generates calculator-specific
input files, and provides statistical-mechanical analysis tools.

**Source code on GitHub:** `github.com/rgraucrespo/sod
<https://github.com/rgraucrespo/sod>`_ — releases, examples, and issue tracker.

Key Features
------------

- **Symmetry-aware enumeration** of inequivalent configurations in supercells
- **Multi-nary and multi-target** substitutions on multiple sites
- **Flexible input generation** for GULP, LAMMPS, VASP, CASTEP, and Quantum ESPRESSO
- **Statistical analysis** in canonical and grand-canonical ensembles
- **Constrained Periodic Motif Expansion (CPME)** effective Hamiltonian fitted from ab initio training data
- **Monte Carlo sampling** of configuration space at finite temperature using the CPME Hamiltonian
- **Special Quasirandom Structures (SQS)** and generalized GQS methods
- **Machine-learning potentials** via MACE, with GPU-batched energies and relaxation
- **Comprehensive examples** covering diverse disorder models

Quick Start
-----------

New users should:

1. Read the :doc:`overview` for a high-level introduction
2. Follow :doc:`installation` to build the code
3. Work through :doc:`enumeration` — the core of every SOD calculation
4. Explore :doc:`examples` — start with ``example01``

How this documentation is organised
-----------------------------------

**Getting started** introduces the method and gets SOD built on your machine.
**Running SOD** follows the standard workflow in order: enumerate the
configurations, generate and run the calculations, then analyse the ensemble.
**Alternatives to enumeration** covers what to do when enumerating the
configurations and averaging over all of them is out of reach: sampling,
quasirandom structures, and effective Hamiltonians. **Reference** holds the
worked examples, terminology and citations.

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   overview
   installation

.. toctree::
   :maxdepth: 2
   :caption: Running SOD

   enumeration
   calculators
   mace
   thermodynamics

.. toctree::
   :maxdepth: 2
   :caption: Alternatives to enumeration

   random
   sqs_gqs
   cpme_mc

.. toctree::
   :maxdepth: 2
   :caption: Reference

   examples
   glossary
   citing

.. toctree::
   :maxdepth: 1
   :caption: Project links

   GitHub repository <https://github.com/rgraucrespo/sod>
