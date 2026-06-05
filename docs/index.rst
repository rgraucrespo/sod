SOD Documentation
=================

**SOD** (Site Occupancy Disorder) is a toolkit for modelling site-occupancy
disorder in periodic solids using the supercell ensemble method. It enumerates
symmetry-inequivalent atomic configurations, generates calculator-specific
input files, and provides statistical-mechanical analysis tools.

Key Features
------------

- **Symmetry-aware enumeration** of inequivalent configurations in supercells
- **Multi-nary and multi-target** substitutions on multiple sites
- **Flexible input generation** for GULP, LAMMPS, VASP, CASTEP, and Quantum ESPRESSO
- **Statistical analysis** in canonical and grand-canonical ensembles
- **Periodic Motif Expansion (PME)** effective Hamiltonian fitted from ab initio training data
- **Monte Carlo sampling** of configuration space at finite temperature using the PME Hamiltonian
- **Special Quasirandom Structures (SQS)** and generalized GQS methods
- **Comprehensive examples** covering diverse disorder models

Quick Start
-----------

New users should:

1. Read the :doc:`overview` for a high-level introduction
2. Follow :doc:`installation` to build the code
3. Explore :doc:`examples` — start with ``example01``

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   overview
   installation
   examples
   sqs_gqs
   pme_mc
   glossary
