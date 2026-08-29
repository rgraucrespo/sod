Calculator input files
======================

Once the inequivalent configurations are known, SOD writes a complete input
file set for the calculator of your choice into each ``nXX/cYY/`` folder.  The
``FILER`` value in :ref:`INSOD <insod>` selects the calculator; a single
template file in the working directory carries everything that is *not*
configuration-specific.

Supported calculators
---------------------

.. list-table::
   :header-rows: 1
   :widths: 10 34 30 26

   * - FILER
     - Output
     - Template file
     - Generated in each ``cYY/``
   * - ``-1``
     - No calculator files
     - —
     - —
   * - ``0``
     - CIF (P1, one per configuration)
     - —
     - ``configuration.cif``
   * - ``1``
     - GULP
     - ``template_input.gin``
     - ``input.gin``
   * - ``2``
     - LAMMPS
     - ``template_in.lammps``
     - ``in.lammps`` + ``conf.data``
   * - ``11``
     - VASP
     - — (``POSCAR`` generated directly)
     - ``POSCAR``
   * - ``12``
     - CASTEP
     - ``template_castep.cell``
     - ``castep.cell``
   * - ``13``
     - Quantum ESPRESSO
     - ``template_pw.in``
     - ``pw.in``

``FILER = 0`` is also the mode used by the MACE calculator, which reads the
per-configuration CIF files directly — see :doc:`mace`.

To regenerate the input files after changing ``FILER`` or editing a template,
run ``sod_gener.sh``: it reuses the existing ``ENSEMBLE`` and skips the
combinatorics.

Template file philosophy
------------------------

For every calculator with ``FILER > 0`` except VASP, you provide one template
file in the working directory.  It should look as much as possible like the
real input file that will later be run, with two SOD-specific additions:

``@configuration_structure@``
   A token on a line by itself, marking where SOD inserts the
   configuration-specific structure block (cell parameters and atomic
   positions).  For GULP, CASTEP and Quantum ESPRESSO it is mandatory and must
   appear exactly once, alone on its line.  **LAMMPS is the exception**: the
   structure is written to the separate ``conf.data``, so ``template_in.lammps``
   neither needs nor uses this token.

``@configuration_number@``
   An optional token that may appear anywhere — in titles or output file names,
   for instance — and is replaced by the zero-padded configuration index.

All calculator-specific settings — force-field parameters, k-points,
convergence thresholds, output directives — stay in the template unchanged.
SOD only inserts the structure.

Templates follow the naming convention ``template_<real_output_name>``, so the
generated file name is immediately obvious from the template name.

Calculator-specific notes
-------------------------

GULP (``template_input.gin``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SOD replaces ``@configuration_structure@`` with the ``cell``, ``frac`` and
atom-type block.  The template must define the force field in one of three
ways: inline, through a ``library`` directive, or through a ``kim_model``
directive (an OpenKIM potential, as in ``example11``); genersod stops if none
is present.  If the template contains a ``library`` directive, SOD reads
the library, extracts only the entries needed for the atom types in the
problem, and writes a reduced copy into each configuration folder, so every
folder is self-contained.  If the force field is defined inline in the
template, no library file is copied.

If a SOD species label differs from the corresponding GULP atom type, add a
mapping line to the template::

   # sod_type_map <SOD_species> <GULP_type>

LAMMPS (``template_in.lammps``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``template_in.lammps`` is a LAMMPS run script containing ``read_data
conf.data`` and all force-field settings.  SOD copies it into each ``cYY/``
folder with only two modifications: ``# sod_type_map`` lines are stripped (they
are SOD directives, not valid LAMMPS syntax), and ``@configuration_number@`` is
substituted if present.  Everything else — coefficients, masses, charges, run
commands — is copied verbatim.  SOD also writes the configuration-specific
``conf.data``.

To build ``conf.data``, SOD must know which LAMMPS numeric type corresponds to
each SOD species; **it cannot infer this**.  Provide explicit mapping comment
lines::

   # sod_type_map <SOD_species> <role> type=<N> [bond_type=<M>]

where ``<role>`` is ``core`` or ``shell``.  One ``core`` line is required for
every SOD species; a ``shell`` line (with ``bond_type=<M>``) is additionally
required for any species represented as a core–shell pair.  SOD stops with an
error if a species has no mapping.

VASP
~~~~

No template is needed: SOD generates ``POSCAR`` in each ``cYY/`` folder.  If an
``INCAR`` sits alongside ``INSOD`` it is copied into every configuration
folder; if it is absent SOD warns and still generates the inputs.  You must
supply ``KPOINTS`` and ``POTCAR`` (and ``INCAR``, if not provided at the
``INSOD`` level) yourself.

CASTEP (``template_castep.cell``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A normal ``.cell`` file.  SOD replaces ``@configuration_structure@`` with the
``%BLOCK LATTICE_CART`` and ``%BLOCK POSITIONS_FRAC`` blocks.  The ``.param``
file must be supplied separately; SOD does not manage it.

Quantum ESPRESSO (``template_pw.in``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A normal ``pw.x`` input.  SOD assumes ``ibrav = 0`` and replaces
``@configuration_structure@`` with the ``CELL_PARAMETERS`` and
``ATOMIC_POSITIONS`` blocks.  ``ATOMIC_SPECIES`` is defined by you in the
template and is not modified; pseudopotential files are not copied.

MACE (``FILER = 0``)
~~~~~~~~~~~~~~~~~~~~

No template, and no external calculator to run.  SOD writes
``configuration.cif`` into each ``cYY/`` folder and ``sod_mace.sh`` evaluates
those structures directly, writing ``nXX/ENERGIES`` itself.  See :doc:`mace`.

Running the calculations: job_sender
------------------------------------

After ``sod_comb.sh`` (or ``sod_gener.sh``) has generated the input files, a
``job_sender`` script is written to the working directory.  Running it executes
the external calculator in each ``nXX/cYY/`` folder in turn:

.. code-block:: bash

   ./job_sender

No ``job_sender`` is created for ``FILER = 0``, where there is no external
calculator to run.

By default ``job_sender`` calls the bare command name for each calculator.  If
your installation uses a different executable name or a wrapper, set the
matching environment variable in your ``~/.bashrc`` — ``job_sender`` inherits
it from your interactive shell:

.. list-table::
   :header-rows: 1
   :widths: 25 25 50

   * - Variable
     - Default
     - Calculator
   * - ``SOD_VASP``
     - ``vasp``
     - VASP
   * - ``SOD_GULP``
     - ``gulp``
     - GULP
   * - ``SOD_LAMMPS``
     - ``lmp``
     - LAMMPS
   * - ``SOD_CASTEP``
     - ``castep``
     - CASTEP
   * - ``SOD_QE``
     - ``pw.x``
     - Quantum ESPRESSO

For example, if VASP is installed as ``vasp_std``:

.. code-block:: bash

   export SOD_VASP=vasp_std

Collecting the results
----------------------

Each calculator has a wrapper script that walks the ``nXX/cYY/`` tree and
writes the two-column ``nXX/ENERGIES`` file consumed by the analysis tools:

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Script
     - What it extracts
   * - ``sod_gulp_ener.sh``
     - Final energies from GULP ``output.gout``
   * - ``sod_gulp_enth.sh``
     - Enthalpies from GULP output (constant pressure) → ``ENTHALPIES``
   * - ``sod_gulp_free.sh``
     - Vibrational free energies from GULP output
   * - ``sod_gulp_single_ener.sh``
     - Single-point energies from GULP output
   * - ``sod_gulp_cell.sh``
     - Cell parameters from GULP output → ``CELL``
   * - ``sod_vasp_ener.sh``
     - Final energies from VASP ``OUTCAR``
   * - ``sod_vasp_enth.sh``
     - Enthalpies from VASP ``OUTCAR`` (``energy(sigma->0) + P V``)
   * - ``sod_vasp_cell.sh``
     - Cell parameters from VASP ``CONTCAR`` → ``CELL``
   * - ``sod_vasp_mag.sh``
     - Magnetic moments from VASP ``OUTCAR``
   * - ``sod_castep_ener.sh``
     - Energies from CASTEP ``castep.castep``
   * - ``sod_qe_ener.sh``
     - Energies from QE ``pw.out`` (converted Ry → eV)

All of them can be run from ``SODPROJECT/`` to process every level at once, or
from a single ``nXX/`` folder to process just that level.  With MACE no
extraction step is needed: ``sod_mace.sh`` writes ``ENERGIES`` directly.

The energy file formats, and what the analysis tools do with them, are
described in :doc:`thermodynamics`.
