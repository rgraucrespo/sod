Installation
============

Requirements
------------

SOD is distributed as a source archive and built with GNU Make and a modern
Fortran compiler.

The minimum requirements are:

- **GNU-compatible ``make``** (available on all Unix-like systems)
- **Fortran 2003 compiler** (``gfortran``, ``ifort``, ``flang``, or equivalent)
  — tested with ``gfortran`` 9.x and later

Older Fortran compilers may not support the allocatable array and derived-type
features used in SOD.

Obtaining the source
--------------------

Download the SOD release archive, for example ``sod0.80.tar.gz``, and copy it
to a directory that will contain the installation. Let this parent directory be
called ``ROOTSOD``.

Unpack the archive with:

.. code-block:: bash

   tar xzvf sod(version).tar.gz

This creates a versioned SOD directory inside ``ROOTSOD``, for example
``ROOTSOD/sod0.80/``.

Building SOD
------------

Change into the unpacked SOD directory and compile all executables with:

.. code-block:: bash

   cd ROOTSOD/sod(version)
   make all

This builds the executables in the ``bin`` directory.

To install executables and scripts system-wide (default prefix ``/usr/local``):

.. code-block:: bash

   make install              # installs to /usr/local/bin
   make install PREFIX=/opt  # installs to /opt/bin

To remove build products, run:

.. code-block:: bash

   cd ROOTSOD/sod(version)
   make clean

Making the executables available
--------------------------------

Add the ``bin`` directory to your shell ``PATH`` so that the SOD executables
and wrapper scripts can be called from any working directory.

For example:

.. code-block:: bash

   export PATH=$PATH:ROOTSOD/sod(version)/bin

To make this persistent, add that line to your ``.bashrc`` or equivalent shell
startup file.

Configuring external calculators
---------------------------------

``sod_comb.sh`` and ``sod_gener.sh`` write a ``job_sender`` script that invokes
the external calculator in each configuration directory.  By default
``job_sender`` calls the bare executable name (``vasp``, ``gulp``, ``lmp``,
``castep``, or ``pw.x``).  To use a different name or wrapper, export the
relevant variable in your ``~/.bashrc``:

.. code-block:: bash

   export SOD_VASP=vasp_std      # VASP
   export SOD_GULP=gulp6         # GULP
   export SOD_LAMMPS=lmp_mpi     # LAMMPS
   export SOD_CASTEP=castep19    # CASTEP
   export SOD_QE=pw.x            # Quantum ESPRESSO

The variable is inherited by ``job_sender`` from your interactive shell; no
sourcing or alias tricks are needed.  The default bare command is used when the
variable is unset.

Main programs and scripts
-------------------------

The build provides the main compiled executables together with shell wrappers in
``bin/``. Typical entry points (always use the ``sod_*.sh`` wrappers, not the
bare executables directly):

- ``sod_comb.sh`` — configuration enumeration and input-file generation
- ``sod_stat.sh`` — canonical statistical analysis
- ``sod_gcstat.sh`` — grand-canonical statistical analysis
- ``sod_pme.sh`` — PME Hamiltonian fitting and evaluation
- ``sod_mc.sh`` — Monte Carlo sampling using the PME Hamiltonian
- ``sod_mcstat.sh`` — thermodynamic integration over MC temperatures (run from ``nXX/PMEx/``)
- ``sod_sqs.sh`` / ``sod_gqs.sh`` — SQS/GQS quasirandom structure identification
- ``sod_gener.sh`` — regenerate calculator input files after changing ``FILER`` or a template

Verifying the installation
--------------------------

After compilation, you can verify that the executables are visible with a
command such as:

.. code-block:: bash

   which sod_comb.sh

If the path is set correctly, this should point to the SOD ``bin`` directory.

Running the regression tests
-----------------------------

Verify the build is correct by running the regression test suite from the
top-level SOD directory:

.. code-block:: bash

   make test
   # or equivalently:
   ./bin/sod_run_tests.sh

This runs 28 tests (combsod, genersod, statsod, pmesod, mcsod, mcstatsod,
gcstatsod, sqssod, and gqssod workflows) against committed reference outputs.
All tests should pass.

Troubleshooting
---------------

**Compilation fails: "gfortran: command not found"**
   Install the GNU Fortran compiler. On macOS with Homebrew: ``brew install gcc``.
   On Linux (Debian/Ubuntu): ``apt-get install gfortran``.

**Compilation fails with "Error: Unrecognized option"**
   Your compiler is too old. SOD requires Fortran 2003 or later. Update your
   compiler or use a different one (e.g., ``ifort``).

**Scripts not found after adding ``bin/`` to ``PATH``**
   Verify your ``.bashrc`` or shell startup file was edited correctly, and run
   ``source ~/.bashrc`` (or equivalent) to reload your shell environment.
   Then test with ``which sod_comb.sh``.

**Regression tests fail**
   Most often due to floating-point tolerance issues on different platforms.
   Check the test output for specific failures. Contact the SOD developers
   if failures persist.

Next steps
----------

For an introduction to the code and workflow, see :doc:`overview`.
For worked examples, see :doc:`examples`.
