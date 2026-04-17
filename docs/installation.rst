Installation
============

Requirements
------------

SOD is distributed as a source archive and built with GNU Make and a Fortran
compiler.

The minimum requirements are:

- GNU-compatible ``make``
- ``gfortran`` or another compatible Fortran compiler

Obtaining the source
--------------------

Download the SOD release archive, for example ``sod0.71.tar.gz``, and copy it
to a directory that will contain the installation. Let this parent directory be
called ``ROOTSOD``.

Unpack the archive with:

.. code-block:: bash

   tar xzvf sod(version).tar.gz

This creates a versioned SOD directory inside ``ROOTSOD``, for example
``ROOTSOD/sod0.71/``.

Building SOD
------------

Change into the unpacked SOD directory and compile all executables with:

.. code-block:: bash

   cd ROOTSOD/sod(version)
   make all

This builds the executables in the ``bin`` directory.

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

Main programs and scripts
-------------------------

The build provides the main compiled executables together with shell wrappers in
``bin/``. Typical entry points include:

- ``sod_comb.sh`` for configuration generation
- ``sod_stat.sh`` for canonical statistical analysis
- ``sod_gcstat.sh`` for grand-canonical statistical analysis
- ``combsod``
- ``genersod``
- ``statsod``
- ``gcstatsod``

Verifying the installation
--------------------------

After compilation, you can verify that the executables are visible with a
command such as:

.. code-block:: bash

   which sod_comb.sh

If the path is set correctly, this should point to the SOD ``bin`` directory.

If the release includes the regression test script, it can be run from the
top-level SOD directory with:

.. code-block:: bash

   ./run_tests.sh

Next steps
----------

For an introduction to the code and workflow, see :doc:`overview`.
For worked examples, see :doc:`examples`.
