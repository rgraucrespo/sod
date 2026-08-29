Configurational averages and thermodynamics
===========================================

Once every configuration in a level has an energy, SOD computes ensemble
averages and thermodynamic quantities, weighting each inequivalent
configuration by its degeneracy.  Two ensembles are available: **canonical**
(fixed composition, ``statsod``) and **grand-canonical** (several compositions
at a common chemical potential, ``gcstatsod``).

Canonical statistics
--------------------

Run, from ``SODPROJECT/`` to process every level, or from a single ``nXX/``
folder to process only that one:

.. code-block:: bash

   sod_stat.sh

It reads up to four input files.

ENSEMBLE
~~~~~~~~

The output of ``sod_comb.sh``, ``sod_random.sh`` or ``sod_mc.sh``, giving the
inequivalent configurations and their degeneracies.  See
:ref:`ensemble-format`.

ENERGIES
~~~~~~~~

The configuration energies, in two-column format ``m  E_nm``, where ``m`` is
the configuration index matching ``ENSEMBLE`` and ``E_nm`` is the energy in eV.
Lines starting with ``#`` are comments, and the explicit ``m`` index says which
configuration each energy belongs to, so rows need not be in order.

.. important::

   ``statsod`` requires an energy for **every** configuration of the level and
   stops with ``missing energies for N configuration(s)`` otherwise; so do
   ``gcstatsod`` and ``mcstatsod``.  The indexed format allows a partial file
   because CPME calibration reads only the configurations named in
   ``calib_config_list`` — not so that the ensemble analysis can run on a
   subset.

``ENERGIES`` is normally produced by one of the extraction scripts described in
:doc:`calculators`, or directly by ``sod_mace.sh`` (:doc:`mace`).  All the
extractors derive ``m`` from the ``cYY`` directory name and warn if energies
are missing for some configurations.

``ENTHALPIES``, written by ``sod_gulp_enth.sh`` and ``sod_vasp_enth.sh`` for
constant-pressure calculations, uses the same indexed two-column format but is
kept in a separate file so that ordinary energy-based post-processing keeps its
existing meaning.

TEMPERATURES
~~~~~~~~~~~~

A single column of temperatures (K) for the Boltzmann statistics::

   300
   600
   1000
   1250
   1500
   1750
   2000

If the file does not exist, ``statsod`` uses T = 0, 300 and 1000 K.
``gcstatsod`` has its own fallback of 300 and 1000 K.

.. note::

   If ``ENSEMBLE`` carries a Metropolis sampling-temperature header from
   ``sod_mc.sh``, ``statsod`` treats the sample as already Boltzmann-biased: it
   ignores ``TEMPERATURES``, uses the sampling temperature from ``ENSEMBLE``,
   and assigns probabilities :math:`P_m = \Omega_m / \sum_m \Omega_m`.  For
   enumerated ensembles and uniform MC output the sampling temperature is
   absent or ``-1``, and the ordinary Boltzmann form
   :math:`P_m = \Omega_m \exp(-E_m/k_\mathrm{B}T) / Z` is used.  For Metropolis
   samples, probabilities and averages are reported at the sampling
   temperature; free energy and entropy are not evaluated from the biased
   sample.

DATA
~~~~

Any additional observables to average.  The first line is the number of
columns *ncol*, followed by one row per configuration::

   2
   34.5   4.34
   37.7   4.35
   35.6   4.38
   38.8   4.41

The data can be cell lengths, volumes (see the SOD papers for strategies to
obtain average cell parameters), or any other calculated observable.  Scripts
such as ``sod_vasp_cell.sh`` help produce it — read them before use.

.. warning::

   ``CELL`` rows have **no** configuration-index column.  The VASP, GULP and
   MACE cell extractors therefore write ``CELL`` only when their results cover
   every configuration ``1..N`` in the level's ``ENSEMBLE``.  A sparse
   calculation — a CPME calibration subset, say — can still produce a valid
   indexed ``ENERGIES``, but the cell extractors warn and skip ``CELL`` rather
   than write misaligned positional data.

Output
~~~~~~

``sod_stat.sh`` writes ``probabilities.dat`` (the configuration probabilities
at each temperature) and ``thermodynamics.dat`` (the averaged quantities).  It
also writes ``ave_data.dat`` when a ``DATA`` file is present and
``ave_spectra.dat`` when ``SPECTRA`` and ``XSPEC`` are present.

The infinite-temperature limit is computed and reported as a final row whether
or not ``TEMPERATURES`` is supplied (it is skipped only for a Metropolis
sample, where free energy and entropy are not evaluated).

.. important::

   Configurational averages — cell parameters, enthalpies and similar — converge
   quickly with supercell size.  Entropies and free energies do not: they are
   not defined by averaging, converge very slowly, and are generally in large
   error at accessible supercell sizes.  We do not recommend using SOD for
   entropies and free energies unless appropriate correcting procedures are
   applied.

Grand-canonical statistics
--------------------------

Since version 0.51, statistics can be done in a grand-canonical ensemble,
combining supercells of different compositions.  See ``example04``
(perovskites) and ``example05`` (pyrochlore).

Create a folder named ``x???`` at the same level as the ``nXX/`` folders — for
example ``x250`` for composition *x* = 0.250 — and run:

.. code-block:: bash

   sod_gcstat.sh

The analysis needs ``ENSEMBLE_00``, ``ENSEMBLE_01``, … and ``ENERGIES_00``,
``ENERGIES_01``, … (one pair per substitution level), optionally ``DATA_00``,
… and optionally ``TEMPERATURES``.  If the ``nXX`` naming convention was
followed and ``x???`` sits alongside those folders, the script copies all of
these automatically — you only need to supply ``INGC``.

The INGC file
~~~~~~~~~~~~~

``INGC`` is short.  To fix the chemical potential:

.. code-block:: text

   # nsubsmin nsubsmax
   0   4
   # Specify x or mu, and provide its value
   mu -0.5
   # Stress corrections flag (lambda=0: no correction; lambda=1: bulk moduli-based correction)
   0

.. warning::

   The stress-corrections flag is not optional.  ``gcstatsod`` always reads it,
   so an ``INGC`` that stops after the ``x``/``mu`` line aborts with a Fortran
   read error.  The two parameter lines below it are read only when the flag is
   non-zero, but every committed example carries all four data lines.

Alternatively, specify the composition (the fraction *x* = nsubs/npos of
substituted sites) and let SOD find the chemical potential that reproduces it
at each temperature, by bisection on the resulting polynomial equation:

.. code-block:: text

   # nsubsmin nsubsmax
   0   4
   # Specify x or mu, and provide its value
   x 0.25
   # Stress corrections flag (lambda=0: no correction; lambda=1: bulk moduli-based correction)
   0

Here *x* = 0.25 corresponds to 2 substitutions in the canonical example, but
the grand-canonical analysis includes all compositions from 0 to 4
substitutions.  The output is ``probabilities.dat`` and ``thermodynamics.dat``,
as in the canonical case.

Stress-volume correction
~~~~~~~~~~~~~~~~~~~~~~~~

A cell with *n* substitutions (with *n* different from *xN*) that contributes
to a grand-canonical ensemble representing composition *x* carries an extra
stress-related energy cost, because the equilibrium volumes at *n/N* and at *x*
differ.  In an approximation based on the second-order Birch-Murnaghan equation
of state, this stress-volume correction is

.. math::

   E_\mathrm{SVC}(n,x) = \frac{9}{8}\, V(x)\, B(x)
   \left[ \left( \frac{V(x)}{V(n/N)} \right)^{2/3} - 1 \right]^2

where :math:`B(x)` and :math:`V(x)` are the bulk modulus and equilibrium volume
at composition *x*, and :math:`V(n/N)` the equilibrium volume at the
composition of the contributing supercell.  The expression is positive definite
for any volume departure, reduces correctly to
:math:`\tfrac{1}{2} B(x) (\Delta V)^2 / V(x)` in the small-strain limit, and
needs only :math:`B_0` as a parameter (implicitly assuming :math:`B' = 4`).

Enable it by adding to ``INGC`` (see ``example05``):

.. code-block:: text

   # Stress corrections flag (lambda=0: no correction; lambda=1: bulk moduli-based correction)
   1
   # Parameters for volume variation with x: v0, v1, bv (Angstrom^3)
   1286.687504  1302.567820  0
   # Parameters for bulk modulus variation with x: bm0, bm1, bb (GPa)
   150 150 0

These settings interpolate the equilibrium volumes and bulk moduli linearly
between the solid-solution endmembers; non-zero bowing parameters *bv* and *bb*
give a quadratic interpolation instead.

.. note::

   This functionality has not been extensively tested.  If you intend to use
   the correction scheme, please contact the SOD developers first.

Averaging spectra
-----------------

Both the canonical and the grand-canonical analysis can average spectra over
the configurational ensemble.

What is calculated for a configuration is often a list of peaks rather than a
spectrum.  In that case run

.. code-block:: bash

   sod_p2s.sh

to broaden the peaks into the spectra that will be averaged.  It reads two
files with fixed names:

``PEAKS``
   One line per configuration, each holding that configuration's list of peaks.

``INP2S``
   The parameters of the broadening.  From ``example05/n04``, which together
   with that level's ``PEAKS`` regenerates its committed ``SPECTRA`` and
   ``XSPEC`` exactly:

   .. code-block:: text

      # nconfs
      22
      # peaks
      12
      # xmin
      -680
      # xmax
      -620
      # npoints
      601
      # broadening (sigma)
      0.85

It writes ``SPECTRA`` (a count line, then one line of intensities per
configuration on the *x* grid) and ``XSPEC`` (the *x* values themselves).
``PEAKS`` and ``INP2S`` are shipped for levels n00, n01, n02, n04, n13, n14 and
n15 of ``example05``; the remaining levels ship only the resulting spectra.

If these files are present in the ``nXX/`` folders when ``sod_stat.sh`` or
``sod_gcstat.sh`` runs, a file ``ave_spectra.dat`` is written with the
configurational averages of the spectra at each temperature.

Where to go next
----------------

- :doc:`cpme_mc` — reaching compositions too large to enumerate
- :doc:`sqs_gqs` — single representative structures instead of full averages
