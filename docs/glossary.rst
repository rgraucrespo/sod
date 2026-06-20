Glossary
========

.. glossary::
   :sorted:

   calibration energies
      Energies of selected configurations at the *target* intermediate
      composition ``nXX``, used to fit the ε scale factors of a :term:`PME`
      Hamiltonian.  They are stored in ``nXX/ENERGIES`` (two-column format
      ``m  E_nm``: configuration index and energy in eV) and are read at the
      row indices listed in the ``calib_config_list`` line of
      :term:`pme.model`.  The number of indices actually consumed is set by
      the ``n_calib`` field (0–9).  Calibration energies are distinct from
      :term:`reference energies`: reference energies define the Hamiltonian's
      V terms; calibration energies tune the ε corrections.

   binary substitution
      The simplest substitution mode: a single new species replaces a
      fraction of atoms on one target site type. For example, Ni/Mg in
      rocksalt MgO. Binary substitutions are a special case of
      :term:`multi-nary substitution` with one new species.

   canonical ensemble
      A statistical-mechanical treatment of a system at fixed composition.
      Each :term:`inequivalent configuration` is weighted by its Boltzmann
      factor :math:`g_i \exp(-E_i / k_\mathrm{B} T)`, where :math:`g_i`
      is the :term:`degeneracy` and :math:`E_i` is the energy. The
      :term:`statsod` program carries out canonical averaging and writes
      ``probabilities.dat`` and ``thermodynamics.dat``. See also
      :term:`grand-canonical ensemble`.

   chemical potential
      The thermodynamic variable that controls composition in a
      :term:`grand-canonical ensemble` calculation. ``gcstatsod`` either
      accepts a fixed chemical potential (*mu* in :term:`INGC`) or solves
      for the value that reproduces a target composition *x* at each
      temperature.

   combsod
      The enumeration executable. Given a parent structure (:term:`INSOD`)
      and symmetry operators (:term:`SGO`), it identifies all
      :term:`inequivalent configuration`\s, computes their
      :term:`degeneracy`, and writes :term:`ENSEMBLE` and :term:`EQMATRIX`.
      Normally invoked via ``sod_comb.sh``.

   configuration
      An ordered arrangement of atoms in the :term:`supercell` consistent
      with the specified substitution pattern. Each distinct arrangement is
      a configuration; SOD identifies the subset of
      :term:`inequivalent configuration`\s from the full combinatorial set
      and assigns a :term:`degeneracy` to each.

   degeneracy
      The number of symmetry-equivalent supercell arrangements that map
      onto a given :term:`inequivalent configuration` under the crystal
      symmetry operators. Degeneracies act as statistical weights in
      thermodynamic averaging: a configuration with degeneracy :math:`g`
      contributes :math:`g` times as much to the partition function as a
      non-degenerate one. Degeneracies are listed in :term:`ENSEMBLE`.

   EQMATRIX
      An output file written by :term:`combsod` that records, for each
      symmetry operator, how it maps atomic positions in the
      :term:`supercell` onto one another. EQMATRIX is used internally by
      :term:`genersod` and the :term:`SQS` tools and does not normally
      require direct inspection.

   FILER
      An integer parameter in :term:`INSOD` that selects which external
      calculator's input files :term:`genersod` produces. The supported
      values are: ``-1`` (no files), ``0`` (CIF), ``1`` (GULP), ``2``
      (LAMMPS), ``11`` (VASP), ``12`` (CASTEP), and ``13`` (Quantum
      ESPRESSO). For most calculators a template file must be present in
      the :term:`SODPROJECT`; VASP is the exception and generates a
      ``POSCAR`` directly.

   gcstatsod
      The grand-canonical statistical analysis executable. Reads
      :term:`ENSEMBLE` and ``ENERGIES`` files for multiple compositions
      together with :term:`INGC`, and produces Boltzmann-weighted averages
      of energies and other observables across the composition range. Normally
      invoked via ``sod_gcstat.sh``. See :term:`grand-canonical ensemble`.

   genersod
      The input-file generation executable. It reads :term:`ENSEMBLE` and a
      calculator template file to produce configuration-specific input files
      in the ``nXX/cYY/`` directory tree (see :term:`SODPROJECT`). Normally
      invoked automatically by ``sod_comb.sh``; can be re-run independently
      via ``sod_gener.sh`` after changing :term:`FILER` or the template.

   grand-canonical ensemble
      A statistical-mechanical treatment that includes configurations from
      supercells with different numbers of substitutions (i.e. different
      compositions). The :term:`chemical potential` or target composition is
      held fixed rather than the substitution count. The grand-canonical
      approach is necessary when comparing energies across compositions or
      when computing phase boundaries. Implemented in :term:`gcstatsod`.
      See also :term:`canonical ensemble`.

   GQS
   Generalized Quasirandom Structures
      An extension of :term:`SQS` to finite temperatures. Rather than
      selecting the single configuration whose :term:`Warren parameters`
      are closest to ideal random values at 0 K, GQS computes the
      thermally weighted average of pair correlations from the full
      Boltzmann ensemble at each temperature. The result reflects how
      short-range order evolves with temperature. Implemented in
      ``gqssod`` and invoked via ``sod_gqs.sh``.

   inequivalent configuration
      A :term:`configuration` that is not related to any other enumerated
      configuration by a crystal symmetry operation. The set of all
      inequivalent configurations is the minimal representative set: each
      is counted once, weighted by its :term:`degeneracy`. SOD's primary
      task is to identify this set efficiently for arbitrary substitution
      patterns and space groups.

   INMC
      Input file for Metropolis Monte Carlo sampling, placed in
      :term:`SODPROJECT`.  Specifies: symmetry reduction flag (``0`` off /
      ``1`` on), number of production steps (``n_prod``), starting
      configuration (``random`` or space-separated site indices), whether to
      write a per-step trace (``write_trace``), number of equilibration steps
      (``n_equil``), restart probability (``restart_prob``), and random seed
      (``-1`` for system clock, any positive integer for a fixed seed).
      Temperature is not stored in ``INMC``; instead, ``TEMPERATURES`` in
      :term:`SODPROJECT` lists one temperature per line for the ``sod_mc.sh``
      loop.  Read by :term:`mcsod` via ``sod_mc.sh``.  A copy is saved next to
      each MC run's output (``nXX/MCT_TTTK/PMEx/``) for reference.  Uniform
      random sampling uses :term:`randomsod` and does not read ``INMC``.

   INGC
      Input file for grand-canonical analysis, placed in an ``x???/``
      working folder. It specifies the composition range
      (``nsubsmin``/``nsubsmax``) and whether to fix the
      :term:`chemical potential` (``mu``) or the composition fraction
      (``x``). Optionally enables the :term:`stress-volume correction`.

   pme.model
      Input/control file for the :term:`PME` Hamiltonian at a given
      composition level.  The editable control file is placed at
      ``SODPROJECT/pme.model`` and applies to the target level specified by
      :term:`INSOD`; SOD copies it to ``nXX/pme.model`` as a run record.
      The file has 7 data lines:

      - **PME choice** — which Hamiltonian variant to use: ``0`` (PME0,
        low-side only), ``1`` (PME1, high-side only), or ``2`` (PMEh,
        weighted hybrid).
      - **PMEorder** — cap on the expansion order (must be 2, 3, or 4); a
        lower effective order is used with a warning if training data is
        insufficient.
      - **n_calib** — how many :term:`calibration energies` to use when
        fitting the ε corrections: ``0`` (no calibration, use manual ε
        values given in this file) or ``1``–``9`` (fit ε₀..ε_{n−1} from the
        first ``n_calib`` indices of ``calib_config_list``).
      - **calib_config_list** — space-separated configuration indices (up
        to 9) selected from ``nXX/ENERGIES`` for calibration.  Auto-filled
        by ``pmesod`` via recursive bisection of the predicted energy range
        when writing ``pme.model.tmp``.
      - **epsilon_low** — ε₀..ε_{PMEorder} for the low-side expansion
        (``PMEorder + 1`` values; defaults 0, 1, 1, ...; ε₀ is an additive
        energy offset in eV, ε₁..ε_K are multiplicative scale factors).
      - **epsilon_high** — ε₀..ε_{PMEorder} for the high-side expansion
        (same format as ``epsilon_low``).
      - **alpha eta** — PMEh hybrid weighting parameters (``alpha`` > 1
        sharpness of the low/high blend, default ``2.0``; ``eta``
        log-scale asymmetry, default ``0.0``; ``> 0`` favours low-end,
        ``< 0`` favours high-end).  Always present but only used for PMEh.

      Optionally followed by ``# PME energy terms`` and the fitted V₀ and
      V_k coefficients (low and high side), which are written by ``pmesod``
      after fitting and are read back on subsequent runs to skip refitting.

      When ``SODPROJECT/pme.model`` is absent, ``pmesod`` defaults to PMEh
      (or PME0 if no high-side reference data is available) with ε = 1,
      ``alpha = 2.0``, ``eta = 0.0``, and writes
      ``SODPROJECT/pme.model.tmp`` plus a copy at ``nXX/pme.model.tmp`` as a
      starting template.

   INSOD
      The main input file for a SOD calculation, placed in the
      :term:`SODPROJECT`. It specifies the parent crystal structure,
      supercell dimensions, symmetry file reference, substitution targets,
      the :term:`nsubs` count(s), and the :term:`FILER` value. The format
      is strict: each data value occupies a fixed line, preceded by a blank
      line and a comment line.

   INSQS
      Input file for :term:`SQS` and :term:`GQS` analysis, placed in the
      :term:`SODPROJECT` or a specific ``nXX/`` folder. Specifies the
      maximum cluster order, cutoff radii, cluster weights, and the
      scoring mode (least-squares or van de Walle).

   MAINFOLDER
      Former name for :term:`SODPROJECT`. Retained here for reference;
      all documentation and scripts now use SODPROJECT.

   MCT
      A Metropolis Monte Carlo output directory, named ``MCT_TTTK`` where
      ``TTT`` is the integer temperature in kelvin (e.g. ``MCT_600K``).
      Lives directly under ``nXX/`` (sampling method first); because the
      Metropolis walk is Hamiltonian-driven, the run's ``ENSEMBLE``,
      ``ENERGIES``, ``OUTMC``, ``INMC`` and optional ``MCTRACE`` are written
      in a ``PMEx`` subdirectory, ``nXX/MCT_TTTK/PMEx/``.  All ``MCT_*K``
      directories are processed together by ``sod_mcstat.sh`` to perform
      thermodynamic integration.

   mcsod
      The Metropolis Monte Carlo sampling executable.  Uses a :term:`PME`
      effective Hamiltonian to drive sampling of the configuration space at a
      single temperature.  The loop over multiple temperatures is handled by
      ``sod_mc.sh``, which invokes one single-temperature MC run per
      ``TEMPERATURES`` entry.  Reads :term:`INSOD`, :term:`INMC`, and
      ``SODPROJECT/pme.model`` when present, and writes output to
      ``nXX/MCT_TTTK/PMEx/``.  Any ε corrections from
      :term:`calibration energies` are already baked into the loaded
      Hamiltonian before sampling begins.  Normally invoked via
      ``sod_mc.sh``.  For energy-free uniform random sampling, see
      :term:`randomsod`.

   randomsod
      The uniform random-sampling executable (wrapper ``sod_random.sh``).
      Draws ``nconfigs`` independent uniform configurations at the target
      level with **no** energy evaluation and writes them to
      ``nXX/random/ENSEMBLE`` (with ``-symmetry on``, the degeneracy column
      holds visit counts).  It is the sampling counterpart of :term:`combsod`
      for levels too large to enumerate; energies, if wanted, are computed a
      posteriori via the usual structure-writer → DFT → ``statsod`` path.
      Reads :term:`INSOD` and ``SGO`` (always), plus :term:`EQMATRIX` with
      ``-symmetry on``.

   mcstatsod
      The Monte Carlo thermodynamics program.  Run from ``nXX/``, it reads
      ``../TEMPERATURES`` and the ``MCT_TTTK/PMEx/ENSEMBLE`` and
      ``MCT_TTTK/PMEx/ENERGIES`` files for each sampled temperature (the
      ``PMEx`` variant is taken from ``../pme.model``, defaulting to
      ``PMEh``), computes
      :math:`E_\mathrm{ave}(T) = \sum \omega E / \sum \omega` from MC
      visit counts, and integrates :math:`d(\beta F)/d\beta = U(\beta)`
      by trapezoidal quadrature.  The exact high-temperature reference
      :math:`S(T \to \infty) = k_\mathrm{B} \ln C(n_\mathrm{pos},
      \mathrm{lev})` anchors the integration; a linear extrapolation
      handles the tail from the highest sampled inverse temperature to
      :math:`\beta = 0`.  Output is written to ``thermodynamics.dat`` in
      the same column format as :term:`statsod`.  Always invoked via
      ``sod_mcstat.sh`` from ``nXX/``.  See also :term:`MCT`.

   molecule substitution
      Substitution of a rigid multi-atom molecular group at a crystal site,
      specified in :term:`INSOD` using the ``@NAME`` prefix in
      ``newsymbol``. SOD reads the geometry from ``NAME.xyz`` (standard XYZ
      format), computes the centre of mass, and places the molecule at the
      substituted site with an independent random orientation for each
      configuration. All calculator output formats expand the molecule into
      its constituent atoms. See also :term:`parent molecule` and
      :term:`vacancy`.

   parent molecule
      A molecular group that is part of the parent structure rather than
      substituted in. The parent keeps a single spherical placeholder species
      at the site — a point that obeys the :term:`SGO` site symmetry and the
      inequivalent-configuration enumeration — and ``genersod`` materialises
      the molecule only when writing the calculation inputs. Declared by
      writing the parent species with an ``@NAME`` prefix in the :term:`INSOD`
      ``symbol`` list (e.g. ``@MA``) — the same ``@NAME`` convention used for
      :term:`molecule substitution`. The ``@`` is stripped, leaving an
      ordinary placeholder species ``NAME``, and ``NAME.xyz`` is recorded as
      the molecule to materialise. The placeholder must not also be a
      substitution target. Each site is expanded with an independent random
      orientation, as for :term:`molecule substitution`; this suits a
      dynamically isotropic (free-rotor) group, represented as a sphere for
      the symmetry analysis and made explicit only for individual
      calculations. Expanded by ``genersod`` across all calculator output
      formats; ``combsod`` treats ``@NAME`` as a plain placeholder species.

   multi-nary substitution
      Simultaneous placement of two or more new species on a single target
      site type, covering ternary, quaternary, and higher disordered alloys.
      The ``nsubs`` field in :term:`INSOD` takes a space-separated list of
      counts, one per new species (e.g. ``2 2 2`` for three species). Can
      be combined with :term:`multi-target substitution`.

   multi-target substitution
      Simultaneous substitution on two or more crystallographically
      distinct site types, with all configurations enumerated jointly under
      the full crystal symmetry. The ``sptarget`` line in :term:`INSOD`
      lists the involved site indices, and ``nsubs`` provides one count
      per target site. Can be combined with :term:`multi-nary substitution`.

   nsubs
      The substitution count field in :term:`INSOD`. Controls how many
      atoms of each new species are placed on the target site(s). Accepts a
      fixed integer, a range (``n1:n2``) for scanning all compositions, a
      space-separated list for :term:`multi-nary substitution`, or multiple
      lines for :term:`multi-target substitution`. The endpoints are valid:
      ``0`` yields the single parent (unsubstituted) configuration, and a
      count equal to the number of target sites yields the single
      fully-substituted configuration.

   OUTGQS
      Output file written by ``gqssod`` (via ``sod_gqs.sh``). Contains the
      thermally averaged :term:`Warren parameters` at each temperature
      requested in the ``TEMPERATURES`` file, computed by Boltzmann
      weighting over the configurational ensemble. See :term:`GQS`.

   ENSEMBLE
      The main configuration-list file. For enumeration by :term:`combsod`,
      it contains one line per :term:`inequivalent configuration` listing its
      index, :term:`degeneracy`, and substituted atom positions. For Monte
      Carlo output from ``mcsod``, rows and Omega have sampler-dependent
      meanings documented in ``OUTMC`` and the PME/MC guide. ENSEMBLE is
      required by :term:`statsod`, :term:`gcstatsod`, and
      ``sqssod``. Known as ``OUTSOD`` in versions of SOD before 0.80.

   OUTSQS
      Output file written by ``sqssod`` (via ``sod_sqs.sh``). A ranked
      list of :term:`inequivalent configuration`\s, from best to worst
      :term:`SQS` candidate, with their :term:`Warren parameters` and
      deviation scores. Rank 0 is the configuration closest to ideal random
      mixing. See :term:`SQS`.

   pair correlations
      See :term:`Warren parameters`.

   PME
   Periodic Motif Expansion
      An effective Hamiltonian for SOD that expresses the energy of a
      configuration as a sum of interaction terms up to 4-body, fitted
      from :term:`reference energies` at low and/or high substitution
      levels.  Three variants are available:

      - **PME0** — low-side expansion only; energy evaluated as
        :math:`E_0^\mathrm{low} + \sum_i \varepsilon_i V_i^\mathrm{low}`.
      - **PME1** — high-side (hole) expansion only; energy evaluated
        symmetrically in terms of the unoccupied sites.
      - **PMEh** — weighted hybrid of PME0 and PME1,
        :math:`w_\mathrm{low}(x) E_\mathrm{low} + w_\mathrm{high}(x) E_\mathrm{high}`,
        where :math:`x` is the substitution fraction and the weights
        follow a piecewise power-law scheme controlled by ``alpha``
        (sharpness, default 2.0) and ``eta`` (log-scale asymmetry,
        default 0.0).  Edge regions within ``PMEorder/N`` of either
        boundary use the corresponding end-member PME with full weight.

      The expansion order and ε scale factors are controlled via
      :term:`pme.model`.  PME energies for all enumerated configurations at
      the target level are written to ``nXX/PMEx/ENERGIES``.  For large
      target levels where full enumeration is intractable, the PME
      Hamiltonian drives Monte Carlo sampling (see :term:`mcsod`).

   pmesod
      The PME fitting and evaluation executable.  Reads :term:`INSOD`, the
      reference energy files (``n00/ENERGIES``, ``n01/ENERGIES``, …), and
      optionally :term:`pme.model` (with any sparse calibration energies in
      ``nXX/ENERGIES`` referenced by ``calib_config_list``).  Fits the
      :term:`PME` V terms from :term:`reference energies`, applies any
      :term:`calibration energies` via the ε corrections, evaluates the
      resulting Hamiltonian for all enumerated configurations at the target
      level, and writes predictions under ``nXX/PME0/``, ``nXX/PME1/``,
      and/or ``nXX/PMEh/``.  When no ``pme.model`` is present it also writes
      ``SODPROJECT/pme.model.tmp`` as a suggested control file, with a
      bisection-selected ``calib_config_list`` spanning the predicted
      energy range.  Normally invoked via ``sod_pme.sh``.

   reference energies
      Energies of low-x (n00–n04) or high-x (n(M)–n(M-4)) configurations
      used to fit the V1–V4 interaction terms that define a :term:`PME`
      Hamiltonian.  They are stored in the respective ``nXX/ENERGIES`` files
      in two-column format ``m  E_nm`` (configuration index and energy in eV),
      and are the *training data* for the Hamiltonian.  Contrast with
      :term:`calibration energies`, which tune the fitted model at an
      intermediate target composition.

   SGO
      The space group operators file, placed in :term:`SODPROJECT`.
      Each line encodes one symmetry operation as three matrix elements
      (one row of the rotation matrix) and one translation component.
      Pre-computed SGO files for common space groups are distributed with
      SOD in the ``sgo/`` library; custom files can be constructed from the
      International Tables of Crystallography or the Bilbao Crystallographic
      Server.

   SODPROJECT
      The top-level working directory for a SOD project (previously called
      MAINFOLDER).  It contains :term:`INSOD`, :term:`SGO`, any calculator
      template files, ``TEMPERATURES``, ``INMC``, and receives
      :term:`EQMATRIX` and ``supercell.cif`` after running
      :term:`combsod`.  Composition-level subdirectories ``nXX/`` hold
      enumeration results, calculator runs, and PME outputs; ``x???/``
      directories hold grand-canonical inputs.

   site-occupancy disorder
      Structural disorder in which two or more distinct atomic species
      share the same crystallographic site type — that is, the atoms are
      mixed randomly (or with partial order) over a set of symmetry-equivalent
      positions. SOD models this type of disorder by explicitly enumerating
      the ordered configurations that represent the disordered solid.

   sod_type_map
      A comment directive used in calculator template files (GULP and
      LAMMPS) to map SOD species names to the type names expected by the
      calculator. Lines of the form ``# sod_type_map <SOD_species>
      <calc_type>`` are parsed by :term:`genersod` and stripped from the
      generated input files; they do not appear in the final calculator
      inputs. Required whenever a SOD species label differs from the
      corresponding calculator atom type.

   SQS
   Special Quasirandom Structures
      The :term:`inequivalent configuration`\s whose :term:`Warren
      parameters` most closely match the ideal random-mixing values at the
      same composition. Identified by ``sqssod`` (via ``sod_sqs.sh``), which
      ranks all configurations by their deviation from ideal pair
      correlations. The rank-0 configuration is the best SQS for use in
      single-configuration calculations that aim to represent the disordered
      solid. See also :term:`GQS`.

   statsod
      The canonical statistical analysis executable. Reads :term:`ENSEMBLE`,
      ``ENERGIES``, an optional ``DATA`` file, and an optional
      ``TEMPERATURES`` file. Enumeration and Uniform Monte Carlo ENSEMBLE files
      are interpreted with Boltzmann weighting. Metropolis-sampled ENSEMBLE
      files are identified from their sampling-temperature header and are
      averaged with ``Omega/sum(Omega)`` because the sample already contains
      the energy bias. Normally invoked via ``sod_stat.sh``. See
      :term:`canonical ensemble`.

   stress-volume correction
      A correction applied in :term:`grand-canonical ensemble` analysis to
      account for the elastic energy penalty that arises when a supercell
      at one composition contributes to the ensemble at a different
      composition. The correction is based on the second-order
      Birch-Murnaghan equation of state and requires only the bulk modulus
      and equilibrium volumes of the two solid-solution end-members as
      input. Enabled in :term:`INGC`.

   supercell
      A periodic cell formed by repeating the primitive unit cell of the
      parent structure in one or more directions. SOD enumerates all
      :term:`inequivalent configuration`\s within a fixed supercell; the
      supercell size controls the accessible compositions and the quality of
      configurational averages and :term:`SQS`.

   supercell ensemble method
      The approach used by SOD to model substitutional disorder: all
      :term:`inequivalent configuration`\s of a fixed :term:`supercell`
      are explicitly enumerated, and thermodynamic properties are computed
      as Boltzmann-weighted averages over this ensemble. As supercell size
      increases, the ensemble converges toward the thermodynamic limit.

   vacancy
      A lattice site from which the host atom is absent. Vacancies are
      specified in :term:`INSOD` using the ``%NAME`` prefix in
      ``newsymbol`` (e.g. ``%O`` for an oxygen vacancy); the atom is
      simply omitted from all generated output files. Multiple vacancy
      types can coexist with ordinary and :term:`molecule substitution`
      modes in the same calculation.

   Warren parameters
      Normalized pair correlation functions that measure the degree of
      chemical short-range order at each neighbour distance. For a binary
      alloy at composition *x*, the Warren parameter for a given pair
      distance equals *x* under ideal random mixing. SOD's :term:`SQS`
      and :term:`GQS` tools rank configurations by how closely their
      Warren parameters match these ideal random values.
