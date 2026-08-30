## Version 0.92 (30 August 2026)

Structure generation gains a way to say *which* configurations you want, symmetry
data can be produced for supercells too large to enumerate, and three defects
found in a pre-release audit are fixed.

- **`genersod`: `-choose` now takes a selection label as well as configuration indices.** The SQS workflow used to end with "read `OUTSQS`, then `sod_gener.sh -choose <index>`" — a manual step whose answer is fixed, and whose file invites one specific silent mistake: `OUTSQS` begins with a **Rank** column followed by a **Config** column, and rank 1 is almost never configuration 1, so `-choose 1` generates a perfectly valid structure that is not the SQS. `-choose` now also accepts `bestSQS`, `lowestENERGY`, `highestENERGY` or `random`, each optionally followed by how many to take:

  ```
  sod_gener.sh -choose bestSQS          # the best-ranked configuration in OUTSQS
  sod_gener.sh -choose bestSQS 10       # the ten best, in rank order
  sod_gener.sh -choose lowestENERGY 5   # the five lowest-energy configurations
  sod_gener.sh -choose random 20        # twenty at random, all different
  ```

  Labels are case-insensitive, and each reads its file (`OUTSQS`, `ENERGIES`, `ENSEMBLE`) from the same directory as the ensemble being generated, so the answer always belongs to that ensemble — including when it lives in `nXX/random/`. `random` draws with probability proportional to degeneracy, because the rows of `ENSEMBLE` stand for `degen` microstates each: drawing them equiprobably would over-represent the rare ones and a property averaged over the subset would not estimate the ensemble average. The draw is without replacement, so the selected configurations are always all different. There is deliberately no `bestGQS`: `gqssod` prints its ranking to the terminal only, and ranks once per temperature, so the label could not name a single structure.

- **`sod_comb.sh -eqmatrix-only`: symmetry data without enumeration.** `randomsod -sym on` needs an `EQMATRIX`, but the only thing that produced one was the enumeration that is infeasible at exactly the supercell sizes where random sampling is the point. The files were in fact already obtainable — on a 4×4×4 halide ensemble of ~3×10¹⁷ configurations, `combsod` writes `supercell.cif`, `EQMATRIX` and `OPERATORS` in about a second and only then aborts in the allocator — but only as the by-product of a run that fails. This makes it a supported operation that exits successfully. The `INSOD` substitution counts are still validated, so an impossible level is still rejected.

- **Structure generation through a symlinked project directory silently produced nothing.** `SODPROJECT` was canonicalised but `LAUNCH_DIR` was not, so on any path containing a symlink — macOS `mktemp` under `/var`, or any symlinked project folder — the level name failed to resolve, the resulting absolute path was truncated to the 40 characters of an internal buffer, the level was reported as "skipped", and `sod_gener.sh` exited 0 having written no structures and then wrote a `job_sender` for them. All three faults are fixed. Reported and fixed by Antonio Rivas.

- **`-choose lowestENERGY`/`highestENERGY` could overrun a buffer.** An `ENERGIES` file with more usable lines than the ensemble has configurations — plausible, since it is written by the `sod_*_ener.sh` scripts and edited by hand — wrote past the end of an array. Energies are now stored by configuration index, and partial or repeated files are reported rather than mishandled.

- **INSOD rejected long molecule names instead of truncating them in silence.** `@WATER` became `WATE` and surfaced much later as `Error: WATE.xyz not found`, a filename the user never wrote. Both the target-site and parent-molecule paths now fail at once, naming the token and the limit.

- **`sod_sqs.sh` no longer reports success having scored nothing**, and distinguishes an ensemble it could not find from one it found but could not score — the two need different advice.

- **New `example20`**: molecules and vacancies together on multi-nary target sites — a 2×2×2 (MA,FA)Pb(I,Br)₃ perovskite with Schottky vacancies, whose A site is `@FA %VA @MA` and halide site `Br %VX I`. This is the case 0.91 enabled but nothing exercised.

- **Everyone in `CONTRIBUTORS.md` is now credited on the public repository.** Releases are squashed into a single commit, so only the person running the release appeared as a contributor; the release commit now carries a `Co-authored-by:` trailer for each contributor, which GitHub's Contributors sidebar counts. `AUTHORS.md` becomes `CONTRIBUTORS.md`, naming the project lead separately from the people who have contributed code.

## Version 0.91 (29 August 2026)

Correctness and documentation release. One silent structure-generation bug, a
restructured manual, and the input files for the spectra workflow put back.

- **`genersod`: structure files were wrong for three or more target sites** (breaking bug fix): with `sptarget` naming three or more sublattices and one new species on each, the atoms of the third and later targets were written with the *second* target's symbols. Driving example13 (Sr on La, Mn on Fe, a vacancy on O) at `FILER = 0` produced a CIF containing 30 Fe, 7 La, 2 Mn and 1 Sr — every oxygen relabelled, no oxygen left — with no warning of any kind. Enumeration, degeneracies and `ENSEMBLE` were never affected: this was purely the writers, so only runs that generated calculator input for three or more targets are involved. Anyone who has run such a case should regenerate the input files with `sod_gener.sh` and check the composition of a configuration before trusting earlier results. Two related gaps are fixed with it: the LAMMPS writer did not expand a molecule (`@NAME`) placed on a second target, writing its symbol truncated to three characters, and molecules and vacancies were written literally (`%O` as an "atom") in multi-target multi-nary substitutions. All six writers now resolve target sites through one shared routine that handles any number of targets and species, so molecules and vacancies also work with multi-nary substitution for the first time.

- **`cpmesod` wrote an `ENERGIES` file that nothing could read**: the CPME predictions in `nXX/CPME0/`, `nXX/CPME1/` and `nXX/CPMEh/` came out as a bare column of energies, while `statsod`, `gcstatsod` and `mcstatsod` all require the indexed two-column `m  E` format documented for `ENERGIES` (and already written by `mcsod`). Running `sod_stat.sh` on a CPME level failed with `Error: could not open or read ENERGIES` — the very step the CPME workflow exists to feed. The file is now written in the standard format; the predicted energies themselves are unchanged.

- **`sod_mace.sh -writerelaxed no`**: relax without keeping a structure file per configuration. `-relax` always wrote `cYY/relaxed.cif` (or `relaxed_POSCAR`), which is right for thousands of configurations and wrong for millions — at 10⁶ configurations that is 16.4 GiB and 3M inodes for 9.7 GB of data. The option keeps `ENERGIES`, `MACE_RELAXATION.dat` and `ENTHALPIES`/`CELL` and writes no per-configuration structure. It is a storage measure, not a speed one: writing the structures is about 1 % of a fixed-cell relaxation.

- **Documentation reorganised, and ten incorrect claims corrected**: the readthedocs manual is now grouped into Getting started / Running SOD / Methods / Reference, with new pages for enumeration (the `INSOD` file field by field), calculator input generation, and configurational averages and thermodynamics — material that previously existed only in `README.md`. The README follows the same order and has a table of contents. Several documented recipes did not work as written: the minimal `INGC` shown for grand-canonical analysis omitted the stress-corrections flag, which `gcstatsod` always reads, so following it produced a Fortran runtime error; the `INMC` listing showed an "MC sampler" first field that `mcsod` does not read, which shifted every later value by one line; `ENERGIES` was described as accepting a subset of configurations, but `statsod`, `gcstatsod` and `mcstatsod` all require every configuration and abort otherwise (the indexed format exists for CPME calibration); `@configuration_structure@` was said to be mandatory for every calculator, though LAMMPS neither needs nor uses it; and `sod_mcstat.sh` is run from `nXX/`, not `nXX/CPMEx/`. Two example descriptions had stale numbers (example02 has 97 inequivalent configurations at n04, not 99; example06 is a 3×3×3 supercell with 24 inequivalent configurations, not 2×2×2 with 9).

- **example05: `PEAKS` and `INP2S` restored**: the spectra-averaging workflow (`sod_p2s.sh`) had no input files anywhere in the package — only the `SPECTRA`/`XSPEC` they produce — so the one worked example of that step could not be run. Both files are back for levels n00, n01, n02, n04, n13, n14 and n15, verified by regenerating each level's committed spectra from them byte for byte.

- **The test suite no longer reports failures in the released package**: `make test` ran two checks that replay committed VASP and GULP output from `example01`, which the release deliberately leaves out to keep the package small. With nothing to replay they failed rather than skipping, so six tests failed for anyone verifying a fresh install of 0.90. They now skip, with the reason given.

- **`bin/sod_check_docs.sh`**: new consistency check over the documentation — the version string across `README.md`, `docs/conf.py`, `VERSIONS.md` and the program banners, that every example is described, that every wrapper script is documented, that the README's internal links resolve, and that the docs build without warnings.

- **`sod_mace` fails early on an incompatible ALCHEMI toolkit**: the backend uses a few private `nvalchemi-toolkit` 0.2.0 APIs, and now verifies them before a run starts, naming the missing attribute and the validated version, instead of raising an `AttributeError` minutes into a relaxation.

- **`genersod` prints its version banner** like the other eleven executables.

## Version 0.90 (28 August 2026)

Feature release. Machine-learning interatomic potentials, species-resolved SQS,
and the PME → CPME rename.

- **MACE machine-learning potentials (`sod_mace.sh`, new `pysod/` package)**: SOD can now evaluate energies with a MACE foundation model instead of a DFT or forcefield code, over the same `nXX/cYY/` layout as every other calculator. `sod_mace.sh` writes `nXX/ENERGIES` in the usual two-column `m  E` format, so `sod_stat.sh` consumes it unchanged. Structures are evaluated in GPU batches through NVIDIA ALCHEMI, with optional cuEquivariance kernels; on the documented 710-configuration benchmark, batching single-point evaluation up to 64 structures at a time is a 16x speedup over one structure at a time, and batched relaxation reaches 9.15x, with `-batchmode refill` worth a further 1.14x-1.34x at equal batch size. `-relax` performs FIRE2 relaxation of the atomic positions, `-relaxcell` also relaxes the cell against the MACE stress at a chosen `-pressure`, and `-lattice` writes `nXX/CELL` in the same columns as `sod_vasp_cell.sh`. Options can be set in a `mace_settings.yaml` file (in `nXX/`, else in SODPROJECT) as well as on the command line. The whole Python stack is optional: nothing in the Fortran build or the shell wrappers requires it. See `docs/mace.rst` and `pysod/README.md`.

- **Species-resolved pair SQS for multinary and multi-target disorder (`sqssod`)**: SQS scoring now handles several substituting species and several substituted sublattices at once, scoring each species pair channel separately rather than collapsing the problem to a binary one. `INSQS` gained a user-configurable `n_top_sqs`; `sqssod` reports a top-10 ranking with a deterministic tie-break instead of a full ordering. New example19 (equimolar fcc CoCrFeNi, a Cantor-alloy subsystem) demonstrates the random-sampling → SQS route for a composition far beyond enumeration.

- **High-entropy alloy compositions (INSOD)**: the per-target cap rose from 3 to 5 substituents (6 species including the parent), which is what quaternary and quinary high-entropy alloys need.

- **PME renamed to CPME (Constrained Periodic Motif Expansion)**: the Hamiltonian, its executable, its wrapper and its input file were renamed throughout — `pmesod` → `cpmesod`, `sod_pme.sh` → `sod_cpme.sh`, and the model file is now `cpme.model` (with `cpme.model.tmp` for the auto-generated suggestion). Output directories are `CPME0/`, `CPME1/` and `CPMEh/`. The name describes what the expansion does: it is constrained to motifs that are periodic in the supercell.

- **Internal energies and enthalpies are separate outputs**: `ENERGIES` always holds the internal energy E. A calculation run under an applied pressure additionally writes `ENTHALPIES` (H = E + PV), through the new `sod_vasp_enth.sh` and `sod_gulp_enth.sh`, and through `sod_mace.sh` at nonzero `-pressure`. Previously a constant-pressure run silently put an enthalpy in `ENERGIES`.

- **`CELL` requires complete ensemble coverage**: `CELL` rows are positional, with no configuration index, so a gap would shift every later row onto the wrong configuration. Every producer (`sod_vasp_cell.sh`, `sod_gulp_cell.sh`, `sod_mace.sh -lattice`) now writes `CELL` only when results cover the level's whole `ENSEMBLE`, and reports why when they do not. Sparse calibration subsets remain valid indexed `ENERGIES` input.

- **ENSEMBLE format version 2 is no longer read** (breaking): all programs require the version 3 `ENSEMBLE` written since 0.80. Regenerate an old file with `sod_comb.sh` or `sod_random.sh`.

- **Energy and free-energy collectors no longer trip over spare directories**: a directory sitting beside the configuration directories whose name begins with `c` but is not a configuration index — `c01_backup`, `cpme_test_ref` and the like — made `sod_gulp_ener.sh` and the other collectors print a shell arithmetic error and exit non-zero on an otherwise successful run. All of them now skip such directories, as the `CELL` extractors already did.

- **`sod_gulp_free.sh` fixed** (data loss): run from SODPROJECT, it deleted every `nXX/ENERGIES` and appended all levels' free energies into a single file in the wrong directory — while exiting 0, so the loss was silent. It also wrote a bare one-column file that `statsod` cannot read at all. It now behaves like the other collectors: one indexed two-column `nXX/ENERGIES` per level.

## Version 0.84 (8 July 2026)

Minor release.

- **Multinary and multisite random sampling (`randomsod` / `sod_random.sh`)**: the uniform random sampler now supports the full substitution model — multiple target sublattices (multisite) and multiple substituting species per site (multinary) — exactly as `combsod`. Previously it handled only single-target binary substitution and silently ignored any extra targets or species. Each draw places a uniform colored assignment on every sublattice, and the output `ENSEMBLE` is byte-compatible with the `combsod` format (one column group per target/species). With `-sym on`, draws fold to colored symmetry orbits using the *same* canonical representative as `combsod`, so the two `ENSEMBLE` files match line-for-line on the sampled orbits — differing only in the degeneracy column (`combsod` writes exact orbit sizes, `randomsod` writes visit counts).
- **`sod_pme.sh` / `sod_mc.sh`: `-model` argument validation**: both wrappers now fail fast with `Error: -model requires a filename argument.` when `-model` is given with no filename, instead of looping forever (`sod_mc.sh`) or forwarding a silent error (`sod_pme.sh`).
- **`sod_comb.sh` / `sod_gener.sh`: robust FILER parsing**: the wrappers now take the FILER value from the last non-blank, non-comment line of `INSOD`, so a trailing blank or comment line no longer triggers a spurious "unknown FILER" error.
- **`combsod`: fixed EQMATRIX corruption for large supercells**: the supercell symmetry-operator arrays are now allocated to the exact operator count, fixing out-of-bounds corruption when `nop = nop1 × na × nb × nc` exceeded the old fixed bound of 10000 (e.g. Fm-3m 5×5×5 = 24000). The corruption previously made `sqssod`/`gqssod` score every configuration as zero.

## Version 0.83 (22 June 2026)

Minor release.

- **Warren-Cowley SRO parameters (`gqssod`)**: `gqssod` now writes `wc_parameters.dat` with thermally averaged Warren-Cowley α_n parameters for each symmetrically distinct pair-cluster shell at every requested temperature, plus the T→∞ (random) limit.
- **`AUTHORS.md`**: contributors file converted to Markdown with GitHub profile links.
- **`sod_random.sh` flags shortened** (breaking): `-nconfigs` → `-nconf`, `-symmetry` → `-sym`. Default symmetry changed to `off`. The shell wrapper now catches a missing `-nconf` and prints its own usage. Error messages point to `sod_random.sh` rather than the bare `randomsod` executable.
- **example18 documentation corrected**: all four description locations now accurately describe the actual committed example (4×4×4 supercell, nsubs=96, random-sampling workflow with 50,000 draws).

## Version 0.82 (20 June 2026)

- **Standalone random sampler (`randomsod` / `sod_random.sh`)**: uniform random sampling is now a separate executable, split out of `mcsod`. It draws configurations with no energy evaluation at all and writes `nXX/random/ENSEMBLE`, ready for the usual structure-writer → DFT → `statsod` path — the sampling counterpart of `combsod` for levels too big to enumerate. Run as `sod_random.sh -nconfigs <N> [-symmetry on|off] [-seed clock|<int>]`. As a result, `mcsod` and `INMC` are now Metropolis-only: the `INMC` sampler line is gone (first value is now the symmetry flag), and uniform workflows move to `sod_random.sh`. The shared geometry/symmetry/sampling helpers now live in a new `config_sampling` module used by both programs.

- **Clearer INSOD error for a wrong `nsp`**: if the `nsp` count does not match the number of species actually listed on the `symbol` (or `natsp0`) line, SOD now stops with an explicit message instead of silently dropping the extra species and running with the wrong composition.

- **Faster Metropolis MC (`mcsod`)**: each step now updates the PME energy *incrementally* for the single occupied↔hole swap instead of recomputing the whole cluster expansion (per-step cost O(L^k)/O(H^k) → O(L^(k-1))/O(H^(k-1)), with L substitutions, H holes, k the expansion order), and draws the swapped sites directly from a maintained hole list (O(1)) rather than by rejection sampling. Results are numerically equivalent to a full recompute at every step; the speedup grows with cell size, expansion order, and dilute or near-full filling, and the high-filling move-proposal failure mode is removed. A `test_pme_delta` driver checks the incremental update against the full evaluator, and the MC regression tests compare energetics within tolerance (the optimised walk explores a different but statistically equivalent fixed-seed trajectory).

- **`AUTHORS` file**: lists all contributors to the SOD package.

## Version 0.81 (14 June 2026)

Minor release.

- **Parent-structure molecules (`@NAME`)**: a molecule can now be part of the parent structure, not only a substitution — mark a parent site `@NAME` in the INSOD `symbol` list and `genersod` materialises `NAME.xyz` (independent random orientation per site) in the calculation inputs. New `example18` (MAPbI3–MAPbBr3 solid solution with a methylammonium parent).
- **Geometry-only uniform MC sampling**: `mcsod` sampler 2 now runs with no PME Hamiltonian or reference energies, writing just the configuration ensemble — the practical route to build large SQS/GQS candidate sets before any DFT.
- **MC output layout reorganised** (sampling method first, Hamiltonian variant second): Metropolis output is now `nXX/MCT_TTTK/PMEx/` and uniform output `nXX/MCU/` (per-variant energies under `nXX/MCU/PMEx/`); `sod_mcstat.sh` is now run from `nXX/`.
- **Minor**: `combsod` accepts `nsubs = 0` (parent) and `nsubs = Nsites` (fully substituted); a missing `INCAR` is now a warning rather than a fatal error for VASP output.

## Version 0.80 (5 June 2026)

### New features 

- **Periodic Motif Expansion (PME)**: new executable `pmesod` and wrapper `sod_pme.sh`.
  Fits an effective Hamiltonian `E(c) = V0 + ε0 + ε1·T1(c) + … + εK·TK(c)` (up to 4-body
  motif interactions) from ab initio reference energies at low substitution levels (n00–n0p,
  low side) and/or high substitution levels (n(M-p)–nM, high side), then predicts energies
  for all configurations at a chosen target level. `ε0` is an additive energy offset (eV,
  default 0); `ε1`–`εK` are multiplicative scale factors per motif order (default 1). Three
  variants: PME0 (low-side only), PME1 (high-side only), and PMEh (hybrid, blending both
  sides with a piecewise power-law scheme controlled by `alpha` (sharpness, default 2.0) and
  `eta` (log-scale asymmetry, default 0.0)). Building the Hamiltonian follows a two-phase workflow:
  Phase 1 (referencing) runs DFT at boundary concentrations (low/high substitution),
  calls `sod_pme.sh` to fit the V terms and produce `pme.model.tmp` with a
  bisection-selected calibration list; Phase 2 (calibrating, optional) runs DFT for a
  small set of configurations at the target concentration, then re-runs `sod_pme.sh` with
  `pme.model` to fit the ε correction terms. Predicted energies are written to `nXX/PME0/ENERGIES`, `nXX/PME1/ENERGIES`,
  and `nXX/PMEh/ENERGIES`; fitted ε values and RMSE (before → after calibration) are
  reported on stdout.

- **Monte Carlo (MC) sampling**: new executable `mcsod` and wrapper `sod_mc.sh`.
  Samples configuration space at finite temperature using a PME Hamiltonian. Supports four
  modes: Metropolis reduced (symmetry-deduplicated) and Metropolis full (explicit trajectory),
  Uniform random reduced and Uniform random full. `sod_mc.sh` reads a `TEMPERATURES` file and
  runs one single-temperature MC job per entry. Output follows a sampling-method-first layout:
  Metropolis writes `nXX/MCT_TTTK/PMEx/`, while uniform sampling writes its (Hamiltonian-
  independent) geometry sample to `nXX/MCU/` and, when reference energies are available, the
  per-variant energies to `nXX/MCU/PMEx/`. Uniform sampling runs geometry-only without any PME
  training data.
  `INMC` specifies the sampler, equilibration/production step counts, restart probability,
  and starting configuration. Block-average standard error of the mean is reported for
  reliability assessment. MC output (ENSEMBLE + ENERGIES) is compatible with `sod_stat.sh`
  for thermodynamic post-processing.

- **Thermodynamic integration**: new executable `mcstatsod` and wrapper `sod_mcstat.sh`.
  Integrates the MC internal energy U(T) over inverse temperature β to obtain the Helmholtz
  free energy F and entropy S at each sampled temperature, using S(T→∞) = kB ln C(npos, lev)
  as the high-temperature anchor. Output format matches `statsod` (T, E, F, S columns per
  supercell). Run from `nXX/` after completing a multi-temperature MC run (it reads
  `../TEMPERATURES` and each `MCT_*K/PMEx/` directory, with `PMEx` taken from `../pme.model`,
  defaulting to `PMEh`).

- **Special Quasirandom Structures (SQS)**: new executable `sqssod` and wrapper `sod_sqs.sh`.
  Identifies the configurations at a given substitution level whose cluster 
  correlations best match ideal random mixing. Controlled via `INSQS`; results ranked in
  `OUTSQS`.

- **Generalized Quasirandom Structures (GQS)**: new executable `gqssod` and wrapper
  `sod_gqs.sh`. Extends SQS to finite temperature by computing Boltzmann-weighted thermal
  averages of pair correlations from `ENERGIES` and `TEMPERATURES`; results written to
  `OUTGQS`. Run `sod_gqs.sh` from SODPROJECT; requires `ENERGIES`.

### Removed

- **SPBE (simple pair-based extrapolation)** retired. The PME Hamiltonian supersedes it with
  higher-order many-body terms, hybrid low/high fitting, and recalibration support.
  `spbesod.f90` and `sod_spbe.sh` are no longer distributed.

### Bug fixes

- **`mcsod` — MC sampling counting (critical)**: rejected Metropolis moves were not counted
  as samples of the current state, giving incorrect energy averages and ENSEMBLE weights. Now
  all steps (accepted and rejected) accumulate correctly. The fix affects E_ave and all
  derived thermodynamic quantities from Metropolis runs.

### Other improvements

- **OUTSOD renamed ENSEMBLE**: the configuration-list file is renamed from `OUTSOD` to
  `ENSEMBLE` in all programs, shell scripts, Makefile, and examples. A new self-documenting
  v3 format is introduced: line 1 names the ensemble type and records the configuration
  count and sum of degeneracies; one target line per substitution target follows. The old v2
  format (with or without `#` header) remains fully readable for backward compatibility.
- **Continuous integration**: GitHub Actions workflow runs `make all && make test` on Ubuntu
  and macOS on every push and pull request.
- **Regression test suite** expanded to 28 tests covering combsod (examples 02–14), genersod
  FILER variants (example01), statsod/gcstatsod (example05), pmesod/mcsod (example15),
  sqssod/gqssod (example16), and pmesod PMEh (example17). The thermodynamics tools are now
  covered end to end: thermodynamic integration (`mcstatsod`) over a multi-temperature MC
  ladder (example15), exact canonical statistics (`statsod`) over the full
  8043-configuration enumeration (example16), and a physics cross-check confirming that
  Monte Carlo plus thermodynamic integration reproduces the exact free energy of the full
  enumeration to within ~2 meV on the same PME Hamiltonian (example15).
- **Robust INSOD reading**: all six programs that read INSOD now tolerate blank lines and
  `#` comment lines anywhere in the file.
- `example15` added: Si/Ge substitution in α-quartz 2×2×2 supercell (24 Si sites);
  demonstrates PME fitting with low-side (n00–n03) and high-side (n21–n24) training data,
  and Metropolis MC sampling at 300 K.
- `example16` added: SQS/GQS workflow for 8 Ni/Mg substitutions in 2×2×2 MgO rocksalt
  (8043 configurations); includes OUTSQS and OUTGQS reference outputs.
- `example17` added: Al/Fe substitutions in LaFeO3 3×3×3 perovskite supercell (27 Fe sites);
  demonstrates third-order PMEh for target level n04 with full DFT reference energies.

## Version 0.71 (April 2026)

### New features

- **Multi-target multi-nary substitution**: multi-species and multi-nary substitution combined — each target site independently supports one or more new species. INSOD `nsubs` uses one line per target, each line with space-separated counts (e.g. site 1 `1 1` for ternary, site 2 `1` for binary). Joint configurations enumerated under full crystal symmetry.

### Other improvements

- `example13` added: 3-target multi-species example La₁₋ₓSrₓFe₁₋ᵧMnᵧO₃₋ᵤ in a 2×2×2 supercell (Pm-3m) — simultaneous binary substitution on La, Fe, and O sites (1 substitution per site), including an O vacancy (`%O` syntax); 6 inequivalent configurations from 1536 total. FILER=-1.
- `example14` added: first multi-target multi-nary example — La₁₋ₓ₋ᵧSrₓBaᵧMnᵤFe₁₋ᵤO₃ in a 2×2×2 supercell (Pm-3m); ternary substitution on the La site (1 Sr + 1 Ba) combined with binary substitution on the Fe site (1 Mn); 3 inequivalent configurations from 448 total. FILER=-1.
- `run_tests.sh` added: formal regression test suite (20 tests); runs in isolated temp directories and diffs against committed references: all n*/OUTSOD for examples 02–14 (combsod), one structure file per calculator for example01/FILER* (genersod), thermodynamics/DATA/SPECTRA averages for example05/n02 (statsod canonical), and thermodynamics/DATA/SPECTRA averages for example05/test_gcstat over n10–n16 (gcstatsod grand-canonical).
- `example05` corrected: substitution direction changed from Sn-in-Zr (La₂Zr₂O₇ parent) to Zr-in-Sn (La₂Sn₂O₇ parent) to match the CASTEP calculations. All n00–n16 OUTSODs regenerated; ENERGIES/DATA/SPECTRA from DFT (archive) renumbered accordingly. Statistical analysis folders renamed to `x250` (x=0.25 Zr) and `x750` (x=0.75 Zr).

## Version 0.70 (April 2026)

### New features

- **Multi-species substitution**: simultaneous binary substitution on up to five distinct crystallographic sites in a single run. Joint configurations enumerated under the full crystal symmetry. INSOD uses one `nsubs` line per target site (e.g. `2` on line 1, `1` on line 2).
- **Multi-nary substitution**: simultaneous placement of up to three new species on a single target site (binary, ternary, or quaternary disorder). INSOD uses space-separated `nsubs` values (e.g. `1 2` for ternary, `2 2 2` for quaternary). Combinatorial rank: C(npos,n_tot) × C(n_tot,n1) × C(n2+n3,n2). Directory named `nXX`, `nXX_YY`, or `nXX_YY_ZZ`.
- **Molecules (`@NAME`) and vacancies (`%NAME`)**: `@NAME` in `newsymbol` reads `NAME.xyz`, places the molecule with its centre of mass at the site with a uniform random orientation (Shoemake quaternion method), and expands into individual atoms in the output. `%NAME` creates a vacancy (atom omitted from all output). Supported for all calculators except LAMMPS.
- **Entropy statistics** printed at each substitution level: maximum ensemble entropy kB ln(Ω) and ideal mixing entropy per cell (both in meV/K), plus the percentage of ideal entropy captured by the supercell size.
- **`kim_model` directive** in `genersod`: supports OpenKIM potential models in GULP input files.
- **Energy extraction scripts**: `sod_qe_ener.sh` and `sod_castep_ener.sh` added; `sod_vasp_ener.sh` updated for `nXX/cYY/` directory structure.

### Bug fixes

- **`sod_spbe1.sh`**: `head -1 ../OUTSOD` retrieved the `# SOD OUTSOD format version 2` comment line, not the substitution line, so `nsites` was set to the string `"format"` and the subsequent arithmetic failed. Fixed by using `grep -v "^#" ../OUTSOD | head -1` to skip comment lines.
- **`combsod.f90` recursive mode**: positions read from the previous-level OUTSOD (global atom indices) were not converted to local indices before use in the candidate-expansion loop, causing spurious out-of-range positions when the target species is not species 1 (atini > 1). Fixed by subtracting `atini − 1` from `indconf_prev` after reading.
- **`spbesod.f90` index offset**: OUTSOD position indices are global; `eqmatrix` uses local (1..npos) indices. Direct use caused out-of-bounds access for any system where the target species is not species 1. Fixed by computing the minimum OUTSOD position and subtracting it as an offset before all `eqmatrix` and `de1`/`de2` lookups.
- `random_seed()` moved to program start in `gcstatsod.f90`.

### Limits expanded

- `natmax` increased to 10,000 across all modules.
- `nsubsmax` increased to 1,000.
- `nconfmax` increased to 1,000,000 and `ntempmax` to 1,000 in `statsod` and `gcstatsod`; large arrays (`p`, `spec`) made allocatable to keep memory on-demand.

### Code quality (Fortran 2003/2008)

- `use iso_fortran_env, only: int64, real64` added to all programs and functions; `real(kind=8)`, `real*8`, `integer(kind=8)` replaced with portable named constants throughout.
- `kb` and other physics/tolerance parameters given correct double-precision literals (`_real64` suffix) where previously assigned single-precision values.
- `implicit none` added to `ccf.f90` (was missing); `character*N` old-style syntax replaced with `character(len=N)` throughout all source files.
- All bare `stop` statements replaced with `stop 1`, so shell scripts can detect failure via `$?`.
- Symmetry-image gather loops converted to `do concurrent` in all four enumeration branches of `combsod.f90`.
- All geometric and coordinate variables in `cell.f90`, `genersod.f90`, and `combsod.f90` (operator matrices, coordinate arrays, cell parameters, molecule geometry) upgraded from default `real` to `real(real64)` for full double-precision consistency.

### Other improvements

- **INSOD format change**: multi-species `nsubs` no longer uses `/` separator; now one line per target site (e.g. `2` on line 1, `1` on line 2). Old `/`-format INSODs will produce a clear error message.
- `invertOUTSOD.f90` output format generalised to support arbitrary number of sites.
- Guard added to `gcstatsod` to abort on multi-species/multi-nary OUTSOD (not supported for grand-canonical).
- `example07` updated to P4(3)32 γ-Fe₂O₃ structure (a=8.344 Å, 29 inequivalent configurations).
- `example11` added: equimolar NiCoFeCr Cantor subsystem alloy (2×2×2 FCC primitive cell, 8 atoms), 23 inequivalent configurations from 2520 total.
- `example12` added: complex perovskite La₀.₇₅Sr₀.₂₅Mn₀.₂₅Fe₀.₇₅O₃ in a 2×2×2 supercell — multi-species example (binary on two sites simultaneously); 13 inequivalent configurations from 784 total; GULP input files (FILER=1).

## Version 0.62 (March 2026)

### Removed
- MAPPER feature removed from `combsod.f90`, `genersod.f90`, and all INSOD example files.
- METADISE output (FILER=3) removed; LAMMPS is now FILER=2.
- Classical-code block in INSOD (`ishell`, `newshell`, `genxtl`, `genarc`) removed from `combsod.f90` and `genersod.f90`. INSOD files no longer require this block.
- `sod_bcs2sgo.sh` and `bcs2sgo.awk` deleted (obsolete).

### New features
- Template-based input generation for all calculators: `template_input.gin` (GULP), `template_in.lammps` (LAMMPS), `template_castep.cell` (CASTEP), `template_pw.in` (QE). VASP generates POSCAR directly.
- Uniform `nXX/cYY/` directory structure for all calculators including CIF.
- `# sod_type_map` lines in template files are parsed by SOD for type mapping and stripped from all generated input files.
- Quantum ESPRESSO support (FILER=13): `genersod` inserts `CELL_PARAMETERS {angstrom}` and `ATOMIC_POSITIONS {crystal}` blocks.
- Recursive enumeration: `nsubs_min` to `nsubs_max` range in INSOD; each level uses OUTSOD from the previous level.
- Output filenames (e.g. `c01.gin`) now use minimum zero-padding for correct sorting, removing the previous 99999-configuration limit.
- Five parallel `example1_*` examples (GULP, LAMMPS, VASP, CASTEP, QE): Ni/Mg in MgO rocksalt, nsubs=4.

### Bug fixes
- Stress-volume correction in `gcstatsod.f90` now uses the correct second-order Birch-Murnaghan expression.
- Monoclinic cell parameter `a` calculation in `sod_gulp_cell.sh` corrected (was dividing by uninitialised awk variables).
- `sod_vasp_cell.sh` rewritten to support all lattice types, matching `sod_gulp_cell.sh`.
- `sod_comb.sh` no longer performs `cd ..` when FILER=0.
- `sod_gener.sh` `cd ..` outside FILER≠0 block corrected; dynamic zero-padding applied.
- `Einf` uninitialised before accumulation in `statsod.f90` corrected.
- `nat=sum(natsp(1:nsp))` fix in `genersod.f90` (previously summed uninitialised elements).
- `sod_gcstat.sh` now copies `DATA` files alongside OUTSOD/ENERGIES.
- Guard added to avoid `log(0)` when x=0 or x=1 in entropy calculations.
- `example1_lammps/INSOD` had `sptarget=2` instead of `sptarget=1`.
- `example5/INSOD` corrected (`natsp0` was `1 1 1` instead of `1 1 2`).

### Other improvements
- `combsod` now only performs combinatorics; `sod_comb.sh` calls `genersod` automatically when FILER≠0.
- GULP library reduction: sections with no relevant interactions are omitted from the copied library.
- Makefile now cross-platform (macOS and Linux).
- Hardcoded pi values replaced with `acos(-1.0)` throughout.

## Version 0.52 (October 2023)
- Made arrays allocatable in gcstatsod.f90, which fixes some issue with gfortran compilar in MacOS.
- Added the sgo/ folder, which was missing in previous release by mistake. 

## Version 0.51 (September 2023)
- Some memory issues in combinatorial generation sorted.
- sod_gener.sh generates VASP or GULP input files from OUTSOD (and INSOD) without having to redo the combinatorics.
- Grand-canonical analysis now posible using sod_gcstat.sh.
- Stress-volume correction  in the grand-canonical analysis implemented.
- Canonical analysis with sod_stat.sh now includes the full disorder limit.
- Other minor improvements in format for statistical analysis. 
- Statistics on spectra (both canonical and grand-canonical)
- Input files for CASTEP can be created. 

## Version 0.47 (January 2019)
- spbe module creates INSPBE template (INSPBE.tmp) to make rescaling step easier. 
- Bug that led to error message when running combinatorics with FILER=0 corrected. 

## Version 0.46 (November 2018)

- Increased precision of reals in spbe to avoid numerical errors. 	
- The spbe calculation can be corrected using two reference energies, via INSPBE file.   
- Example 1 was changed to illustrate the spbe procedure (including rescaling).

## Version 0.44 (October 2018)

- Added a simple pair-based extrapolation (spbe) calculator to predict energies for n>2 substitutions from the energies for n=0, 1, and 2 substitutions. 
- Makefile added; just type ```make all``` from ROOTSOD to compile all the fortran code.
- All scripts converted to bash; they are all in the ROOTSOD/bin folder, and have extension .sh
- Resolved issue in scripts related to number of columns from ```ls -al``` command.
- ```sod_stat``` now calculates thermodynamic properties at T=1K, 300K, 1000K and in the limit of a very high temperature, if the TEMPERATURES file does not exist.
