## Version 0.44 (October 2018)

- Added a simple pair-based extrapolation (spbe) calculator to predict energies for n>2 substitutions from the energies for n=0, 1, and 2 substitutions. 
- Makefile added; just type ```make all``` from ROOTSOD to compile all the fortran code.
- All scripts converted to bash; they are all in the ROOTSOD/bin folder, and have extension .sh
- Resolved issue in scripts related to number of columns from ```ls -al``` command.
- ```sod_stat``` now calculates thermodynamic properties at T=1K, 300K, 1000K and in the limit of a very high temperature, if the TEMPERATURES file does not exist.

## Version 0.46 (November 2018)

- Increased precision of reals in spbe to avoid numerical errors. 	
- The spbe calculation can be corrected using two reference energies, via INSPBE file.   
- Example 1 was changed to illustrate the spbe procedure (including rescaling).

## Version 0.47 (January 2019)
- spbe module creates INSPBE template (INSPBE.tmp) to make rescaling step easier. 
- Bug that led to error message when running combinatorics with FILER=0 corrected. 

## Version 0.51 (September 2023)
- Some memory issues in combinatorial generation sorted.
- sod_gener.sh generates VASP or GULP input files from OUTSOD (and INSOD) without having to redo the combinatorics.
- Grand-canonical analysis now posible using sod_gcstat.sh.
- Stress-volume correction  in the grand-canonical analysis implemented.
- Canonical analysis with sod_stat.sh now includes the full disorder limit.
- Other minor improvements in format for statistical analysis. 
- Statistics on spectra (both canonical and grand-canonical)
- Input files for CASTEP can be created. 

## Version 0.52 (October 2023)
- Made arrays allocatable in gcstatsod.f90, which fixes some issue with gfortran compilar in MacOS.
- Added the sgo/ folder, which was missing in previous release by mistake. 

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

## Version 0.71 (April 2026)

### New features

- **Multi-target multi-nary substitution**: multi-species and multi-nary substitution combined — each target site independently supports one or more new species. INSOD `nsubs` uses one line per target, each line with space-separated counts (e.g. site 1 `1 1` for ternary, site 2 `1` for binary). Joint configurations enumerated under full crystal symmetry.

### Other improvements

- `example13` added: 3-target multi-species example La₁₋ₓSrₓFe₁₋ᵧMnᵧO₃₋ᵤ in a 2×2×2 supercell (Pm-3m) — simultaneous binary substitution on La, Fe, and O sites (1 substitution per site), including an O vacancy (`%O` syntax); 6 inequivalent configurations from 1536 total. FILER=-1.
- `example14` added: first multi-target multi-nary example — La₁₋ₓ₋ᵧSrₓBaᵧMnᵤFe₁₋ᵤO₃ in a 2×2×2 supercell (Pm-3m); ternary substitution on the La site (1 Sr + 1 Ba) combined with binary substitution on the Fe site (1 Mn); 3 inequivalent configurations from 448 total. FILER=-1.
- `run_tests.sh` added: formal regression test suite (20 tests); runs in isolated temp directories and diffs against committed references: all n*/OUTSOD for examples 02–14 (combsod), one structure file per calculator for example01/FILER* (genersod), thermodynamics/DATA/SPECTRA averages for example05/n02 (statsod canonical), and thermodynamics/DATA/SPECTRA averages for example05/test_gcstat over n10–n16 (gcstatsod grand-canonical).
- `example05` corrected: substitution direction changed from Sn-in-Zr (La₂Zr₂O₇ parent) to Zr-in-Sn (La₂Sn₂O₇ parent) to match the CASTEP calculations. All n00–n16 OUTSODs regenerated; ENERGIES/DATA/SPECTRA from DFT (archive) renumbered accordingly. Statistical analysis folders renamed to `x250` (x=0.25 Zr) and `x750` (x=0.75 Zr).
