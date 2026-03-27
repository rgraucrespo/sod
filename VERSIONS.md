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

## Version 0.60 (March 2026)
- Removed MAPPER feature entirely from `combsod.f90`, `genersod.f90`, and all INSOD example files.
- Bug fix: stress-volume correction in `gcstatsod.f90` now uses the correct second-order Birch-Murnaghan expression.
- Memory guards added to all `allocate` calls in `combsod.f90`; prints an informative error and exits if memory is insufficient.
- Edge cases handled gracefully: nsubs=0 and nsubs=npos now produce a single valid configuration and exit cleanly; nsubs>npos raises an error.
- Guard added to entropy calculation to avoid `log(0)` when x=0 or x=1.
- `sod_vasp_cell.sh` rewritten to support all lattice types (cubic, tetragonal, orthorhombic, hexagonal, rhombohedral, monoclinic, triclinic), matching the functionality of `sod_gulp_cell.sh`.
- Bug fix: monoclinic cell parameter `a` calculation in `sod_gulp_cell.sh` corrected (was dividing by uninitialised awk variables).
- Bug fix: `sod_comb.sh` no longer performs `cd ..` when FILER=0 (no CALCS directory is created in that case).
- Optional xtl and arc output directives for GULP input files, controlled by two new `genxtl` and `genarc` flags in INSOD.
- Output filenames (e.g. `c01.gin`) now use the minimum zero-padding needed for correct sorting, rather than a fixed 5-digit width. This also removes the previous limit of 99999 configurations.
- Deleted `sod_bcs2sgo.sh` and `bcs2sgo.awk` (obsolete).
- Bug fix: example5/INSOD corrected (`natsp0` was `1 1 1` instead of `1 1 2` for the two-oxygen-site pyrochlore structure).
- Bug fix: `Einf` uninitialised before accumulation in `statsod.f90` (caused wrong infinite-temperature energy).
- Bug fix: `nat=sum(natsp)` in `genersod.f90` corrected to `nat=sum(natsp(1:nsp))` to avoid summing uninitialised array elements.
- Bug fix: `sod_gcstat.sh` now copies `DATA` files alongside OUTSOD/ENERGIES for grand-canonical data averaging.
- Bug fix: `sod_gener.sh` `cd ..` outside the FILER≠0 block corrected; dynamic zero-padding applied (consistent with `sod_comb.sh`).
- Refactoring: `combsod` now only performs combinatorics; `sod_comb.sh` calls `genersod` automatically when FILER≠0, removing duplicated file-generation code.
- Code cleanup: removed dead code and commented-out blocks; hardcoded pi values replaced with `acos(-1.0)` in `cell.f90` and `peaks2spec.f90`; obsolete `bcs2sgo` files deleted.
- Makefile now cross-platform (macOS and Linux).

## Version 0.61 (March 2026)
- Graceful failure when configuration count overflows integer range; 
- Blank line after FILER value no longer required in INSOD when no classical-code block follows (filer >= 10).
- Calculation input folder renamed from `CALCS` to `nXX` (XX = zero-padded number of substitutions), with existence check to prevent overwriting.
- Classical-code block in INSOD (ishell, newshell, genxtl, genarc) is now read for all filer < 10, including filer = 0.
- Bug fix: stale `indcount > 99999` guards removed from CASTEP case in `genersod.f90`.
- New feature: Quantum ESPRESSO input files (filer = 13, extension `.pwi`). `genersod` writes `CELL_PARAMETERS {angstrom}` and `ATOMIC_POSITIONS {crystal}`, sandwiched between user-supplied `top.qe` and `bottom.qe` templates.

