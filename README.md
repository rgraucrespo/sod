# SOD 0.71 - Notes for users

SOD (standing for Site-Occupancy Disorder) is a package of tools for the computer modelling of periodic systems with site disorder, using the supercell ensemble method. 

The package is distributed under the [GPL licence](https://github.com/ypriverol/sod/blob/master/LICENSE.md). 

You can find below the essential info needed to use SOD. Please note that SOD authors can give only limited support to users.


## Functionalities

- Identification of all inequivalent configurations of site substitutions in an arbitrary supercell of an  initial target structure with any space group symmetry.
- Calculation of the degeneracies of configurations.
- Generation of input files for GULP, LAMMPS, VASP, CASTEP and Quantum ESPRESSO.
- Simple extrapolation of energies from low to high concentrations within a supercell.
- Statistical mechanics processing of output using either canonical or grand-canonical ensembles.


## Content of the folders

- sod(version)/src contains the source files.
- sod(version)/sgo is a library of space group operators (e.g. 131.sgo contains the operators of the space group 131).
- sod(version)/bin contains the executables. Compile for your platform using `make all`.
- sod(version)/examples contains fourteen examples, covering a range of structure types and substitution modes: binary, multi-nary, multi-species, multi-target multi-nary, vacancies and molecules.

## Examples

The **example01/\*** series all use the same physical system — Ni/Mg substitutions in a 2×2×2 supercell of rocksalt MgO (space group Fm-3m) — and illustrate input generation for five different calculators using the same INSOD file (differing only in FILER):

- **example01/FILER1_gulp**: Ni/Mg substitutions in a 2×2×2 MgO rocksalt supercell (32 Mg sites), nsubs=4. GULP (FILER=1) with `catlow.lib` (Buckingham potentials, core-shell model for Ni and O). Demonstrates the `library` directive and `sod_type_map` (SOD species `O` maps to GULP type `O2`).

- **example01/FILER2_lammps**: Same system, nsubs=4. LAMMPS (FILER=2). Demonstrates `template_in.lammps` with `sod_type_map` lines and core-shell representation in the generated `conf.data`.

- **example01/FILER11_vasp**: Same system, nsubs=4. VASP (FILER=11). SOD generates `POSCAR` into each `cYY/` folder; the user must supply `INCAR`, `KPOINTS` and `POTCAR`.

- **example01/FILER12_castep**: Same system, nsubs=4. CASTEP (FILER=12). Demonstrates `template_castep.cell`; the user must supply the `.param` file.

- **example01/FILER13_qe**: Same system, nsubs=4. Quantum ESPRESSO (FILER=13). Demonstrates `template_pw.in`; replace `Ni.upf`, `Mg.upf`, `O.upf` with your pseudopotential files before running.

- **example02**: Al/Fe substitution in magnetite Fe₃O₄ (1×1×1 cubic spinel unit cell, space group Fd-3m). Al substitutions from 0 to 16 across 24 mixed Fe sites (8 tetrahedral + 16 octahedral) are enumerated (FILER=-1); the `n04/` folder contains the OUTSOD for the 4-substitution case (99 inequivalent configurations). The two structurally distinct Fe environments make this a non-trivial site-disorder problem.

- **example03**: Fe/Sb disorder in a 2×2×2 supercell of rutile FeSbO₄ (space group P4₂/mnm). 8 M(IV) sites replaced by equal numbers of Fe(III) and Sb(V) give 180 inequivalent configurations out of 12870 total (FILER=-1; no calculation input files generated). The OUTSOD is pre-computed in `n08/`. See [Grau-Crespo et al., Chemistry of Materials (2004)](https://pubs.acs.org/doi/abs/10.1021/cm035271y).

- **example04**: Al/Fe substitution in a 2×2×2 supercell of LaFeO₃ perovskite (cubic approximation, space group Pm-3m). 4 Al in 8 Fe sites give 6 inequivalent configurations. GULP input files (FILER=1) using the Bush et al. Buckingham potentials with core-shell model, defined inline in `template_input.gin`. Also includes SPBE energy extrapolation (subfolders `n04/spbe0` and `n04/spbe1`) and a grand-canonical statistics example (`x250`, composition x=0.25).

- **example05**: Zr/Sn substitution in the La₂Sn₂O₇–La₂Zr₂O₇ solid solution (pyrochlore, space group Fd-3m). INSOD enumerates all compositions (nsubs=0:16, FILER=-1), spanning pure La₂Sn₂O₇ (n00) to pure La₂Zr₂O₇ (n16). Pre-generated CASTEP structure files for nsubs=2 are in `n02/c1-c3/` as a usage illustration. Pre-computed DFT results (ENERGIES, DATA, SPECTRA) across all n00–n16 are provided to demonstrate canonical and grand-canonical averaging of ¹³⁹La NMR spectra with stress-volume correction (see `x250` and `x750`).

- **example06**: Li/Mg substitution coupled with an H vacancy in a 2×2×2 supercell of rutile MgH₂ (space group P4₂/mnm, #136). The charge-neutral defect pair — Li⁺ on an Mg²⁺ site plus one H⁻ vacancy — is enumerated as a **multi-species** substitution: `sptarget: 1 2` with `nsubs` = 1 (Li) on line 1 and 1 (%H vacancy) on line 2. The 2×2×2 supercell contains 16 Mg and 32 H sites, giving 512 total Li–vacancy arrangements of which **9 are inequivalent** under the full tetragonal symmetry. VASP input files (FILER=11) are generated; the user must supply `INCAR`, `KPOINTS` and `POTCAR`. See: Smith, K.C., Fisher, T.S., Waghmare, U.V., Grau-Crespo, R., *Dopant-vacancy binding effects in Li-doped magnesium hydride*, Phys. Rev. B **82**, 134109 (2010).

- **example07**: Fe vacancies in maghemite (γ-Fe₂O₃) starting from the P4₃32 cubic structure (space group 212, a=8.344 Å). A 1×1×3 supercell has 12 vacant Fe sites; SOD fills 4 of them with Fe, leaving 8 vacancies to reach the γ-Fe₂O₃ stoichiometry. 29 inequivalent configurations out of C(12,4)=495. Demonstrates the `%NAME` vacancy syntax (`%Fe` in `newsymbol`). GULP input files (FILER=1) generated using the Catlow library (`catlow.lib`).

- **example08**: Methylammonium (MA = CH₃NH₃⁺) substitution in a 4×4×4 supercell of cubic CsPbI₃ perovskite (Pm-3m, 320 atoms). 2 of the 64 Cs A-sites are replaced by MA molecules, yielding 9 inequivalent configurations. Demonstrates the `@NAME` molecule syntax (`@MA` in `newsymbol`) with `MA.xyz` geometry. LAMMPS data files (FILER=2) are generated; each `conf.data` has 334 atoms (62 Cs + 64 Pb + 192 I + 2 C + 2 N + 12 H).

- **example09**: Simultaneous Mg/La substitution (2 sites) and O vacancy (1 site) in a 2×2×2 supercell of LaFeO₃ perovskite — a **multi-species** example. Two target sites are specified (`sptarget: 1 3`), with `nsubs` given as `2` (line 1) and `1` (line 2). Configurations are enumerated jointly under the full crystal symmetry and written to `n02_01/`. FILER=-1.

- **example10**: $Ti_{50}Zr_{25}Nb_{25}$ ($Ti_2ZrNb$) alloy with biomedical interest — a **multi-nary** example. Two species, Zr (4 atoms) and Nb (4 atoms), substitute Ti in a 2×2×2 supercell of the BCC structure (16 atoms total). No output files are requested (FILER=-1).

- **example11**: Equimolar NiCoFeCr Cantor subsystem alloy in a 2×2×2 supercell of the FCC primitive cell (8 atoms total, 2 of each species) — a **multi-nary** example with three new species. The FCC primitive cell (a=b=c=2.491 Å, α=β=γ=60°) is used with operators from `225_primitive.sgo`. `nsubs: 2 2 2` places 2 Co, 2 Fe and 2 Cr on Ni sites, leaving 2 Ni. Out of 2520 total arrangements (8!/(2!2!2!2!)), 23 are inequivalent under the full FCC symmetry. GULP input files (FILER=1) are generated using an OpenKIM EAM potential (`kim_model` directive).

- **example12**: Complex perovskite La₀.₇₅Sr₀.₂₅Mn₀.₂₅Fe₀.₇₅O₃ in a 2×2×2 supercell (space group Pm-3m) — a **multi-species** example. Two target sites are substituted simultaneously: 2 Sr replace La on 8 La-sites (binary, site 1) and 2 Mn replace Fe on 8 Fe-sites (binary, site 2). `nsubs` given as `2` (line 1) and `2` (line 2). Out of 784 total joint arrangements, 13 are inequivalent under the full cubic symmetry. FILER=-1.

- **example13**: La₁₋ₓSrₓFe₁₋ᵧMnᵧO₃₋ᵤ in a 2×2×2 supercell (space group Pm-3m) — a **multi-species** example with three target sites. One Sr replaces La (8 La-sites), one Mn replaces Fe (8 Fe-sites), and one O vacancy is created (24 O-sites), all enumerated simultaneously. `sptarget: 1 2 3`, `nsubs` on three lines: `1`, `1`, `1`. Vacancy on the O site uses the `%O` syntax. 6 inequivalent configurations from 1536 total. FILER=-1.

- **example14**: La₁₋ₓ₋ᵧSrₓBaᵧMnᵤFe₁₋ᵤO₃ in a 2×2×2 supercell (space group Pm-3m) — the first **multi-target multi-nary** example, combining multi-site and multi-nary substitution. Target 1 (La site, 8 atoms): ternary disorder with 1 Sr + 1 Ba (2 new species). Target 2 (Fe site, 8 atoms): binary disorder with 1 Mn. `nsubs` on two lines: `1 1` (line 1) and `1` (line 2). 3 inequivalent configurations from 448 total. FILER=-1.

## Molecules (@NAME) and vacancies (%NAME)

Two special prefixes extend the `newsymbol` field in INSOD beyond simple atomic substitution:

- **`@NAME`** — molecule: SOD reads `NAME.xyz` from the working directory (standard XYZ format: natoms, comment, then symbol x y z per line in Ångström), computes the centre of mass, and places the molecule at the substituted site with a uniformly random orientation. Each site gets an independent rotation. All output formats (CIF, GULP, VASP, CASTEP, QE, LAMMPS) expand the molecule into its individual atoms.

- **`%NAME`** — vacancy: the atom at the substituted site is simply omitted from all output files. `NAME` is informational only (e.g. `%Fe` or `%FeB`). All output formats including LAMMPS support vacancies.

Both can appear simultaneously in `newsymbol(1:2)`, and multiple molecule types can be used (up to 10 types). The `newsymbol` field accepts up to 5 characters (prefix + 4-character name), e.g. `@MA`, `@FA`, `@CO2`, `%FeB`.

## Compiling & installing SOD

- Download the file sod(version).tar.gz (e.g. sod0.71.tar.gz) and copy to a directory, say ROOTSOD:
 
```bash
tar xzvf sod(version).tar.gz
```

- Make compile all the executables into the **bin** folder:
 
```bash 
make all
```

- Add ROOTSOD/sod(version)/bin to your executables path 

```bash 
# add the bin folder to the executables path in your .bashrc file
export PATH=$PATH:ROOTSOD/sod(version)/bin
```

## Running SOD

- We recommend to create a new folder (say MAINFOLDER) for each sod project. This will be referred to as the working directory.

- In MAINFOLDER, you must create a file named *INSOD* which contains all the information for running the combinatorics part of the program. Use the *INSOD* file given in one of the examples as a template. The file is self-explanatory. The format of this file is rigid, so keep the same number of blank lines.

- The `nsubs` field in INSOD controls how many atoms of each new species are placed at the target site(s). Several formats are supported:
  - **Fixed count** (e.g. `4`): enumerate configurations with exactly 4 substitutions.
  - **Range** (e.g. `1:8`): SOD loops over all integer values from 1 to 8 in sequence, creating one `nXX/` folder per concentration. Only valid when a single target site and a single new species are specified.
  - **Multi-nary** (e.g. `1 2` or `2 2 2`): place the specified numbers of each new species simultaneously on a single target site. Up to 3 new species are supported (quaternary disorder, e.g. `nsubs: 2 2 2` for NiCoFeCr). Extension to higher orders (quinary and beyond) is planned but not yet implemented.
  - **Multi-species** (e.g. two lines `2` then `1`): one line per target site, in the order listed in `sptarget`. The example places 2 substitutions on the first target site and 1 on the second, enumerating all joint configurations simultaneously under the full crystal symmetry.
  - **Multi-target multi-nary** (e.g. line 1: `1 1`, line 2: `1`): combines multi-species and multi-nary — each target site can independently have multiple new species. Each line contains space-separated counts for the new species on that site, and enumerates all joint configurations simultaneously.

- In MAINFOLDER, you must also include a file named SGO with the matrix-vector representations of the symmetry operators. First check if your space group is included in the ROOTSOD/sod(version)/sgo library; if this is the case, just copy the file into your working directory, under the name SGO:

```bash
cp ROOTSOD/sod(version)/sgo ./SGO
```

otherwise you have to create the file using the International Tables of Crystallography, or from the website of the Bilbao Crystallographic Server <www.cryst.ehy.es>. The first three numbers in each line are one row of the operator matrix and the fourth number is the component of the operator translation vector.

The FILER value in INSOD controls which type of calculation input files are generated:

| FILER | Code | Template file needed |
|------:|------|----------------------|
| -1 | (none — no files generated) | — |
| 0 | CIF (P1, one per config) | — |
| 1 | GULP | `template_input.gin` |
| 2 | LAMMPS | `template_in.lammps` |
| 11 | VASP | — (no template; SOD generates `POSCAR` directly) |
| 12 | CASTEP | `template_castep.cell` |
| 13 | Quantum ESPRESSO | `template_pw.in` |

### Template file philosophy

For each supported calculator (FILER > 0, except for FILER = 11 - VASP), the user provides a single template file in the working directory. The template should look as much as possible like the real input file that will later be run, with two kinds of SOD-specific additions:

- **`@configuration_structure@`** — a token on a line by itself, marking where SOD inserts the configuration-specific structure block (cell parameters and atomic positions). This token is mandatory and must appear exactly once.
- **`@configuration_number@`** — an optional token that can appear anywhere (e.g. in titles or output file names) and is replaced by the zero-padded configuration index.

All calculator-specific settings — force-field parameters, k-points, convergence thresholds, output directives — remain in the template unchanged. SOD only inserts the structure.

For GULP, CASTEP and QE, template files follow the naming convention `template_<real_output_name>`, so the generated file name is immediately obvious:

| Input file | Generated file in each `cYY/` |
|------------|-------------------------------|
| `template_input.gin` | `input.gin` |
| `template_in.lammps` | `in.lammps` (+ `conf.data`) |
| `template_castep.cell` | `castep.cell` |
| `template_pw.in` | `pw.in` |

### Calculator-specific notes

**GULP** (`template_input.gin`): The template must contain `@configuration_structure@` on a line by itself; SOD replaces it with the `cell`, `frac` and atom-type block. If the template contains a `library` directive, SOD reads the library, extracts only the entries needed for the atom types in the problem, and writes a reduced copy of the library into each configuration folder so that each folder is self-contained. If the force field is defined inline in the template, no library file is copied. If any SOD species label differs from the corresponding GULP atom type, add mapping comment lines to the template of the form `# sod_type_map <SOD_species> <GULP_type>`.

**LAMMPS** (`template_in.lammps`): The user provides `template_in.lammps`, a LAMMPS run script containing `read_data conf.data` and all force-field settings. SOD copies it into each `cYY/` folder with only two modifications: `# sod_type_map` lines are stripped (they are SOD directives, not valid LAMMPS syntax), and the optional token `@configuration_number@` is replaced if present. Everything else — force-field coefficients, masses, charges, run commands — is copied verbatim. SOD also generates a `conf.data` file in each folder containing the configuration-specific structural information.

To generate `conf.data`, SOD needs to know which LAMMPS numeric type corresponds to each SOD species. **SOD cannot infer this automatically.** The user must provide explicit mapping comment lines in `in.lammps` of the form:

```
# sod_type_map <SOD_species> <role> type=<N> [bond_type=<M>]
```

where `<role>` is `core` or `shell`. One `core` line is required for every SOD species. A `shell` line (with `bond_type=<M>`) is additionally required for any species represented as a core–shell pair. SOD will stop with an error if any species has no mapping.

**VASP**: No template is needed. SOD generates `POSCAR` into each `cYY/` folder. The user must supply `INCAR`, `KPOINTS` and `POTCAR` in each folder separately.

**CASTEP** (`template_castep.cell`): The template is a normal `.cell` file. SOD replaces `@configuration_structure@` with the `%BLOCK LATTICE_CART` and `%BLOCK POSITIONS_FRAC` blocks. The user must supply the `.param` file separately (it is not managed by SOD).

**Quantum ESPRESSO** (`template_pw.in`): The template is a normal `pw.x` input. SOD assumes `ibrav = 0` and replaces `@configuration_structure@` with `CELL_PARAMETERS` and `ATOMIC_POSITIONS` blocks. The user defines `ATOMIC_SPECIES` in the template; SOD does not modify it. Pseudopotential files are not copied by SOD.

## Directory structure and script calling levels

SOD uses the following directory hierarchy:

```
MAINFOLDER/               ← working directory: INSOD, SGO, template files, OUTSOD, EQMATRIX
  n01/              ← one folder per substitution level (zero-padded)
    OUTSOD
    ENERGIES
    c01/  c02/ ...  ← one folder per inequivalent configuration
  n02/  n04/ ...
  x250/             ← grand-canonical working folder (user-created, one per composition)
    INGC
  n04/spbe0/        ← SPBE subfolder (user-created)
  n04/spbe1/        ← SPBE subfolder (user-created)
```

The table below specifies from which directory each post-processing script should be called. Scripts marked **MAINFOLDER/ or nXX/** detect their calling level automatically: run from MAINFOLDER/ to process all substitution levels at once, or from a specific `nXX/` folder to process only that level.

| Script | Call from | What it does |
|---|---|---|
| `sod_comb.sh` | MAINFOLDER/ | Runs combinatorics and generates calculator input files |
| `sod_gulp_ener.sh` | MAINFOLDER/ or nXX/ | Extracts final energies from GULP `output.gout` files |
| `sod_vasp_ener.sh` | MAINFOLDER/ or nXX/ | Extracts final energies from VASP `OUTCAR` files |
| `sod_castep_ener.sh` | MAINFOLDER/ or nXX/ | Extracts final energies from CASTEP `castep.castep` files |
| `sod_qe_ener.sh` | MAINFOLDER/ or nXX/ | Extracts final energies from QE `pw.out` files (converts Ry→eV) |
| `sod_gulp_free.sh` | MAINFOLDER/ or nXX/ | Extracts vibrational free energies from GULP output |
| `sod_gulp_single_ener.sh` | MAINFOLDER/ or nXX/ | Extracts single-point energies from GULP output |
| `sod_gulp_cell.sh` | MAINFOLDER/ or nXX/ | Extracts cell parameters from GULP output → `CELL` file |
| `sod_vasp_cell.sh` | MAINFOLDER/ or nXX/ | Extracts cell parameters from VASP `CONTCAR` files → `CELL` file |
| `sod_vasp_mag.sh` | MAINFOLDER/ or nXX/ | Extracts magnetic moments from VASP `OUTCAR` files |
| `sod_stat.sh` | MAINFOLDER/ or nXX/ | Runs canonical statistical mechanics (`statsod`) |
| `sod_gcstat.sh` | x???/ | Runs grand-canonical statistical mechanics (`gcstatsod`) |
| `sod_spbe0.sh` | nXX/spbe0/ | SPBE energy extrapolation from the x=0 end |
| `sod_spbe1.sh` | nXX/spbe1/ | SPBE energy extrapolation from the x=1 end |

- To run the combinatorics program, just type:


```bash
sod_comb.sh
```

## Output of sod_comb.sh 

- When the programme finishes, it writes to the standard output the total number of configurations and the number of independent configurations according to the crystal symmetry, plus some other useful information.

- It writes the data file *OUTSOD*, which contains information on the independent configurations (one line for each configuration). The first number is the index of the configuration, the second is its degeneracy, and the next numbers are the substitution sites.

- It also writes the file *EQMATRIX*, which gives information about  how each supercell operator transforms each atom position. 

- For all calculators, SOD creates a folder `nXX/` (where XX is the zero-padded number of substitutions, e.g. `n04`) and inside it one folder `cYY/` per configuration (e.g. `c1`, `c01`, `c001` depending on total count). Each `cYY/` folder contains the complete input for that configuration. For GULP, LAMMPS, CASTEP, QE, and VASP, a `job_sender` script is written in the working directory to run all configurations in sequence. For CIF (FILER=0), only the `configuration.cif` files are written; no `job_sender` is created. This naming convention is used by the other SOD scripts for statistics and energy extrapolation.


## Configurational averages and thermodynamics:

In order to calculate configurational averages and obtain thermodynamic quantities, you need to execute the script:

```bash
sod_stat.sh
```

which requires 4 input files:

- *OUTSOD*, which was the output from sod_comb

- *TEMPERATURES*, a list of temperatures for the Boltzmann statistics, in one column, e.g.:

```bash
300
600
1000
1250
1500
1750
2000
```

(if the TEMPERATURES file does not exist, sod_stat calculates thermodynamic properties at T=1K, 300K, 1000K and in the limit of a very high temperature). 

- *ENERGIES*, which contains (in one column) the energies of all the configurations, in the same order that they were generated by SOD (like in the OUTSOD file). There are some scripts in ROOTSOD/sod(version)/bin/  to help you do this:

   1. If you are using GULP, the script  ```sod_gulp_ener.sh``` will extract all the energies, assuming the output file is named `output.gout` in each `cYY/` folder. If you have calculated vibrational free energies for each configuration, ```sod_gulp_free.sh``` will extract these.

   2. If you are using VASP, the script ```sod_vasp_ener.sh``` will extract all the energies, assuming the output file is named `OUTCAR` in each `cYY/` folder.

   3. If you are using CASTEP, the script ```sod_castep_ener.sh``` will extract all the energies, assuming the output file is named `castep.castep` in each `cYY/` folder.

   4. If you are using Quantum ESPRESSO, the script ```sod_qe_ener.sh``` will extract all the energies, assuming the output file is named `pw.out` in each `cYY/` folder.

- *DATA*, which contains *ncol* colums of data to average. The first line contains just the number *ncol* of columns to read. For example:

```bash
2
34.5   4.34
37.7   4.35
35.6   4.38
38.8   4.41
```

The data can be cell lenghts, or volumes (please see SOD papers for strategies on how to obtain average cell parameters) or any other observable obtained from the calculations. Scripts like ```sod_vasp_cell.sh``` can help you do this, please edit carefully before using them.

```sod_stat.sh``` will generate two files: probabilities.dat and thermodynamics.dat, whose content is self-explanatory.


Important note: While configurational averages (e.g. of cell parameters and enthalpies) tend to converge very quickly with supercell size, entropies and free energies, which are not defined by averaging, converge very slowly with supercell size, and are generally in large error when using the SOD method. We therefore do not recommend using SOD for the calculation of entropies and free energies, unless appropriate correcting procedures have been used.


## Grand-canonical analysis 

From version 0.51, it is possible to do statistics in a grand-canonical ensemble, i.e. including results from supercells with different compositions. Please see example04 (perovskites) and example05 (pyrochlore).

We recommend to create a file with name x??? at the same level as the n?? files. For example x250 is used for a grandcanonical analyis at composition x=0.250. 

 The *OUTSOD_00*, *OUTSOD_01*, *OUTSOD_02*, *OUTSOD_03*, and *OUTSOD_04* files are the *OUTSOD* files for 0, 1, 2, 3, and 4 substitutions, respectively. You also need the *ENERGIES_00* ... *ENERGIES_04* files there. Optionally, you can add *DATA_00*, ..., *DATA_04*. The *TEMPERATURES* file can also be provided (optional, as for the canonical statistics). 

In order to do the grand-canonical analysis, you need the grand-canonical input file *INGC*, which has a very simple structure. For example, to set the chemical potential, the first few lines of the file look like this:

```bash
# nsubsmin nsubsmax
0   4
# Specify x or mu, and provide its value
mu -0.5
```

But it is possible to specify the composition (fraction x=nsubs/npos of sites that are substituted) and the chemical potential will be calculated automatically for each temperature. The chemical potential is obtained from a solution to a polynomial equation, using a simple bisection method. To do this the INGC file should look like:

```bash
# nsubsmin nsubsmax
0   4
# Specify x or mu, and provide its value
x 0.25
```

In this example, x=0.25 corresponds to 2 substitutions in the canonical example, but in the grand-canonical analysis of the example, all compositions from 0 to 4 substitutions are included. 

To run the grand-canonical analysis, type:

```sod_gcstat.sh```

If the naming convention *n??* for the different compositions was followed, and the *x???* file is at the same level of those, the script will copy all the necessary files automatically, so you only need the *INGC* file. The analysis produces a probabilities.dat and a thermodynamics.dat file as in the canonical analysis.   

Finally, it is possible to make a "stress-volume" correction to the energies in the grand-canonical configurational ensemble. This correction is a simple way to account for the fact that if a cell with number of substitutions n (different from xN) contributes to a grand-canonical ensemble representing composition x, there is an additional stress-related energy cost. This is due to the difference in equilibrium volumes at compositions n/N and x. In a physically grounded approximation based on the second-order Birch-Murnaghan (BM2) equation of state, the energy of the stress-volume correction ($ESVC$) is

$ESVC(n,x) = \frac{9}{8} V(x) B(x) \left[ \left(\frac{V(x)}{V\left(\frac{n}{N}\right)}\right)^{2/3} - 1 \right]^2$

where $B(x)$ is the bulk modulus and $V(x)$ is the equilibrium volume at composition $x$, and $V(n/N)$ is the equilibrium volume at the composition of the contributing supercell. This expression is positive definite for any volume departure, reduces correctly to $\frac{1}{2} B(x) (\Delta V)^2 / V(x)$ in the small-strain limit, and requires only $B_0$ as a parameter (implicitly assuming $B'=4$).

It is possible to introduce this correction in the grand-canonical analysis by adding the following lines to INGC (see example5): 

```bash
# Stress corrections flag (lambda=0: no correction; lambda=1: bulk moduli-based correction)
1
# Parameters for volume variation with x: v0, v1, bv (Angstrom^3)
1302.567820  1286.687504  0
# Parameters for bulk modulus variation with x: bm0, bm1, bb (GPa)
150 150 0
```
This setting leads to a simple linear interpolation of the equilibrium volumes and bulk moduli between the solid solution endmembers. A quadratic interpolation is also possible by using non-zero values of the bowing parameters $bv$ and/or $bb$. This functionality has not been well tested yet. If interested in using this correction scheme, please contact the SOD developers for further information. 


## Averaging spectra

Both in the canonical and in the grand-canonical analysis, it is possible to evaluate averages of spectra in the corresponding configurational ensemble. 

Often, what is originally calculated for a given configuration is a list of peaks. In that case, running the "peaks2spectra" code, with the script: 

```sod_p2s.sh```

will generate the broadened spectra that will be averaged. This requires two input files (with fixed names):

- *PEAKS* contains a list of peaks in each line; each line represents a different configuration (a different spectrum is generated for each configuration)
- *INP2S* contains other info needed to generate the spectra, e.g. xmin, xmax, broadening(sigma), etc. See in example5, where we use:
```
 # nconfs
 22
 # peaks
 4
 # xmin
 -680
 # xmax
 -620
 # npoints
 601
 # broadening (sigma)
 0.85
```

That generates two output files:

- *SPECTRA*, where each line contains the generated intensity values in the x grid, for each configuration.
- *XSPEC*, which contains the list of x values at which the intensities are given.

If these files are present within the n??/ folders when running the statistics codes (either ```sod_stat.sh``` or ```sod_gcstat.sh```), then a file with name ```ave_spectra.dat``` will be created with the configurational averages of the spectra at different temperatures.  


## Extrapolating energies from low to high concentrations

SOD has a simple pair-based extrapolation (SPBE) method, which uses the energies from *n*= 0, 1 and 2 substitutions to predict the energies for *n*>2 (equation 1 in [Arce-Molina et al. PCCP 2018, 20, 18047-18055](https://pubs.rsc.org/en/content/articlehtml/2018/cp/c8cp01369a)).   

In order to run this program, you will need the following files:

- *EQMATRIX* obtained from running sod_comb at any composition
- *OUTSOD* for *n* substitutions
- *ENERGIES0*: a file containing a single number, which is the energy for *n*=0
- *ENERGIES1* and *OUTSOD1*: the *ENERGIES* and *OUTSOD* files for *n*=1 
- *ENERGIES2* and *OUTSOD2*: the *ENERGIES* and *OUTSOD* files for *n*=2 
- *INSPBE* file if you want to introduce some rescaling parameters (optional, see below)

If all the above files are present in a folder, you can run the spbe module by running the ```spbesod``` executable. 

However,  the easiest way to run the spbe module is like this:

- Make sure that the folders n00, n01 and n02 contain an ENERGIES and an OUTSOD file each (OUTSOD is not necessary for n00)
- If you want to use spbe, say, for n=3, first run ```sod_comb.sh``` for n=3 substitutions (this automatically creates the folder n03), and create a subfolder within it, say n03/spbe/
- From the n03/spbe folder, just run the script ```sod_spbe0.sh```, which will copy the relevant input files into the current folder and will call ```spbesod```
- It is also possible to run the spbe program using data from the other end of the solid solution (i.e. *x*=1). In that case, run the script ```sod_spbe1.sh```, which will copy the files from the folders with *N*, *N*-1, *N*-2 substitutions, will "invert" the OUTSOD files as needed, and call ```spbesod```. 

Finally, it is possible to introduce some rescaling in the first-order and second-order terms to improve the match with a reference set of calculations. You need to give two reference energies in the INSPBE file. The recommended procedure is to run spbe first without rescaling, pick the minimum-energy and maximum-energy configurations (they are identified at the end of the OUTSPBE file) and run them with DFT (or whatever method provides the reference/target values), then input these two values as reference energies in INSPBE, and run the sod_spbe0.sh script again. See example04 (inside n04/spbe0), where the reference energies for configurations 1 and 6 are given as input. 


## Citing SOD

If you use SOD in your research work, please include a citation to this article:

*Grau-Crespo, R., Hamad, S., Catlow, C. R. A., & De Leeuw, N. H. (2007). Symmetry-adapted configurational modelling of fractional site occupancy in solids. Journal of Physics: Condensed Matter, 19(25), 256201.*
[Original Paper](http://iopscience.iop.org/article/10.1088/0953-8984/19/25/256201/meta) 


Happy SODing!!!

Ricardo Grau-Crespo (r.grau-crespo@qmul.ac.uk) and Said Hamad (said@upo.es)
