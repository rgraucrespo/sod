# SOD 0.61 - Notes for users

SOD (standing for Site-Occupation Disorder) is a package of tools for the computer modelling of periodic systems with site disorder, using the supercell ensemble method. 

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
- sod(version)/examples contains six examples, covering rocksalt (MgO), inverse spinel (magnetite), rutile, perovskite (LaFeO3), pyrochlore (La2Zr2O7) and wurtzite (ZnO) structures.

## Examples

The **example1_\*** series all use the same physical system — Ni/Mg substitutions in a 2×2×2 supercell of rocksalt MgO (space group Fm-3m) — and illustrate input generation for five different calculators using the same INSOD file (differing only in FILER):

- **example1_gulp**: 1 Ni substitution (nsubs=1). GULP (FILER=1) with `catlow.lib` (Buckingham potentials, core-shell model for Ni and O). Demonstrates the `library` directive and `sod_type_map` (SOD species `O` maps to GULP type `O2`).

- **example1_lammps**: 4 Ni substitutions (nsubs=4). LAMMPS (FILER=2). Demonstrates `template_in.lammps` with `sod_type_map` lines and core-shell representation in the generated `conf.data`.

- **example1_vasp**: 4 Ni substitutions. VASP (FILER=11). SOD generates `POSCAR` into each `cYY/` folder; the user must supply `INCAR`, `KPOINTS` and `POTCAR`.

- **example1_castep**: 4 Ni substitutions. CASTEP (FILER=12). Demonstrates `template_castep.cell`; the user must supply the `.param` file.

- **example1_qe**: 4 Ni substitutions. Quantum ESPRESSO (FILER=13). Demonstrates `template_pw.in`; replace `Ni.upf`, `Mg.upf`, `O.upf` with your pseudopotential files before running.

- **example2**: Al/Fe substitution in magnetite Fe₃O₄ (1×1×1 cubic spinel unit cell, space group Fd-3m). 4 Al substitutions across 24 mixed Fe sites (8 tetrahedral + 16 octahedral) give 99 inequivalent configurations. VASP input files (FILER=11) pre-generated in `n04/`. The two structurally distinct Fe environments make this a non-trivial site-disorder problem.

- **example3**: Fe/Sb disorder in a 2×2×2 supercell of rutile FeSbO₄ (space group P4₂/mnm). 8 M(IV) sites replaced by equal numbers of Fe(III) and Sb(V) give 182 inequivalent configurations out of 12870 total. VASP input files (FILER=11) pre-generated in `n08/`. See [Grau-Crespo et al., Chemistry of Materials (2004)](https://pubs.acs.org/doi/abs/10.1021/cm035271y).

- **example4**: Al/Fe substitution in a 2×2×2 supercell of LaFeO₃ perovskite (space group Pnma). 4 Al in 8 Fe sites give 6 inequivalent configurations. GULP input files (FILER=1) using the Bush et al. Buckingham potentials with core-shell model, defined inline in `template_input.gin`. Also includes SPBE energy extrapolation (subfolders `n04/spbe0` and `n04/spbe1`) and a grand-canonical statistics example (`x250`, composition x=0.25).

- **example5**: Sn/Zr substitution in the La₂Zr₂O₇–La₂Sn₂O₇ solid solution (pyrochlore, space group Fd-3m). 2 Sn substitutions in 16 Zr sites give 5 inequivalent configurations. CASTEP input files (FILER=12) are generated via `template_castep.cell`; DFT outputs are not included to save space. Demonstrates canonical and grand-canonical averaging of ¹³⁹La NMR spectra across compositions n00–n16, with optional stress-volume correction (see `x125` and `x250`).

- **example6**: Mg/Zn substitution in a 2×2×1 supercell of wurtzite ZnO (space group P6₃mc). 2 Mg substitutions in 8 Zn sites reduce to 3 inequivalent configurations under the hexagonal symmetry. Quantum ESPRESSO input files (FILER=13) generated via `template_pw.in`. Replace `Mg.upf`, `Zn.upf`, `O.upf` with your pseudopotential files before running.

## Compiling & installing SOD

- Download the file sod(version).tar.gz (e.g. sod0.61.tar.gz) and copy to a directory, say ROOTSOD:
 
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

- We recommend to create a new folder (say FOLDER_NAME) for each sod project. This will be referred to as the working directory.

- In FOLDER_NAME, you must create a file named *INSOD* which contains all the information for running the combinatorics part of the program. Use the *INSOD* file given in one of the examples as a template. The file is self-explanatory. The format of this file is rigid, so keep the same number of blank lines.

- In FOLDER_NAME, you must also include a file named SGO with the matrix-vector representations of the symmetry operators. First check if your space group is included in the ROOTSOD/sod(version)/sgo library; if this is the case, just copy the file into your working directory, under the name SGO:

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

- To run the combinatorics program, just type:


```bash
sod_comb.sh
```

## Output of sod_comb.sh 

- When the programme finishes, it writes to the standard output the total number of configurations and the number of independent configurations according to the crystal symmetry, plus some other useful information.

- It writes the data file *OUTSOD*, which contains information on the independent configurations (one line for each configuration). The first number is the index of the configuration, the second is its degeneracy, and the next numbers are the substitution sites.

- It also writes the file *EQMATRIX*, which gives information about  how each supercell operator transforms each atom position. 

- For calculators that use template files (GULP, LAMMPS, CASTEP, QE), SOD creates a folder `nXX/` (where XX is the zero-padded number of substitutions, e.g. `n04`) and inside it one folder `cYY/` per configuration (e.g. `c1`, `c01`, `c001` depending on total count). Each `cYY/` folder contains the complete, ready-to-run input for that configuration. A `job_sender` script is written in the working directory to run all configurations in sequence. For VASP and CIF, a single `nXX/` folder is created containing one file per configuration. This naming convention is used by the other SOD scripts for statistics and energy extrapolation.


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

   1. If you are using GULP, the script  ```sod_gulp_ener.sh``` will extract all the energies, assuming all output files,  with extension .gout, are in the same folder. If you have calculated vibrational free energies for each configuration, ```sod_gulp_free.sh``` will extract these. 

   2. If you are using VASP, the script ```sod_vasp_ener.sh``` will extract all the energies, assuming you have separate folders for each configuration. 

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

From version 0.51, it is possible to do statistics in a grand-canonical ensemble, i.e. including results from supercells with different compositions. Please see example4 (perovskites) and example5 (pyrochlore).

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

Finally, it is possible to introduce some rescaling in the first-order and second-order terms to improve the match with a reference set of calculations. You need to give two reference energies in the INSPBE file. The recommended procedure is to run spbe first without rescaling, pick the minimum-energy and maximum-energy configurations (they are identified at the end of the OUTSPBE file) and run them with DFT (or whatever method provides the reference/target values), then input these two values as reference energies in INSPBE, and run the sod_spbe0.sh script again. See example4 (inside n04/spbe0), where the reference energies for configurations 1 and 6 are given as input. 


## Citing SOD

If you use SOD in your research work, please include a citation to this article:

*Grau-Crespo, R., Hamad, S., Catlow, C. R. A., & De Leeuw, N. H. (2007). Symmetry-adapted configurational modelling of fractional site occupancy in solids. Journal of Physics: Condensed Matter, 19(25), 256201.*
[Original Paper](http://iopscience.iop.org/article/10.1088/0953-8984/19/25/256201/meta) 


Happy SODing!!!

Ricardo Grau-Crespo (r.grau-crespo@qmul.ac.uk) and Said Hamad (said@upo.es)
