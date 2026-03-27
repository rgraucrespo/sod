# SOD 0.61 - Notes for users

SOD (standing for Site-Occupation Disorder) is a package of tools for the computer modelling of periodic systems with site disorder, using the supercell ensemble method. 

The package is distributed under the [GPL licence](https://github.com/ypriverol/sod/blob/master/LICENSE.md). 

You can find below the essential info needed to use SOD. Please note that SOD authors can give only limited support to users.


## Functionalities

- Identification of all inequivalent configurations of site substitutions in an arbitrary supercell of an  initial target structure with any space group symmetry.
- Calculation of the degeneracies of configurations.
- Generation of input files for codes like GULP and VASP.
- Simple extrapolation of energies from low to high concentrations within a supercell.
- Statistical mechanics processing of output using either canonical or grand-canonical ensembles.


## Content of the folders

- sod(version)/src contains the source files.
- sod(version)/sgo is a library of space group operators (e.g. 131.sgo contains the operators of the space group 131).
- sod(version)/bin contains the executables. Compile for your platform using `make all`.
- sod(version)/examples contains six examples, covering rocksalt (MgO), inverse spinel (magnetite), rutile, perovskite (LaFeO3), pyrochlore (La2Zr2O7) and wurtzite (ZnO) structures.

## Examples

- **example1**: Ca/Mg substitutions in a 2×2×2 supercell of rocksalt MgO. VASP input files created.

- **example2**: Al/Fe substitution in magnetite Fe₃O₄ (cubic unit cell of the spinel structure). The substitution occurs in both tetrahedral and octahedral sites, allowing the study of the distribution of substitutions over the two different sites. VASP input files created.

- **example3**: Fe/Sb disorder in a 2×2×2 supercell of rutile FeSbO₄. Starting from a generic MO₂ composition, M(IV) is replaced by equal parts of Fe(III) and Sb(V). VASP input files created. See [Grau-Crespo et al., Chemistry of Materials (2004)](https://pubs.acs.org/doi/abs/10.1021/cm035271y).

- **example4**: Al/Fe substitution in a 2×2×2 supercell of LaFeO₃ perovskite, using GULP. Also contains examples of energy extrapolation using SPBE0 and SPBE1 (see folders inside n04) and of grand-canonical statistics (x0.25 folder).

- **example5**: Averaging NMR spectra in the canonical or grand-canonical configurational space of the La₂Zr₂O₇–La₂Sn₂O₇ solid solution. CASTEP input files are generated (FILER=12). The DFT calculation files are excluded to save space; only the statistical analysis files (ENERGIES, PEAKS, SPECTRA, etc.) are included in the n00–n16 folders.

- **example6**: Mg/Zn substitution in a 2×2×1 supercell of wurtzite ZnO (space group 186, P6₃mc). The 8 cation sites with 2 Mg substitutions (Mg₀.₂₅Zn₀.₇₅O) reduce to 3 inequivalent configurations under the hexagonal symmetry. Quantum ESPRESSO input files are generated (FILER=13). Includes `top.qe` and `bottom.qe` templates; replace `Mg.upf`, `Zn.upf`, `O.upf` with your pseudopotential files before running.

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

- If you want to generate GULP input files for all the independent configurations found by SOD, in addition to setting FILER=1 in the INSOD file, you must provide two files in the working directory:

top.gulp contains the heading of the gulp input file (until the keyword cell).

bottom.gulp contains the tail of the gulp input file (everything after the list of coordinates, including species, potentials, etc).

Optionally, you can request that GULP write output structure files by setting the `genxtl` and `genarc` flags in INSOD (see the GULP examples for the format). Setting `genxtl 1` adds an `output xtl` directive, and `genarc 1` adds an `output arc` directive to each GULP input file.

- To run the combinatorics program, just type:


```bash
sod_comb.sh
```

## Output of sod_comb.sh 

- When the programme finishes, it writes to the standard output the total number of configurations and the number of independent configurations according to the crystal symmetry, plus some other useful information.

- It writes the data file *OUTSOD*, which contains information on the independent configurations (one line for each configuration). The first number is the index of the configuration, the second is its degeneracy, and the next numbers are the substitution sites.

- It also writes the file *EQMATRIX*, which gives information about  how each supercell operator transforms each atom position. 

- A folder named *nXX* is generated (where XX is the zero-padded number of substitutions, e.g. *n04*), which contains the calculation input files, a copy of *OUTSOD*, and a script that sends the jobs. This naming convention is used by the other SOD scripts for statistics and energy extrapolation.


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

From version 0.44, we have implemented a simple pair-based extrapolation (SPBE) method, which uses the energies from *n*= 0, 1 and 2 substitutions to predict the energies for *n*>2 (equation 1 in [Arce-Molina et al. PCCP 2018, 20, 18047-18055](https://pubs.rsc.org/en/content/articlehtml/2018/cp/c8cp01369a)).   

In order to run this program, you will need the following files:

- *EQMATRIX* obtained from running sod_comb at any composition
- *OUTSOD* for *n* substitutions
- *ENERGIES0*: a file containing a single number, which is the energy for *n*=0
- *ENERGIES1* and *OUTSOD1*: the *ENERGIES* and *OUTSOD* files for *n*=1 
- *ENERGIES2* and *OUTSOD2*: the *ENERGIES* and *OUTSOD* files for *n*=2 
- *INSPBE* file if you want to introduce some rescaling parameters (optional, see below)

If all the above files are present in a folder, you can run the spbe module by running the ```spbesod``` executable. 

However,  the easiest way to run the spbe module is like this:

- Use the names n00 n01 n02 n03 etc for the folders containing the calculations for 0, 1, 2, 3... substitutions. 
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
