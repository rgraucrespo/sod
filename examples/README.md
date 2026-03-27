## Explanation of Examples

---

### **Example 1**: **Ca/Mg Substitutions in a 2x2x2 Supercell of Rocksalt MgO**
- **Description**: This example simulates the substitution of Ca and Mg in a 2x2x2 supercell of rocksalt MgO.
- **Objective**: Investigate the effects of Ca/Mg substitution in the MgO crystal structure.
- **Method**: VASP input files were created to perform DFT calculations for this substitution process.

---

### **Example 2**: **Al/Fe Substitution in Magnetite Fe₃O₄ (Cubic Unit Cell of the Spinel Structure)**
- **Description**: This example involves Al and Fe substitutions in the magnetite Fe₃O₄ structure, which has a cubic spinel structure.
- **Objective**: Study the distribution of Al/Fe substitutions over tetrahedral and octahedral sites within the spinel structure.
- **Method**: VASP input files were created for this substitution process, enabling the analysis of substitutions on both sites.

---

### **Example 3**: **Fe/Sb Disorder in a 2x2x2 Supercell of Rutile FeSbO₄**
- **Description**: This example models Fe/Sb disorder in a 2x2x2 supercell of rutile FeSbO₄.
- **Approach**: 
  - Due to the chemical incompatibility of starting with Fe₂O₄ or Sb₂O₄ in terms of oxidation states, the model begins with a generic MO₂ composition.
  - M(IV) is replaced with equal parts of Fe(III) and Sb(V) to simulate the disorder.
- **Reference**: For further details, see *Grau-Crespo et al.*, *Chemistry of Materials* (2004): [DOI link](https://pubs.acs.org/doi/abs/10.1021/cm035271y).
- **Method**: VASP input files were created to simulate the disorder in this system.

---

### **Example 4**: **Al/Fe Substitution in a 2x2x2 Supercell of LaFeO₃ Perovskite**
- **Description**: This example involves Al/Fe substitution in the LaFeO₃ perovskite structure.
- **Features**: 
  - Includes **energy extrapolation** using SPBE0 and SPBE1.
  - Uses **grand-canonical statistics** to model different configurations, with a specific example in the `x0.25` folder.
- **Method**: GULP is used to create input files for this Al/Fe substitution study in the perovskite structure.

---

### **Example 5**: **Averaging NMR Spectra in the Configurational Space of the La₂Zr₂O₇-La₂Sn₂O₇ Solid Solution**
- **Description**: This example demonstrates how to average NMR spectra in the configurational space of the La₂Zr₂O₇-La₂Sn₂O₇ solid solution.
- **Method**: 
  - NMR peaks were originally calculated using CASTEP.
  - Configurational statistics files are provided (excluding DFT files to save space).
- **Objective**: To study the NMR spectra by accounting for different configurations of the solid solution.
