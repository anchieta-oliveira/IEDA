# IEDA
IEDA is a tool for analyzing, quantifying, and visualizing intermolecular electronic density (IED) in biomolecular systems. The methodology uses quantum approaches derived from Mulliken overlap populations and orbital-based (OB), allowing the identification of regions where there is electron sharing between non-covalently bonded molecules, such as in hydrogen bonds. The goal of IEDA is to provide an electronic descriptor capable of assisting in the study of binding affinity, molecular recognition, and detailed analysis of intermolecular interactions in computational models and experimental data.
![](https://github.com/anchieta-oliveira/IED_/blob/main/doc/gallery/figure_1.png)

## Installation 
### Step 1: Clone the Repository
Use Git to clone the repository to your local machine:
```bash
git clone https://github.com/anchieta-oliveira/IEDA
```

### Step 2: Create a virtual environment dedicated to IEDA
Navigate to the cloned Git folder:
```bash
cd IEDA/
```
Create a new Conda environment:
```bash
conda create -n IEDA python=3.12
```
Activate environment conda IEDA:
```bash
conda activate IEDA
```

### Step 3: Then install the package
Install install the package
```bash
pip install .
```

### Step 4: Quickstart
Displaying IEDA usage options:
```bash
IEDA
```
Calculate IED matrix:
```bash
IEDA matrix --pdb for/your/pdbfile.pdb --qm for/your/qmfile.out --qm_sof orca 
```
Assemble the 3D map in a .pdb file with the IED per atom in the beta column, which can be viewed in any visualization program (pyMOL, VMD, Chimera, etc.):
```bash
IEDA map_3D --pdb for/your/pdbfile.pdb --ied for/your/IED_mulliken.npy --intrachain False
```
OBS: The "--intrachain" parameter only considers the IED between different chains, for example, the protein in chain A and the ligand in chain X. The map will only represent the IED between the protein and the ligand. 

Calculate the IED between two selections:
```bash
IEDA two_sel --pdb for/your/pdbfile.pdb --qm for/your/qmfile.out --qm_sof orca --sel_a "chain A" --sel_b "resname LIG" 
```
The IEDA selection pattern can be seen [here](https://anchieta-oliveira.github.io/IEDA/tutorial/selection_syntax.html).

### Step 4.2: Quickstart Coogle Colab
IEDA can be used within Jupyter Notebook and Google Colab environments. A simple usage example is available [here](https://colab.research.google.com/drive/1PvCzlCzPxuIH9JdD_qH3dX13WGOK9imU?usp=sharing).


For more information, please refer to the [documentation](https://anchieta-oliveira.github.io/IEDA/). 