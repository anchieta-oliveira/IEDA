# IEDA
IEDA is a tool for analyzing, quantifying, and visualizing intermolecular electronic density (IED) in biomolecular systems. The methodology uses quantum approaches derived from Mulliken overlap populations and atom-in-molecule (AIM) theory, allowing the identification of regions where there is electron sharing between non-covalently bonded molecules, such as in hydrogen bonds. The goal of IEDA is to provide an electronic descriptor capable of assisting in the study of binding affinity, molecular recognition, and detailed analysis of intermolecular interactions in computational models and experimental data.
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
conda env create -n IEDA python=3.12
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