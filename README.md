# Nuclear Landmarks	

LAMMPS input files and source codes to execute chromatin simulations in the presence of lamina, nucleolus and speckles. The cell type being used is HFF combined with the 
SPIN state analysis carried out for the same cell type. 

The files in this repository are useful for reproducing calculations reported in https://www.biorxiv.org/content/10.1101/2021.11.12.468401v1.abstract

The /src folder contains files that must be added in the src folder of LAMMPS before building LAMMPS.
The /input_files contains all the necessary LAMMPS input files and dependency files to perform Molecular Simulations. 

The atom type naming convention we used is
- Types 1-4 are chromatin
- Type 5 is nucleolus 
- Type 6 is lamina 
- Type 7 is speckles