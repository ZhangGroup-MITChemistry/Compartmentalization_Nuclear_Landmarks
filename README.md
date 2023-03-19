# Compartmentalization with Nuclear Landmarks	

LAMMPS input files and source codes to execute chromatin simulations in the presence of particle based description of lamina, nucleolus and speckles. The cell type is HFF (Human Foreskin Fibroblasts) combined with the SPIN (Spatial Position Inference of the Nuclear genome, Wang et al., <i>Genome Biology</i>, 2021) state data for the same cell type. 

The files in this repository are useful for reproducing calculations reported in the manuscript "[Compartmentalization with nuclear landmarks yields random, yet precise, genome organization](https://www.cell.com/biophysj/pdf/S0006-3495(23)00156-X.pdf)" by Kamat et al., <i> Biophysical Journal</i>, 2023 and [Preprint](https://www.biorxiv.org/content/10.1101/2021.11.12.468401v1.abstract).
<hr>
The /src folder contains files that must be added in the src folder of LAMMPS before building LAMMPS.
The /input_files contains all the necessary LAMMPS input files and dependency files to perform Molecular Simulations. 

The atom type naming convention we used is
- Types 1-4 are chromatin
- Type 5 is nucleolus 
- Type 6 is lamina 
- Type 7 is speckles
<hr>
If the code or data is useful please cite Kamat et al., <i>Biophysical Journal</i>, 2023. Additional simulation data and scripts is available upon request. Please reach out to Kartik Kamat (kartikk(at)mit.edu) with any questions.
