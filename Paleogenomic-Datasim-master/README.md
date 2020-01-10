# Paleogenomic-Datasim
Simulation and analysis of paleogenomic data.

All of these scripts are written for use with qsub and SGE

Example of execution to simulate and phase:

* Samples with 25 generations of age
* 100 of these aged samples
* 500 modern reference samples (for reference panel)
* With sequence lengths of 20Mbp
* Branch scaling factor for seq-gen of 7.6e-08
* Simulate a split between aged and reference samples (OPTIONAL)
* Such that the split took place 200 generations ago (OPTIONAL)

`./run_scripts.sh 25 100 500 20000000 7.6e-08 split 200`

Look at `run_scripts.sh` for more configuration details.

This work received support from Luis Aguilar, Alejandro de León, Carlos S. Flores, and Jair García of the Laboratorio Nacional de Visualización Científica Avanzada.
