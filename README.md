# Two-component Bose-Hubbard model in 1D

This C++ code uses tensor networks Matrix Product State ansatz to study ground state properties of a two-component Bose-Hubbard model in 1D (DMRG/imaginary time TDVP using ITensor v3 C++ implementation). 

In particular, you can study exotic phases such as Bose-Bose lattice quantum droplets.

The code provides (details in Jupyten notebook "model.ipynb")

  1. MPS for the ground state of the system for given set of Hamiltonian parameters. MPS are obtained within within desired accuracy of the entropy at the central bond
  1. The expectation value of the total Hamiltonian, and its separated parts in the ground state.
  2. Atom densiteis and one-body correlation matrix.
  3. Two-body correlation matrix for atomic pairs $a-b$.
  4. Density-density correlations for atoms $a$ and $b$.

Using:

  1. Install ITensor v3 library (C++ implementation): https://itensor.org/docs.cgi
  2. Install TDVP library for ITensor: https://github.com/ITensor/TDVP
  4. Enter Bose-Hubbard model parameters in bash script run.sh (number of bosons, number of sites, interactions strenght, hopping, etc.).
  5. Execture "bash run.sh".

Repo contains:

  1. Source code.
  2. Jupyter notebook with two-component Bose-Hubbard Hamiltonian and calculated quantities in "model.ipynb".
