# Bose-Bose lattice droplets

This C++ code uses Matrix Product State ansatz to study ground state properties of a  Bose-Hubbard model in 1D (DMRG/imaginary time TDVP using ITensor v3 C++ implementation).

The code provides (details in Jupyten notebook "model.ipynb")

  1. MPS for the ground state of the system for given set of Hamiltonian parameters. MPS are obtained within within desired accuracy of the entropy at the central bond
  1. Expectation value of the total Hamiltonian, and its separated parts in the ground state
  2. one- and two-body correlators
  3. density-density correlations

Using:

  1. Install ITensor v3 library (C++ implementation): https://itensor.org/docs.cgi
  2. Enter Bose-Hubbard model parameters in bash script run.sh (number of bosons, number of sites, interactions strenght, hopping, etc.)
  3. Execture "bash run.sh"

Repo contains:

  1. Source code
  2. Jupyter notebook with two-component Bose-Hubbard Hamiltonian and calculated quantities in "model.ipynb"
