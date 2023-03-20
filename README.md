# Bose-Bose lattice droplets

This C++ code is devoted to studies of the ground state properties of two-component Bose-Hubbard mode in 1D within the Matrix Product State (ansatz using 
ITensor v3 C++ implementation).

The code provides:
  1. MPS for the ground state of the system for given set of Hamiltonian parameters. MPS are obtained within within desired accuracy of the entropy at the central bond
  1. Expectation value of the total Hamiltonian, and its separated parts in the ground state
  2. one- and two-body correlators
  3. density-density correlations

Using:
  1. Install ITensor v3 library (C++ implementation): https://itensor.org/docs.cgi
  2. enter Bose-Hubbard model parameters in bash script run.sh
  3. execture "bash run.sh"

Repo contains:
  1. Source code
  2. Jupyter notebook with two-component Bose-Hubbard Hamiltonian and calculated quantities in "model.ipynb"
