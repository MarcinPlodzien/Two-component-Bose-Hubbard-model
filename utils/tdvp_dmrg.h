#include "itensor/all.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace itensor;
using namespace std;


MPS imag_time_evol(Boson sites,
                   MPS state,
                   MPO Hamiltonian,
                   int NoOfSteps = 2,
                   int nosweeps = 2,
                   Real dt_bysweep = 0.01,
                   int MaxBondDim = 200);

MPS dmrg_sequence(Boson sites,
                  MPS state,
                  vector<MPO> H_terms, // {H_total, H_hop_a, H_hop_b, H_aa, H_bb, H_ab}
                  int MaxBondDim,
                  std::string densities_entropies,
                  std::string convergence_params,
                  std::string sites_file,
                  std::string mps_file);
