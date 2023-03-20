#include "itensor/all.h"
#include "tdvp.h"
#include "basisextension.h"
#include "general.h"
#include <string>
#include <iostream>
#include <fstream>

using namespace itensor;
using namespace std;


MPS imag_time_evol(Boson sites,
                   MPS state,
                   MPO H_total,
                   int NoOfSteps,
                   int nosweeps,
                   Real dt_bysweep,
                   int MaxBondDim)
{
  Real tstep = nosweeps*dt_bysweep; // time evolved during the single round of sweeps

  auto sweepst = Sweeps(nosweeps);
  sweepst.maxdim() = MaxBondDim;
  sweepst.cutoff() = 1E-18;
  sweepst.niter() = 50;

  int cut0 = 0;
  int cut1 = 0;
  auto psi = state;

  for(int n = 1; n <= NoOfSteps ; ++n){
    printf("\n ================= \n n = ", n);
    printf("\n time form ", (n-1)*tstep, " to ", n*tstep, "\n ================= \n");
    if(maxLinkDim(psi) < MaxBondDim){ // first, expand basis if the bond dimension does not exceed MaxBondDim
      cut1 = 0;
      std::vector<Real> epsilonK = {1E-10, 1E-10};
      // Play with the numbers if needed. Here, the evolution procedure is just for 'warming up' before performing DMRG
      addBasis(psi, H_total, epsilonK, {"Cutoff", 1E-10,
                                        "Method", "DensityMatrix",
                                        "KrylovOrd", 3,
                                        "DoNormalize", true,
                                        "Quiet", true});
      }else{cut0 = 1;} // if the bond dimension is already >= MaxBondDim then change cut0 --> 1
                       // so that the basis wont be expanded again
      // TDVP sweep
      if(cut0 == 0){ // evolve when the basis was expanded
        auto energy = tdvp(psi, H_total, -dt_bysweep, sweepst, {"Truncate", false, // and do not truncate
                                                                "DoNormalize", true,
                                                                "Quiet", true,
                                                                "NumCenter", 1});
      }else{
        if(cut1 == 0){
        auto energy = tdvp(psi, H_total, -dt_bysweep, sweepst, {"Truncate", true, // truncate in such a case
                                                                "DoNormalize", true,
                                                                "Quiet", true,
                                                                "NumCenter", 1});
        cut1 = 1; // after evolution with the truncation set cut1 --> 1
      }else{ // the case where cut1 = 1
        auto energy = tdvp(psi, H_total, -dt_bysweep, sweepst, {"Truncate", false, // no truncation
                                                                "DoNormalize", true,
                                                                "Quiet", true,
                                                                "NumCenter", 1});
           }
           }

    printfln("Maximum MPS bond dimension after time evolution is %d", maxLinkDim(psi));
  }

  return psi;
}


MPS dmrg_sequence(Boson sites,
                  MPS state,
                  vector<MPO> H_terms, // {H_total, H_hop_a, H_hop_b, H_aa, H_bb, H_ab}
                  int MaxBondDim,
                  std::string densities_entropies,
                  std::string convergence_params,
                  std::string sites_file,
                  std::string mps_file){

  MPO H_total = H_terms[0];

  int dim = 64;
  vector<int> BondDim_truncation = {};
  while(dim < 0.75*MaxBondDim){
    BondDim_truncation.push_back(dim);
    dim = 2*dim;
  }
  BondDim_truncation.push_back(MaxBondDim);
  int size = BondDim_truncation.size();

  // Initial DMRG sweeps
  auto psi = state;
  for(int n = 0; n < size; n++){
    int nswep = 4;
    auto sweeps = Sweeps(nswep);
    sweeps.maxdim() = BondDim_truncation[n];
    sweeps.cutoff() = 1E-16;
    for(int m = 0; m < 10; m++){
      std::cout << "\n ======= step " << m+1 << "   truncation = " << BondDim_truncation[n] << "\n";
      auto [energy0, psi0] = dmrg(H_total, psi, sweeps, {"Quiet",true});
      psi = psi0;
    }
    // collect densities and entropies
    collect_densities_entropies(sites, psi, densities_entropies);
    // collect convergence_parameters
    collect_convergence_parameters(sites, psi, H_terms, convergence_params);
  }

  // calculate entropies at central bonds
  int L = int(length(state)/2);
  int central_a = L;
  int central_b = L;
  if(L % 2 == 0){ central_a = L - 1; }
  else{ central_b = L + 1; }
  double central_entropy_a = entanglement_entropy(psi, central_a);
  double central_entropy_b = entanglement_entropy(psi, central_b);
  double postdmrg_central_entropy_a = 0.1*central_entropy_a; // initialization
  double postdmrg_central_entropy_b = 0.1*central_entropy_b; // initialization
  double relative_diff_a = fabs(central_entropy_a - postdmrg_central_entropy_a)/central_entropy_a;
  double relative_diff_b = fabs(central_entropy_b - postdmrg_central_entropy_b)/central_entropy_b;

  while((relative_diff_a > 1e-3) || (relative_diff_b > 1e-3)){
    int nswep = 4;
    auto sweeps = Sweeps(nswep);
    sweeps.maxdim() = MaxBondDim;
    sweeps.cutoff() = 1E-16;
    for(int m = 0; m < 10; m++){
      std::cout << "\n ======= step " << m+1 << "   max truncation = " << MaxBondDim << "\n";
      auto [energy0, psi0] = dmrg(H_total, psi, sweeps, {"Quiet",true});
      psi = psi0;
    }
    // collect densities and entropies
    collect_densities_entropies(sites, psi, densities_entropies);
    // collect convergence_parameters
    collect_convergence_parameters(sites, psi, H_terms, convergence_params);

    // save sites and mps`
    writeToFile(sites_file, sites);
    writeToFile(mps_file, psi);

    // recalculate entropy difference
    central_entropy_a = postdmrg_central_entropy_a;
    central_entropy_b = postdmrg_central_entropy_b;

    postdmrg_central_entropy_a = entanglement_entropy(psi, central_a);
    postdmrg_central_entropy_b = entanglement_entropy(psi, central_b);

    relative_diff_a = fabs(central_entropy_a - postdmrg_central_entropy_a)/central_entropy_a;
    relative_diff_b = fabs(central_entropy_b - postdmrg_central_entropy_b)/central_entropy_b;
  }

  printfln("Ground State Found!");

return psi;
}
