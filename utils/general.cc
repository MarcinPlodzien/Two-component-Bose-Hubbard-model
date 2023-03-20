#include "itensor/all.h"
#include <string>
#include <iostream>
#include <fstream>

using namespace itensor;
using namespace std;


tuple<MPO, MPO, MPO, MPO, MPO, MPO, MPO> get_H(SiteSet& sites,
                                               double t,
                                               double U,
                                               double U_ab,
                                               int Na,
                                               int Nb){
  int M_sites = length(sites);
  int L = M_sites/2;

  auto total = AutoMPO(sites);
  // kinetic terms
  auto hopping_a = AutoMPO(sites);
  auto hopping_b = AutoMPO(sites);
  for (int j_site = 1; j_site < L; j_site++){
    int j_site_a = 2*j_site - 1;
    int j_site_b = 2*j_site;
    hopping_a += -t, "Adag", j_site_a, "A", j_site_a + 2;
    hopping_a += -t, "Adag", j_site_a + 2, "A", j_site_a;
    hopping_b += -t, "Adag", j_site_b, "A", j_site_b + 2;
    hopping_b += -t, "Adag", j_site_b + 2, "A", j_site_b;

    total += -t, "Adag", j_site_a, "A", j_site_a + 2;
    total += -t, "Adag", j_site_a + 2, "A", j_site_a;
    total += -t, "Adag", j_site_b, "A", j_site_b + 2;
    total += -t, "Adag", j_site_b + 2, "A", j_site_b;
  }
  auto H_hop_a = toMPO(hopping_a);
  auto H_hop_b = toMPO(hopping_b);

  // intra-component interactions
  auto interaction_aa = AutoMPO(sites);
  auto interaction_bb = AutoMPO(sites);
  for (int j_site = 1; j_site <= L; j_site++){
    int j_site_a = 2*j_site - 1;
    int j_site_b = 2*j_site;
    interaction_aa +=  0.5*U, "N", j_site_a, "N", j_site_a;
    interaction_aa += -0.5*U, "N", j_site_a;
    interaction_bb +=  0.5*U, "N", j_site_b, "N", j_site_b;
    interaction_bb += -0.5*U, "N", j_site_b;
    
    total +=  0.5*U, "N", j_site_a, "N", j_site_a;
    total += -0.5*U, "N", j_site_a;
    total +=  0.5*U, "N", j_site_b, "N", j_site_b;
    total += -0.5*U, "N", j_site_b;
  }
  auto H_aa = toMPO(interaction_aa);
  auto H_bb = toMPO(interaction_bb);

  // inter-component interactions
  auto interaction_ab = AutoMPO(sites);
  for (int j_site = 1; j_site <= L; j_site++){
    int j_site_a = 2*j_site - 1;
    int j_site_b = 2*j_site;
    interaction_ab += U_ab, "N", j_site_a, "N", j_site_b;

    total += U_ab, "N", j_site_a, "N", j_site_b;
  }
  auto H_ab = toMPO(interaction_ab);

  // Hard wall at the edges choose value different than 0 if you want to
  // have additional potential at the first and last sites in the system.
  auto edge_potential = AutoMPO(sites);
  double wall = 0.;
  edge_potential += wall, "N", 1;
  edge_potential += wall, "N", M_sites-1;
  edge_potential += wall, "N", 2;
  edge_potential += wall, "N", M_sites;

  total += wall, "N", 1;
  total += wall, "N", M_sites-1;
  total += wall, "N", 2;
  total += wall, "N", M_sites;
  auto H_edge = toMPO(edge_potential);

  // to punish any numerical leak of particles between odd and even sites
  // we add terms fine*(Na - sum_ja op(n_ja))^2 + 1000*(Nb - sum_jb op(n_jb))^2 
  // so we need terms Nx^2, -2*Nx sum_jx n_jx, (sum_jx n_jx)^2, with x in [a, b]
  double fine = 0.1;
  // -2*sum_x Nx sum_jx n_jx term
  total += fine*Na*Na, "Id", L;
  total += fine*Nb*Nb, "Id", L;
  // sum_x sum_jx n_jx)^2 term
  for (int j_site = 1; j_site <= L; j_site++){
    int j_site_a = 2*j_site - 1;
    int j_site_b = 2*j_site;
    total += fine*(-2*Na), "N", j_site_a;
    total += fine*(-2*Nb), "N", j_site_b;
  }
  // sum_x sum_jx n_jx)^2 term
  for (int j_site_1 = 1; j_site_1 <= L; j_site_1++){
    int j_site_a_1 = 2*j_site_1 - 1;
    int j_site_b_1 = 2*j_site_1;
    for (int j_site_2 = 1; j_site_2 <= L; j_site_2++){
      int j_site_a_2 = 2*j_site_2 - 1;
      int j_site_b_2 = 2*j_site_2;      
      total += fine, "N", j_site_a_1, "N", j_site_a_2;
      total += fine, "N", j_site_b_1, "N", j_site_b_2;
    }
  }

  auto H_total = toMPO(total);

  return make_tuple(H_total, H_hop_a, H_hop_b, H_aa, H_bb, H_ab, H_edge);
}


MPS initial_state(Boson sites, int Na, int Nb){
  auto state = InitState(sites);
  int M_sites = length(sites);

  for(int i = 1; i <= M_sites; ++i)
  { state.set(i,"0"); }

  // a-atoms occupy odd [1,3,5,...] sites in the code while b-atoms occupy the even ones [2,4,6,...]
  int Qa = int(M_sites/2.) - Na;
  int Qb = int(M_sites/2.) - Nb;
  if(Qa % 2 == 0){Qa += 1;} // must be odd
  if(Qb % 2 != 0){Qb += 1;} // must be even

  for(int s = Qa; s < Qa + 2*Na;  s += 2){
    state.set(s, "1");
  }
  for(int s = Qb; s < Qb + 2*Nb;  s += 2){
    state.set(s, "1");
  }
  auto psi = MPS(state);
  return psi;
}


string get_string_filename(double g){
	char g_str_format[40];
  sprintf(g_str_format,"%05.3f",g);
	return g_str_format;
}


string parameters_to_filename(int Na,
                              int Nb,
                              int L,
                              int maxOccupation,
                              double t,
                              double U,
                              double r,
                              int MaxBondDim){

  string str_Natoms = "Na." + str(Na) + "_Nb." + str(Nb);
  string str_L = "_L." + str(L);
  string str_maxOccupation = "_MaxOcc." + str(maxOccupation);
  string str_MaxBondDim = "_MaxBondDim." + str(MaxBondDim);

  string str_t = "_t." + get_string_filename(t);
  string str_U = "_U." + get_string_filename(U);
  string str_r = "_r." + get_string_filename(r);

  return str_Natoms + str_L + str_t + str_U + str_r + str_maxOccupation + str_MaxBondDim + ".txt";
}


MPS apply_operator(Boson sites, string OpName, MPS ket, int site){
  ket.position(site);
  auto newket = noPrime(ket(site) * op(sites, OpName, site));
  ket.set(site, newket);
  return ket;
}


double density_a(Boson sites, MPS state, int site){
  auto M = length(state);
  int L = M/2;
  int site_a = 2*site - 1;
  auto ket = apply_operator(sites, "N", state, site_a);
  auto density_boson_a = innerC(state, ket);
  return density_boson_a.real();
}


double density_b(Boson sites, MPS state, int site){
  auto M = length(state);
  int L = M/2;
  int site_b = 2*site;
  auto ket = apply_operator(sites, "N", state, site_b);
  auto density_boson_b = innerC(state, ket);
  return density_boson_b.real();
}


double entanglement_entropy(MPS state, int lattice_site){
  auto psi = state;
  int bond = lattice_site;
  psi.position(bond);
  auto l = leftLinkIndex(psi, bond);
  auto s = siteIndex(psi, bond);
  auto [U,S,V] = svd(psi(bond), {l, s});
  auto u = commonIndex(U, S);
  double SvN = 0.;
  for(auto n : range1(dim(u)))
  {
    auto Sn = elt(S, n, n);
    auto p = sqr(Sn);
    if(p > 1E-12){ SvN += -p*log(p); }
  }
  // std::cout << "\n\nS at bond: " << bond << " = " << SvN;
  auto entanglement_entropy = SvN;
  return entanglement_entropy;
}


tuple<vector<double>, vector<double>> entropies(MPS state){
  auto M = length(state);
  int L = M/2;
  vector<double> entropies_a(L);
  vector<double> entropies_b(L);
  for(int site = 1; site <= L; site++){
    int site_a = 2*site - 1;
    int site_b = 2*site;
    if((site == 1) || (site == L)){
      entropies_a[site-1] = 0.;
      entropies_b[site-1] = 0.;
    }else{
      entropies_a[site-1] = entanglement_entropy(state, site_a);
      entropies_b[site-1] = entanglement_entropy(state, site_b);
    }
  }
  return make_tuple(entropies_a, entropies_b);
}


tuple<vector<double>, vector<double>> particle_densities(Boson sites, MPS state){
  auto M = length(state);
  int L = M/2;
  vector<double> densities_a(L);
  vector<double> densities_b(L);
  for(int site = 1; site <= L; site++){
    densities_a[site-1] = density_a(sites, state, site);
    densities_b[site-1] = density_b(sites, state, site);
  }
  return make_tuple(densities_a, densities_b);
}


void prepare_file(vector<std::string> column_names, std::string path){
  // open file and fill the first row with names of variables (column_names)
  std::ofstream file (path);
  int n_variables = column_names.size();
  for(int i = 0; i < n_variables; i++){
    file << column_names[i] << " ";
  }file << "\n";
  file.close();
}


void collect_densities_entropies(Boson sites,
                                 MPS state,
                                 std::string densities_entropies){
  vector<double> densities_a;
  vector<double> densities_b;
  tie(densities_a, densities_b) = particle_densities(sites, state);
  vector<double> entropies_a;
  vector<double> entropies_b;
  tie(entropies_a, entropies_b) = entropies(state);

  // print total densities
  double N_a = 0;
  double N_b = 0;
  for(int i = 0; i <= int(length(state)/2)-1; i++){
    N_a += densities_a[i];
    N_b += densities_b[i];
  }
  printfln("Na = %.5f", N_a );
  printfln("Nb = %.5f", N_b );

  // write densities and entropies to file
  std::ofstream dens_entrs (densities_entropies, std::ios::app);
  for(int i = 0; i <= int(length(state)/2)-1; i++){
    int site_number = i+1;
    dens_entrs << site_number << " ";
    dens_entrs << densities_a[i] << " " << densities_b[i] << " ";
    dens_entrs << entropies_a[i] << " " << entropies_b[i] << "\n";
  }
  dens_entrs << "\n"; // separate results obtained in consecutive steps
  dens_entrs.close();
}


void collect_convergence_parameters(Boson sites,
                                    MPS state,
                                    vector<MPO> H_terms, // {H_total, H_hop_a, H_hop_b, H_aa, H_bb, H_ab}
                                    std::string convergence_params
                                    ){
  int L = int(length(state)/2);
  int central_a = L;
  int central_b = L;
  if(L % 2 == 0){ central_a = L - 1; }
  else{ central_b = L + 1; }

  std::ofstream conv_params (convergence_params, std::ios::app);

  // entropy at central bond
  double central_entropy_a = entanglement_entropy(state, central_a);
  double central_entropy_b = entanglement_entropy(state, central_b);
  conv_params << central_entropy_a << " " << central_entropy_b << " ";

  //calculate energies
  double E_total = innerC(state, H_terms[0], state).real();
  double E_hop_a = innerC(state, H_terms[1], state).real();
  double E_hop_b = innerC(state, H_terms[2], state).real();
  double E_aa = innerC(state, H_terms[3], state).real();
  double E_bb = innerC(state, H_terms[4], state).real();
  double E_ab = innerC(state, H_terms[5], state).real();
  conv_params << E_total << " ";
  conv_params << E_hop_a << " ";
  conv_params << E_hop_b << " ";
  conv_params << E_aa << " ";
  conv_params << E_bb << " ";
  conv_params << E_ab << " ";
  conv_params << E_hop_a + E_hop_b + E_aa + E_bb + E_ab<<" ";
  conv_params << "\n";
  conv_params.close();
}
