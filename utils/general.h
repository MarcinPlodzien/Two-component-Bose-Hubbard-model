#include "itensor/all.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace itensor;
using namespace std;


tuple<MPO, MPO, MPO, MPO, MPO, MPO, MPO> get_H(SiteSet& sites,
                                               double t,
                                               double U,
                                               double U_ab,
                                               int Na,
                                               int Nb);

MPS initial_state(Boson sites, int Na, int Nb);

string parameters_to_filename(int Na,
                              int Nb,
                              int L,
                              int maxOccupation,
                              double t,
                              double U,
                              double r,
                              int MaxBondDim);

MPS apply_operator(Boson sites, string OpName, MPS ket, int site);

double density_a(Boson sites, MPS state, int site);
double density_b(Boson sites, MPS state, int site);

double entanglement_entropy(MPS state, int lattice_site);

tuple<vector<double>, vector<double>> entropies(MPS state);

tuple<vector<double>, vector<double>> particle_densities(Boson sites, MPS state);

void prepare_file(vector<string> column_names, std::string path);

void collect_densities_entropies(Boson sites,
                                 MPS state,
                                 std::string densities_entropies);

void collect_convergence_parameters(Boson sites,
                                    MPS state,
                                    vector<MPO> H_terms, // {H_total, H_hop_a, H_hop_b, H_aa, H_bb, H_ab}
                                    std::string convergence_params);

string get_string_filename(double g);
