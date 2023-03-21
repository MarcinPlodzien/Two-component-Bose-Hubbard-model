#include "itensor/all.h"
#include <string>
#include <iostream>
#include <fstream>
#include "general.h"

using namespace itensor;
using namespace std;


double density_density_correlation_a(Boson sites, MPS state, int i_site, int j_site){
  auto M = length(state);
  int L = M/2;
  int i_site_a = 2*i_site - 1;
  int j_site_a = 2*j_site - 1;
  auto ket = apply_operator(sites, "N", state, j_site_a);
       ket = apply_operator(sites, "N", ket, i_site_a);
  auto N_a_i_N_a_j = innerC(state, ket);
  return N_a_i_N_a_j.real();
}


double density_density_correlation_b(Boson sites, MPS state, int i_site, int j_site){
  auto M = length(state);
  int L = M/2;
  int i_site_b = 2*i_site;
  int j_site_b = 2*j_site;
  auto ket = apply_operator(sites, "N", state, j_site_b);
       ket = apply_operator(sites, "N", ket, i_site_b);
  auto N_b_i_N_b_j = innerC(state, ket);
  return N_b_i_N_b_j.real();
}


double one_body_correlation_a(Boson sites, MPS state, int i_site, int j_site){
  auto M = length(state);
  int L = M/2;
  int i_site_a = 2*i_site - 1;
  int j_site_a = 2*j_site - 1;
  auto ket = apply_operator(sites, "A", state, j_site_a);
  auto bra = apply_operator(sites, "A", state, i_site_a);
  auto a_dag_i_a_j = innerC(bra, ket);
  return a_dag_i_a_j.real();
}


double one_body_correlation_b(Boson sites, MPS state, int i_site, int j_site){
  auto M = length(state);
  int L = M/2;
  int i_site_b = 2*i_site;
  int j_site_b = 2*j_site;
  auto ket = apply_operator(sites, "A", state, j_site_b);
  auto bra = apply_operator(sites, "A", state, i_site_b);
  auto b_dag_i_b_j = innerC(bra, ket);
  return b_dag_i_b_j.real();
}


double pair_correlation_ab(Boson sites, MPS state, int i_site, int j_site){
  auto M = length(state);
  int L = M/2;
  int i_site_a = 2*i_site - 1;
  int j_site_a = 2*j_site - 1;
  int i_site_b = 2*i_site;
  int j_site_b = 2*j_site;

  auto ket_tmp = apply_operator(sites, "Adag", state, j_site_a);
  auto ket = apply_operator(sites, "Adag", ket_tmp, j_site_b);
  auto bra_tmp = apply_operator(sites, "Adag", state, i_site_b);
  auto bra = apply_operator(sites, "Adag", bra_tmp, i_site_a);
  auto a_dag_b_dag_b_a = innerC(bra, ket);
  return a_dag_b_dag_b_a.real();
}


tuple<vector<vector<double>>,
      vector<vector<double>>,
      vector<vector<double>>,
      vector<vector<double>>,
      vector<vector<double>>> correlations(Boson sites, MPS state){
// Compute all the correlations
  auto M = length(state);
  int L = M/2;
  vector<vector<double>> one_body_correlations_a(L, vector<double>(L));
  vector<vector<double>> one_body_correlations_b(L, vector<double>(L));
  vector<vector<double>> pair_correlations_ab(L, vector<double>(L));
  vector<vector<double>> density_density_a(L, vector<double>(L));
  vector<vector<double>> density_density_b(L, vector<double>(L));

  cout << "\ncalculate correlations:";
  for(int i_site = 1; i_site <= L; i_site++){
    for(int j_site = 1; j_site <= L; j_site++){
      double adagi_aj = one_body_correlation_a(sites, state, i_site, j_site);
      double bdagi_bj = one_body_correlation_b(sites, state, i_site, j_site);
      double ai_bi_bdagj_adagj = pair_correlation_ab(sites, state, i_site, j_site);
      double nai_naj = density_density_correlation_a(sites, state, i_site, j_site);
      double nbi_nbj = density_density_correlation_b(sites, state, i_site, j_site);
      double nai = density_a(sites, state, i_site);
      double naj = density_a(sites, state, j_site);
      double nbi = density_b(sites, state, i_site);
      double nbj = density_b(sites, state, j_site);

      one_body_correlations_a[i_site-1][j_site-1] = adagi_aj;
      one_body_correlations_b[i_site-1][j_site-1] = bdagi_bj;
      pair_correlations_ab[i_site-1][j_site-1] = ai_bi_bdagj_adagj - adagi_aj*bdagi_bj;
      density_density_a[i_site-1][j_site-1] = nai_naj - nai*naj;
      density_density_b[i_site-1][j_site-1] = nbi_nbj - nbi*nbj;
    }cout << "\n" << round(i_site*100./L) << "%" << std::flush;
  }
  return make_tuple(one_body_correlations_a,
                    one_body_correlations_b,
                    pair_correlations_ab,
                    density_density_a,
                    density_density_b);
}


vector<complex<double>> get_fourier_transform(vector<vector<double>> A){

  int L = A.size();

  vector<complex<double>> A_fourier_transform(L, 0);
  double pi = 3.14159265359;
  for (int k_int = 0; k_int < L; k_int++){
      double k = k_int*2.0*pi/L;
      for(int i = 0; i < L; i++){
          for(int j = 0; j < L; j++){
              A_fourier_transform[k_int] = A_fourier_transform[k_int] + exp(Cplx_i*k*(i-j))*A[i][j]/L/L;
          }
      }
  }

  return A_fourier_transform;
}


void save_correlations(Boson sites, MPS state, 
					   string path_correlations, 
					   string path_fourier_transforms, 
					   string path_scalars){

  vector<vector<double>> one_body_correlations_a;
  vector<vector<double>> one_body_correlations_b;
  vector<vector<double>> pair_correlations_ab;
  vector<vector<double>> density_density_a;
  vector<vector<double>> density_density_b;

  tie(one_body_correlations_a,
      one_body_correlations_b,
      pair_correlations_ab,
      density_density_a,
      density_density_b) = correlations(sites, state);

  int n_rows = one_body_correlations_a.size();
  int n_cols = n_rows;
  std::ofstream file_correlations (path_correlations);
  // set labels
  file_correlations << "# row" << " " << "col" << " ";
  file_correlations << "one_body_correlations_a" << " ";
  file_correlations << "one_body_correlations_b" << " ";
  file_correlations << "pair_correlations_ab" << " ";
  file_correlations << "density_density_a" << " ";
  file_correlations << "density_density_a" << "\n";

  for(int row = 0; row < n_rows; row++){
    for(int col = 0; col < n_cols; col++){
      file_correlations << row+1 << " " << col+1 << " ";
      file_correlations << one_body_correlations_a[row][col] << " ";
      file_correlations << one_body_correlations_b[row][col] << " ";
      file_correlations << pair_correlations_ab[row][col] << " ";
      file_correlations << density_density_a[row][col] << " ";
      file_correlations << density_density_b[row][col] << "\n";
    }
    file_correlations<<"\n";
  }
  file_correlations.close();

  // Get scalars
  int L = n_rows;
  double h_a = 0;
  double h_b = 0;
  double h_ab = 0;
  double f_SF_a = 0;
  double f_SF_b = 0;
  double N_a = 0;
  double N_b = 0;
  for(int i_site = 1; i_site <= L; i_site++){
    for(int j_site = 1; j_site <= L; j_site++){
       
        h_a = h_a + one_body_correlations_a[i_site-1][j_site-1];
        h_b = h_b + one_body_correlations_b[i_site-1][j_site-1];
        h_ab = h_ab + pair_correlations_ab[i_site-1][j_site-1];
    
        f_SF_a = f_SF_a + one_body_correlations_a[i_site-1][j_site-1];
        f_SF_b = f_SF_b + one_body_correlations_b[i_site-1][j_site-1];
        if(i_site == j_site){
          N_a = N_a + one_body_correlations_a[i_site-1][j_site-1];
          N_b = N_b + one_body_correlations_b[i_site-1][j_site-1];
        }
    }
  }
  f_SF_a = f_SF_a/(L*N_a);
  f_SF_b = f_SF_b/(L*N_b);
  h_a = h_a/(L*N_a);
  h_b = h_b/(L*N_b);
  h_ab = h_ab/L/(N_a*N_b);

  cout << "N_a = "<< N_a << endl;
  cout << "N_b = "<< N_b << endl;
  //

  std::ofstream file_scalars (path_scalars);
  file_scalars << "# f_SF_a f_SF_b h_a h_b h_ab\n";
  file_scalars << f_SF_a << " " << f_SF_b << " " << h_a << " " << h_b << " " << h_ab << endl;
  file_scalars.close();

  // Get Fourier transforms
  vector<complex<double>> S_fourier;
  vector<complex<double>> R_fourier;
  vector<complex<double>> Q_fourier;

  S_fourier = get_fourier_transform(density_density_a);
  R_fourier = get_fourier_transform(one_body_correlations_a);
  Q_fourier = get_fourier_transform(pair_correlations_ab);

  std::ofstream file_fourier_transforms(path_fourier_transforms);
  file_fourier_transforms<<"# S R Q \n";
  for(int k_i = 0; k_i < L; k_i++){
    file_fourier_transforms << k_i << " " << S_fourier[k_i].real() << " " << R_fourier[k_i].real() << " " << Q_fourier[k_i].real() << endl;
  }
  file_fourier_transforms.close();
}
