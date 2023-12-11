#include <RcppArmadillo.h>   
#include <iostream>

using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Mat<double> calc_precision_cpp(arma::mat kron_sum_theta, arma::mat theta,
                      arma::mat theta_t, arma::mat sigma2_vec,
                      double ni, arma::vec times) {
  
  //Rcpp::Rcout << t_jminus1 << std::endl; 
  

  // some preliminary calculatuions
  arma::mat V;
  V = arma::inv(kron_sum_theta);
  V = V * sigma2_vec; 
  V = arma::reshape(V, theta.n_rows, theta.n_rows); // stationary variance
  arma::mat Vinv;
  Vinv = arma::inv(V); // inverse of stationary variance
  arma::mat diagonal;
  arma::mat offdiagonal;
  arma::mat tempterm1;
  arma::mat tempterm2;
  arma::mat tempterm3;
  arma::mat tempterm4;
  arma::mat tempterm5;
  arma::mat Omega(theta.n_rows*ni, theta.n_rows*ni, fill::zeros);
  
  double t_j;
  double t_jplus1;
  double t_jminus1;
  
  for(int j = 0; j < ni; j++){
    t_j = times[j];
    if (j == 0){
      // Omega_11
      t_jplus1 = times[j+1];
      tempterm1 = arma::expmat(-theta_t * (t_jplus1 - t_j));
      tempterm2 = arma::expmat(-theta * (t_jplus1 - t_j));
      diagonal = arma::inv(V - V * tempterm1 * Vinv * tempterm2 * V);
    } else if (j == (ni-1)){
      // Omega_nn:
      t_jminus1 = times[j-1];
      tempterm1 = arma::expmat(-theta * (t_j - t_jminus1));
      tempterm1 = tempterm1 * V;
      tempterm2 = arma::expmat(-theta_t * (t_j - t_jminus1));
      tempterm2 = tempterm2 * Vinv;
      
      // Omega_nn
      diagonal = Vinv + Vinv * tempterm1 * arma::inv(V - V * tempterm2 * tempterm1) * V * tempterm2;
    
      // Omega_n-1,n
      offdiagonal = -arma::inv(V - V * tempterm2 * tempterm1) * V * tempterm2;
    
    } else {
      // Omega_j,j-1, Omega_j,j, Omega_j,j+1 where Omega_j,j-1 = t(Omega_j,j+1)
      t_jminus1 = times[j-1];
      t_jplus1 = times[j+1];
      tempterm1 = arma::expmat(-theta_t * (t_jplus1 - t_j));
      tempterm2 = tempterm1 * Vinv;
      tempterm3 = arma::expmat(-theta * (t_j - t_jminus1));
      tempterm4 = tempterm3 * V;
      tempterm5 = arma::trans(tempterm3) * Vinv;

      // Omega_j,j
      diagonal = Vinv + Vinv * tempterm4 * arma::inv(V - V * tempterm5 * tempterm4) *
        V * tempterm5 + arma::inv(V - V * tempterm2 * arma::trans(tempterm1) * V) *
        V * tempterm2 * arma::trans(tempterm1);
        
      // Omega_j-1,j
      offdiagonal = -arma::inv(V - V * tempterm5 * tempterm4) * V * tempterm5;
    }
    
    // Update diagonals: Omega[j,j]
    Omega.submat(j*theta.n_rows, j*theta.n_rows,
                 j*theta.n_rows + theta.n_rows-1, j*theta.n_rows + theta.n_rows-1) = diagonal;
    if (j != 0){
      // Update off-diagonals: Omega[j-1,j] and Omega[j,j-1]
      Omega.submat((j-1)*theta.n_rows, j*theta.n_rows,
                   (j-1)*theta.n_rows + theta.n_rows-1, j*theta.n_rows + theta.n_rows-1) = offdiagonal;
      Omega.submat(j*theta.n_rows, (j-1)*theta.n_rows,
                   j*theta.n_rows + theta.n_rows-1, (j-1)*theta.n_rows + theta.n_rows-1) = arma::trans(offdiagonal);
    }
  }
  return Omega;
}

