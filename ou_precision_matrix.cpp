#include <RcppArmadillo.h>   
#include <iostream>

using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Mat<double> calc_precision_cpp(arma::mat kron_sum_theta, arma::mat theta,
                      arma::mat theta_t, arma::mat sigma2_vec,
                      double ni, arma::vec times) {
  
  double t_jminus1;
  
  //Rcpp::Rcout << t_jminus1 << std::endl; 
  

  // some preliminary calculatuions
  arma::mat var1;
  var1 = arma::inv(kron_sum_theta);
  var1 = var1 * sigma2_vec;
  var1 = arma::reshape(var1, theta.n_rows, theta.n_rows);
  arma::mat temp;
  arma::mat offdiagonal;
  arma::mat tempterm1;
  arma::mat tempterm2a;
  arma::mat tempterm2b;
  arma::mat tempterm2c;
  arma::mat tempterm2;
  arma::mat Omega(theta.n_rows*ni, theta.n_rows*ni, fill::zeros);
  
  for(int j = 0; j < ni; j++){
    double t_j = times[j];
    if (j == 0){
      //cout << "if " << t_j << endl;
      double t_jplus1 = times[j+1];
      temp = arma::inv(var1 - arma::expmat(-theta * (t_jplus1 - t_j)) * var1 * arma::expmat(-theta_t * (t_jplus1 - t_j)));
      tempterm1 = arma::inv(var1 - arma::expmat(-theta * (t_jplus1 - t_j))
                     * var1 * arma::expmat(-theta_t*(t_jplus1-t_j)));
      tempterm2 = arma::expmat(-theta * (t_jplus1 - t_j));
      offdiagonal = -(tempterm1 * tempterm2);
        
    } else if (j == (ni-1)){
      //cout << "else if " << t_j << endl;
      t_jminus1 = times[j-1];
      tempterm2a = arma::expmat(-theta_t * (t_j - t_jminus1));
      tempterm1 = arma::inv(var1);
      tempterm2b = arma::inv(var1 - arma::expmat(-theta * (t_j - t_jminus1)) * var1 * arma::expmat(-theta_t*(t_j - t_jminus1)));
      tempterm2c = arma::expmat(-theta * (t_j-t_jminus1));
      temp = tempterm2a * (tempterm2b * tempterm2c) + tempterm1;
    } else {
      //cout << "else " << t_j << endl;
      t_jminus1 = times[j-1];
      double t_jplus1 = times[j+1];
      tempterm2a = arma::expmat(-theta_t * (t_j - t_jminus1));
      tempterm1 = arma::inv(var1 - arma::expmat(-theta * (t_jplus1 - t_j)) * var1 * arma::expmat(-theta_t * (t_jplus1 - t_j)));
      tempterm2b = arma::inv(var1 - arma::expmat(-theta * (t_j - t_jminus1)) * var1 * arma::expmat(-theta_t*(t_j - t_jminus1)));
      tempterm2c = arma::expmat(-theta * (t_j-t_jminus1));
      tempterm2 = tempterm2a * (tempterm2b * tempterm2c);
      temp = tempterm1 + tempterm2;
      
      tempterm1 = arma::inv(var1 - arma::expmat(-theta * (t_jplus1 - t_j))
                              * var1 * arma::expmat(-theta_t*(t_jplus1-t_j)));
      tempterm2 = arma::expmat(-theta * (t_jplus1 - t_j));
      offdiagonal = -(tempterm1 * tempterm2);
    }
    //update precision matrix Omega
    Omega.submat(j*theta.n_rows, j*theta.n_rows,
                 j*theta.n_rows + theta.n_rows-1, j*theta.n_rows + theta.n_rows-1) = temp;
    if (j != (ni-1)){
      Omega.submat(j*theta.n_rows, (j+1)*theta.n_rows,
                   j*theta.n_rows + theta.n_rows-1, (j+1)*theta.n_rows + theta.n_rows-1) = offdiagonal;
      Omega.submat((j+1)*theta.n_rows, j*theta.n_rows,
                   (j+1)*theta.n_rows + theta.n_rows-1, j*theta.n_rows + theta.n_rows-1) = arma::trans(offdiagonal);
    }

  }
  
  //cout << theta.n_rows << endl;
  
  return Omega;
  
}






