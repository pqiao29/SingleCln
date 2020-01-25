#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "grid_ring.h"
#include "PoissonRegression.h"

// [[Rcpp::plugins(cpp11)]]

Rcpp::List data_indpt(const arma::mat& dat, const int n){
  
  std::map<int, std::vector<int> > nbr = grid_ring(n, 1);
  
  int t_size = dat.col(0).max() + 1, K = dat.col(1).max(), n_lattice = n*n;
  
  // count
  arma::cube obs_cnt(K, n_lattice, t_size, arma::fill::zeros);
  for(int i = 0; i != dat.n_rows; ++i){
    ++obs_cnt((dat(i, 1) - 1), (dat(i, 2) - 1), dat(i, 0));
  }
  
  arma::mat response((n_lattice * (t_size - 1)), K); //uninitialized
  arma::mat covariate((n_lattice * (t_size - 1)), K + 1); covariate.col(0).ones();
  int ind_cov = 0, ind_res = 0; 
  for(int t = 1; t != t_size; ++t){
    // response
    response.rows(ind_res, (ind_res + n_lattice - 1)) = obs_cnt.slice(t).t();
    ind_res += n_lattice;
    // covariate
    for(int i = 0; i != n_lattice; ++i){
      std::vector<int> tmp_nb = nbr[i + 1];
      double tmp; 
      for(int k = 0; k != K; ++k){
        tmp = 0; 
        std::vector<int>::const_iterator iter = tmp_nb.cbegin();
        while(iter != tmp_nb.cend()){
          tmp += log(obs_cnt(k, (*iter++ - 1), t - 1) + 1);
        }
        tmp /= tmp_nb.size();
        covariate(ind_cov, k + 1) = tmp; 
      }
      ++ind_cov;
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("response") = response,
                            Rcpp::Named("covariates") = covariate, 
                            Rcpp::Named("T") = t_size, 
                            Rcpp::Named("K") = K);
  
} 

// [[Rcpp::export]]
Rcpp::List indptSTM_cpp3(const arma::mat& dat1, const arma::mat& dat2, const arma::mat& dat3, const int n){
  
  /*
   * Requirement of data columns: Timepoint (0, 1, 2, ...), group (1, 2, 3, ...), tile (1, 2, ... n*n)
   */
  
  auto in_data1 = data_indpt(dat1, n);
  auto in_data2 = data_indpt(dat2, n);
  auto in_data3 = data_indpt(dat3, n);
  
  int K = in_data1["K"], t_size = in_data1["T"], n_lattice = n*n, 
      d = n_lattice * (t_size - 1);
  
  arma::mat x((3*d), K + 1), y((3*d), K); 
  x.head_rows(d) = Rcpp::as<arma::mat>(in_data1["covariates"]); 
  x.rows(d, (2*d - 1)) = Rcpp::as<arma::mat>(in_data2["covariates"]);
  x.tail_rows(d) = Rcpp::as<arma::mat>(in_data3["covariates"]);
  y.head_rows(d) = Rcpp::as<arma::mat>(in_data1["response"]); 
  y.rows(d, (2*d - 1)) = Rcpp::as<arma::mat>(in_data2["response"]);
  y.tail_rows(d) = Rcpp::as<arma::mat>(in_data3["response"]);
  
  double lik = 0; 
  
  arma::rowvec beta0(K); arma::mat beta(K, K);
  int ind_b0 = 0;
  for(int k = 0; k != K; ++k){
    auto glm_res = PoissonRegression(x, y.col(k));
    arma::vec tmp_res = glm_res["parameters"];
    beta0(ind_b0++) = tmp_res(0);
    beta.col(k) = tmp_res.tail(K);
    lik += Rcpp::as<double>(glm_res["likelihood"]);
  }

  return Rcpp::List::create(Rcpp::Named("likelihood") = lik,
                            Rcpp::Named("intercept") = beta0, 
                            Rcpp::Named("main_effects") = beta);
  
}


