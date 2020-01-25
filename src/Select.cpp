#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <vector>
#include <map>

#include "PoissonRegression.h"
#include "indptSTM.h"


// [[Rcpp::export]]
Rcpp::List Select3(const arma::mat& dat1, const arma::mat& dat2, const arma::mat& dat3, 
                       const int n, const arma::mat& all_models){
  
  // Organize data 
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
  
  // Moedl selction ======================================================
  int model_count = all_models.n_cols, nn = y.n_rows;
  arma::mat selected_beta(K, K);
  arma::vec selected_beta0(K);
  
  for(int k = 0; k != K; ++k){
    
    // Fit full model, get initial estimate
    auto full = PoissonRegression(x, y.col(k));
    double lik = full["likelihood"];
    arma::vec theta = full["parameters"];
    double criterion = 2*lik + log(nn)*x.n_cols;
    
    // model selection (exhaustive search)
    int opt_ind = 0; 
    for(int md = 1; md != model_count; ++md){
      arma::uvec v_u = arma::find(arma::conv_to<arma::vec>::from(all_models.col(md)) > 0);
      arma::mat x_sub = x.cols(v_u); 
      arma::vec theta_sub = theta.elem(v_u);
      lik = Poisson_Newton(x_sub, y.col(k), theta_sub);
      double criterion_new = 2*lik + log(nn)*x_sub.n_cols;
      if(criterion_new < criterion){
        opt_ind = md;
        criterion = criterion_new;
      }
    }
  
    // write result 
    arma::uvec v_u = arma::find(arma::conv_to<arma::vec>::from(all_models.col(opt_ind)) > 0);
    arma::mat x_sub = x.cols(v_u); 
    arma::vec theta_sub = theta.elem(v_u);
    Poisson_Newton(x_sub, y.col(k), theta_sub);
    arma::vec theta_k(K + 1, arma::fill::zeros);
    theta_k.elem(v_u) = theta_sub;
    selected_beta0(k) = theta_k(0);
    selected_beta.col(k) = theta_k.tail(K);
  }
  
  return Rcpp::List::create(Rcpp::Named("intercept") = selected_beta0,
                            Rcpp::Named("main_effects") = selected_beta);
  
}