
#ifndef _POISSONREGRESSIONCPP_H
#define _POISSONREGRESSIONCPP_H
 
#include <RcppArmadillo.h>

double Poisson_Newton(const arma::mat& xx, const arma::vec& y, arma::vec& theta, bool happy = true);

Rcpp::List PoissonRegression(const arma::mat& xx, const arma::vec& y);

#endif // _POISSONREGRESSIONCPP_H
