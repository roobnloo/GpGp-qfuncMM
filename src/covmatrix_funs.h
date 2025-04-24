#ifndef ARMACOVMATRIX_FUNS_H
#define ARMACOVMATRIX_FUNS_H

#include <RcppArmadillo.h>

#include "qfunc_covariance.h"

using namespace Rcpp;
using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]

void get_covfun(
    std::string covfun_name_string,
    std::function<arma::mat(const arma::vec &, const arma::mat &)> p_covfun[1],
    std::function<arma::cube(const arma::vec &, const arma::mat &)>
        p_d_covfun[1],
    const mat &additional_info) {

  if (covfun_name_string.compare("qfuncmm") == 0) {
    p_covfun[0] = QFuncCovariance(additional_info);
    p_d_covfun[0] = DQFuncCovariance(additional_info);
  } else {
    Rcpp::Rcout << "Unrecognized Covariance Function Name \n";
  }
}

#endif
