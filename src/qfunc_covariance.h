#ifndef QFUNC_COVARIANCE_H
#define QFUNC_COVARIANCE_H

#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

class QFuncCovariance {
private:
  // 2 x 4 matrix of stage 1 parameters.
  // Row 0 contains region 1 parameters and row 1 contains region 2 parameters.
  // Parameter list: (k_gamma, nugget_gamma, tau_gamma, phi_gamma)
  arma::mat stage1_parms;

public:
  QFuncCovariance(arma::mat stage1_parms);
  arma::mat operator()(const arma::vec &covparms,
                       const arma::mat &locsub) const;
};

class DQFuncCovariance {
private:
  arma::mat stage1_parms;

public:
  DQFuncCovariance(arma::mat stage1_parms);
  arma::cube operator()(const arma::vec &covparms,
                        const arma::mat &locsub) const;
};

#endif // QFUNC_COVARIANCE_H