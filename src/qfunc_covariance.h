#ifndef QFUNC_COVARIANCE_H
#define QFUNC_COVARIANCE_H

#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

class QFuncCovariance {
private:
  // 2 x 5 matrix of stage 1 info.
  // Row 0 contains region 0 info and row 1 contains region 1 info.
  // Columns: (k_gamma, nugget_gamma, tau_gamma, phi_gamma, sigma_ep)
  // If sigma_ep is zero, then that region is noiseless.
  const arma::mat &stage1_parms;

public:
  QFuncCovariance(const arma::mat &stage1_parms);
  bool IsNoisy(uint region) const;
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