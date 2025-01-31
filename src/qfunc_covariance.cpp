#include "qfunc_covariance.h"

QFuncCovariance::QFuncCovariance(arma::mat stage1_parms)
    : stage1_parms(stage1_parms) {}

arma::mat QFuncCovariance::operator()(const arma::vec &covparms,
                                      const arma::mat &locsub) const {
  if (locsub.n_cols != 4) {
    stop("3 spatial and 1 temporal column expected, but %d columns found.",
         locsub.n_cols);
  }
  return arma::eye(locsub.n_rows, locsub.n_cols);
}

DQFuncCovariance::DQFuncCovariance(arma::mat stage1_parms)
    : stage1_parms(stage1_parms) {}

arma::cube DQFuncCovariance::operator()(const arma::vec &covparms,
                                        const arma::mat &locsub) const {
  stop("Not implemented");
  return arma::zeros(locsub.n_rows, locsub.n_cols, covparms.n_elem);
}