#include "qfunc_covariance.h"

QFuncCovariance::QFuncCovariance(const arma::mat &stage1_parms)
    : stage1_parms(stage1_parms) {}

bool QFuncCovariance::IsNoisy(uint region) const {
  return stage1_parms(region, 4) > datum::eps;
}

// covparms is a vector with covariance parameters
// of the form (rho, k_eta0, k_eta1, tau_eta, nugget_eta)
arma::mat QFuncCovariance::operator()(const arma::vec &covparms,
                                      const arma::mat &locs) const {
  if (locs.n_cols != 5) {
    stop("3 spatial, 1 temporal column, and 1 region column expected, but %d "
         "columns found.",
         locs.n_cols);
  }
  if (covparms.n_elem != 5) {
    stop("5 covariance parameters expected, but %d found.", covparms.n_cols);
  }
  double rho = covparms(0);
  double k_eta[2] = {covparms(1), covparms(2)};
  double tau_eta = covparms(3);
  double nugget_eta = covparms(4);

  uint n = locs.n_rows;
  span voxel = span(0, 2);
  span time = span(3);
  span region = span(4);
  mat covmat(n, n);
  for (uint i0 = 0; i0 < n; i0++) {
    for (uint i1 = 0; i1 <= i0; i1++) {
      uint r0 = round(as_scalar(locs(i0, region)));
      uint r1 = round(as_scalar(locs(i1, region)));
      double time_diff =
          std::abs(as_scalar(locs(i0, time)) - as_scalar(locs(i1, time)));

      // Start by computing the time RBF (A matrix)
      double cov = exp(-pow(tau_eta * time_diff, 2) / 2);
      if (time_diff <= datum::eps) {
        cov += nugget_eta;
      }
      cov *= sqrt(k_eta[r0] * k_eta[r1]);

      rowvec v0 = locs(i0, voxel);
      rowvec v1 = locs(i1, voxel);
      if (r0 != r1) {
        // r0 and r1 belong to different regions
        cov *= rho;
      } else {
        // r0 and r1 belong to the same region
        const uint &r = r0;
        // compute the intra-regional time RBF (B matrix)
        double b = stage1_parms(r, 0) *
                   exp(-pow(stage1_parms(r, 2) * time_diff, 2) / 2);
        if (time_diff <= datum::eps) {
          b += stage1_parms(r, 1);
        }

        // compute the intra-regional spatial matern 5/2 (C matrix)
        double dist = norm(v0 - v1);
        double phi_dist = stage1_parms(r, 3) * dist;
        double c = exp(-phi_dist * sqrt(5)) *
                   (1 + phi_dist * sqrt(5) + pow(phi_dist, 2) * 5.0 / 3.0);
        cov += b * c;
        if (IsNoisy(r) && time_diff <= datum::eps && dist <= datum::eps) {
          cov += 1; // add idiosyncratic voxel noise if applicable
        }
      }
      if (IsNoisy(r0)) {
        cov *= stage1_parms(r0, 4);
      }
      if (IsNoisy(r1)) {
        cov *= stage1_parms(r1, 4);
      }
      covmat(i0, i1) = cov;
      covmat(i1, i0) = cov;
    }
  }
  return covmat;
}

DQFuncCovariance::DQFuncCovariance(const arma::mat &stage1_parms)
    : stage1_parms(stage1_parms) {}

bool DQFuncCovariance::IsNoisy(uint region) const {
  return stage1_parms(region, 4) > datum::eps;
}

// covparms is a vector with covariance parameters
// of the form (rho, k_eta0, k_eta1, tau_eta, nugget_eta)
// These correspond to slices in the cube that is returned.
arma::cube DQFuncCovariance::operator()(const arma::vec &covparms,
                                        const arma::mat &locs) const {
  if (locs.n_cols != 5) {
    stop("3 spatial, 1 temporal column, and 1 region column expected, but %d "
         "columns found.",
         locs.n_cols);
  }
  if (covparms.n_elem != 5) {
    stop("5 covariance parameters expected, but %d found.", covparms.n_cols);
  }
  double rho = covparms(0);
  double k_eta[2] = {covparms(1), covparms(2)};
  double tau_eta = covparms(3);
  double nugget_eta = covparms(4);

  uint n = locs.n_rows;
  span time = span(3);
  span region = span(4);
  uword rho_slice = 0, k_eta0_slice = 1, k_eta1_slice = 2, tau_eta_slice = 3,
        nugget_eta_slice = 4;
  cube dcovmat(n, n, covparms.n_elem, fill::zeros);
  for (uint i0 = 0; i0 < n; i0++) {
    for (uint i1 = 0; i1 <= i0; i1++) {
      uint r0 = round(as_scalar(locs(i0, region)));
      uint r1 = round(as_scalar(locs(i1, region)));
      double time_diff =
          std::abs(as_scalar(locs(i0, time)) - as_scalar(locs(i1, time)));
      double sigep0 = stage1_parms(r0, 4), sigep1 = stage1_parms(r1, 4);
      if (!IsNoisy(r0)) {
        sigep0 = 1;
      }
      if (!IsNoisy(r1)) {
        sigep1 = 1;
      }
      double sigep01 = sigep0 * sigep1;

      // Start by computing the time RBF (A matrix)
      double a01 = exp(-pow(tau_eta * time_diff, 2) / 2);
      // Next compute dA/dtau_eta
      double dadt = a01 * (-tau_eta * pow(time_diff, 2));
      // Add nugget_eta to A if applicable
      if (time_diff <= datum::eps) {
        a01 += nugget_eta;
      }
      double rho01 = r0 == r1 ? 1 : rho; // if same region, no rho
      double sqrt_keta = sqrt(k_eta[r0] * k_eta[r1]);
      double rho_keta = rho01 * sqrt_keta;

      dcovmat(i0, i1, tau_eta_slice) = rho_keta * dadt * sigep01;

      if (time_diff <= datum::eps) {
        dcovmat(i0, i1, nugget_eta_slice) = rho_keta * sigep01;
      }

      if (r0 != r1) {
        // r0 and r1 belong to different regions
        dcovmat(i0, i1, rho_slice) = sqrt_keta * a01 * sigep01;
        dcovmat(i0, i1, k_eta0_slice) =
            0.5 * rho * sqrt(k_eta[r1]) * a01 * sigep01 / sqrt(k_eta[r0]);
        dcovmat(i0, i1, k_eta1_slice) =
            0.5 * rho * sqrt(k_eta[r0]) * a01 * sigep01 / sqrt(k_eta[r1]);
      } else {
        // r0 and r1 in same region: d_rho = 0 and d_keta = 0 for other region
        uword active_k_eta_slice = r0 == 0 ? k_eta0_slice : k_eta1_slice;
        dcovmat(i0, i1, active_k_eta_slice) = a01 * sigep01;
      }

      // Symmetrize each slice
      for (uint islice = 0; islice < covparms.n_elem; islice++) {
        dcovmat(i1, i0, islice) = dcovmat(i0, i1, islice);
      }
    }
  }
  return dcovmat;
}

//' Compute QFunc Covariance Matrix.
//'
//' @param stage1_parms A 2 x 5 matrix of stage 1 paramerters.
//' @param covparms 5-vector of stage1 parameters.
//' @param locs A matrix with \code{n} rows and 5 columns.
// [[Rcpp::export]]
arma::mat test_qfunc_cov(arma::mat stage1_parms, arma::vec covparms,
                         arma::mat locs) {
  QFuncCovariance qfunc_cov(stage1_parms);
  return qfunc_cov(covparms, locs);
}

//' Compute Derivative of QFunc Covariance Matrix.
//'
//' @param stage1_parms A 2 x 5 matrix of stage 1 paramerters.
//' @param covparms 5-vector of stage1 parameters.
//' @param locs A matrix with \code{n} rows and 5 columns.
// [[Rcpp::export]]
arma::cube test_d_qfunc_cov(arma::mat stage1_parms, arma::vec covparms,
                            arma::mat locs) {
  DQFuncCovariance d_qfunc_cov(stage1_parms);
  return d_qfunc_cov(covparms, locs);
}