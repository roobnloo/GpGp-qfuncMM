#define BOOST_DISABLE_ASSERTS

#include "onepass.h"
#include <RcppArmadillo.h>
#include <math.h>
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(BH)]]

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
List vecchia_profbeta_loglik_grad_info(
    NumericVector covparms, StringVector covfun_name, NumericVector y,
    NumericMatrix X, const NumericMatrix locs, NumericMatrix NNarray,
    NumericMatrix additional_info = NumericMatrix()) {

  NumericVector ll(1);
  NumericVector grad(covparms.length());
  NumericVector betahat(X.ncol());
  NumericMatrix info(covparms.length(), covparms.length());
  NumericMatrix betainfo(X.ncol(), X.ncol());

  // this function calls arma_onepass_compute_pieces
  // then synthesizes the result into loglik, beta, grad, info, betainfo
  synthesize(covparms, covfun_name, locs, NNarray, y, X, &ll, &betahat, &grad,
             &info, &betainfo, true, true, additional_info);

  List ret = List::create(Named("loglik") = ll, Named("betahat") = betahat,
                          Named("grad") = grad, Named("info") = info,
                          Named("betainfo") = betainfo);
  return ret;
}

// [[Rcpp::export]]
List vecchia_profbeta_loglik(NumericVector covparms, StringVector covfun_name,
                             NumericVector y, NumericMatrix X,
                             const NumericMatrix locs, NumericMatrix NNarray) {

  NumericVector ll(1);
  NumericVector grad(covparms.length());
  NumericVector betahat(X.ncol());
  NumericMatrix info(covparms.length(), covparms.length());
  NumericMatrix betainfo(X.ncol(), X.ncol());

  // this function calls arma_onepass_compute_pieces
  // then synthesizes the result into loglik, beta, grad, info, betainfo
  synthesize(covparms, covfun_name, locs, NNarray, y, X, &ll, &betahat, &grad,
             &info, &betainfo, true, false);

  List ret = List::create(Named("loglik") = ll, Named("betahat") = betahat,
                          Named("betainfo") = betainfo);
  return ret;
}

List vecchia_meanzero_loglik(NumericVector covparms, StringVector covfun_name,
                             NumericVector y, const NumericMatrix locs,
                             NumericMatrix NNarray) {

  NumericMatrix X(1, 1);
  NumericVector ll(1);
  NumericVector grad(covparms.length());
  NumericVector betahat(X.ncol());
  // NumericVector betahat;
  NumericMatrix info(covparms.length(), covparms.length());
  NumericMatrix betainfo(X.ncol(), X.ncol());
  // NumericMatrix betainfo;

  // this function calls arma_onepass_compute_pieces
  // then synthesizes the result into loglik, beta, grad, info, betainfo
  synthesize(covparms, covfun_name, locs, NNarray, y, X, &ll, &betahat, &grad,
             &info, &betainfo, false, false);

  List ret = List::create(Named("loglik") = ll);
  return ret;
}

// [[Rcpp::export]]
List vecchia_grouped_profbeta_loglik_grad_info(
    NumericVector covparms, StringVector covfun_name, NumericVector y,
    NumericMatrix X, const NumericMatrix locs, List NNlist,
    NumericMatrix additional_info = NumericMatrix()) {

  NumericVector ll(1);
  NumericVector grad(covparms.length());
  NumericVector betahat(X.ncol());
  NumericMatrix info(covparms.length(), covparms.length());
  NumericMatrix betainfo(X.ncol(), X.ncol());

  // this function calls arma_onepass_compute_pieces
  // then synthesizes the result into loglik, beta, grad, info, betainfo
  // maybe the synthesize functions should take in an argument that
  // says which compute_pieces function to use
  synthesize_grouped(covparms, covfun_name, locs, NNlist, y, X, &ll, &betahat,
                     &grad, &info, &betainfo, true, true, additional_info);

  List ret = List::create(Named("loglik") = ll, Named("betahat") = betahat,
                          Named("grad") = grad, Named("info") = info,
                          Named("betainfo") = betainfo);
  return ret;
}

// [[Rcpp::export]]
List vecchia_grouped_profbeta_loglik(NumericVector covparms,
                                     StringVector covfun_name, NumericVector y,
                                     NumericMatrix X, const NumericMatrix locs,
                                     List NNlist) {

  NumericVector ll(1);
  NumericVector grad(covparms.length());
  NumericVector betahat(X.ncol());
  NumericMatrix info(covparms.length(), covparms.length());
  NumericMatrix betainfo(X.ncol(), X.ncol());

  // this function calls arma_onepass_compute_pieces
  // then synthesizes the result into loglik, beta, grad, info, betainfo
  // maybe the synthesize functions should take in an argument that
  // says which compute_pieces function to use
  synthesize_grouped(covparms, covfun_name, locs, NNlist, y, X, &ll, &betahat,
                     &grad, &info, &betainfo, true, false);

  List ret = List::create(Named("loglik") = ll, Named("betahat") = betahat,
                          Named("betainfo") = betainfo);
  return ret;
}

// [[Rcpp::export]]
List vecchia_grouped_meanzero_loglik(NumericVector covparms,
                                     StringVector covfun_name, NumericVector y,
                                     const NumericMatrix locs, List NNlist) {

  NumericMatrix X(1, 1);
  NumericVector ll(1);
  NumericVector grad(covparms.length());
  NumericVector betahat(X.ncol());
  // NumericVector betahat;
  NumericMatrix info(covparms.length(), covparms.length());
  NumericMatrix betainfo(X.ncol(), X.ncol());
  // NumericMatrix betainfo;

  // this function calls arma_onepass_compute_pieces
  // then synthesizes the result into loglik, beta, grad, info, betainfo
  // maybe the synthesize functions should take in an argument that
  // says which compute_pieces function to use
  synthesize_grouped(covparms, covfun_name, locs, NNlist, y, X, &ll, &betahat,
                     &grad, &info, &betainfo, false, false);

  List ret = List::create(Named("loglik") = ll);
  return ret;
}
