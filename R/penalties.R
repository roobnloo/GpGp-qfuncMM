# penalty functions


#' expit function and integral of expit function
#'
#' @noRd
expit <- function(x) {
  stats::plogis(x)
}

intexpit <- function(x) {
  matrixStats::logSumExp(c(0, x))
}

#' penalize large values of parameter: penalty, 1st deriative, 2nd derivative
#'
#' @param x argument to penalty
#' @param tt scale parameter of penalty
#' @param aa location parameter of penalty
#' @noRd
pen_hi <- function(x, tt, aa) {
  -tt * intexpit(x - aa)
}

dpen_hi <- function(x, tt, aa) {
  -tt * expit(x - aa)
}

ddpen_hi <- function(x, tt, aa) {
  -tt * expit(x - aa) / (1 + exp(x - aa))
}

#' penalize small values of parameter: penalty, 1st deriative, 2nd derivative
#'
#' @param x argument to penalty
#' @param tt scale parameter of penalty
#' @param aa location parameter of penalty
#' @noRd
pen_lo <- function(x, tt, aa) {
  -tt * intexpit(-x + aa)
}

dpen_lo <- function(x, tt, aa) {
  +tt * expit(-x + aa)
}

ddpen_lo <- function(x, tt, aa) {
  -tt * expit(-x + aa) / (1 + exp(-x + aa))
}


#' penalize small values of log parameter: penalty, 1st deriative, 2nd derivative
#'
#' @param x argument to penalty
#' @param tt scale parameter of penalty
#' @param aa location parameter of penalty
#' @noRd
pen_loglo <- function(x, tt, aa) {
  if (x == 0) {
    return(0.0)
  } else {
    return(pen_lo(log(x), tt, aa))
  }
}

dpen_loglo <- function(x, tt, aa) {
  if (x == 0) {
    return(0.0)
  } else {
    return(dpen_lo(log(x), tt, aa) / x)
  }
}


ddpen_loglo <- function(x, tt, aa) {
  if (x == 0) {
    return(0.0)
  } else {
    return(ddpen_lo(log(x), tt, aa) / x^2 - dpen_lo(log(x), tt, aa) / x^2)
  }
}
