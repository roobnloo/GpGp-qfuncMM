#' Fit the QFunC Mixed Effects Model
#'
#' Given a response, set of locations, (optionally) a design matrix,
#' and a specified covariance function, return the maximum
#' Vecchia likelihood estimates, obtained with a Fisher scoring algorithm.
#'
#' @param region1 \eqn{M\times L_1} matrix for region1
#' @param region2 \eqn{M\times L_2} matrix for region2
#' @param r1_coords \eqn{L_1\times 3} matrix of coordinates for region1
#' @param r2_coords \eqn{L_2\times 3} matrix of coordinates for region2
#' @param start_parms \eqn{5} vector of starting values
#' @param stage1_parms \eqn{2 \times 5} matrix of stage1 starting values
#' @inheritParams fit_model
fit_qfuncmm <- function(
    region1, region2, r1_coords, r2_coords, start_parms, stage1_parms,
    NNarray = NULL, reorder = TRUE, group = TRUE,
    m_seq = c(10, 30), max_iter = 40, fixed_parms = NULL,
    silent = FALSE, st_scale = NULL, convtol = 1e-4) {
  covfun_name <- "qfuncmm"
  stopifnot("nrow(region1) and nrow(region2) must be equal" = nrow(region1) == nrow(region2))
  stopifnot("ncol(region1) and nrow(r1_coords) must be equal" = ncol(region1) == nrow(r1_coords))
  stopifnot("ncol(region2) and nrow(r2_coords) must be equal" = ncol(region2) == nrow(r2_coords))
  stopifnot("length(start_parms) must be equal to 5" = length(start_parms) == 5)

  if (start_parms[1] < -1 || start_parms[1] > 1) {
    stop("start_parms[1] (rho) should be in the range [-1, 1]")
  }
  if (any(start_parms[-1] < 0)) {
    stop("start_parms[2:5] should be nonnegative")
  }

  ntime <- nrow(region1)
  l1 <- nrow(r1_coords)
  l2 <- nrow(r2_coords)
  nl1 <- ntime * l1
  nl2 <- ntime * l2

  # The region column is a bit of a hack in order to track the regions in the approximation.
  # We must be able to map voxels to regions in order to calculate the covariance matrix.
  # The region column is ignored (for now) when determining the ordering and the neighbors.
  locs1 <- cbind(
    r1_coords[rep(seq_len(l1), each = ntime), ],
    time = rep(1:ntime, times = l1),
    region = 0L
  )
  locs2 <- cbind(
    r2_coords[rep(seq_len(l2), each = ntime), ],
    time = rep(1:ntime, times = l2),
    region = 1L
  )
  region_loc_id <- 5L

  locs <- rbind(locs1, locs2)
  y <- c(region1, region2)
  n <- length(y)

  # Create design matrix
  X <- Matrix::bdiag(rep(1, nl1), rep(1, nl2))

  # TODO: handle fixed_parms
  start_parms <- get_qfuncmm_start_parms(start_parms)
  active <- rep(TRUE, length(start_parms)) # this says all variables are estimated

  # get link functions
  linkfuns <- get_qfuncmm_linkfun()
  link <- linkfuns$link
  dlink <- linkfuns$dlink
  invlink <- linkfuns$invlink
  invlink_startparms <- invlink(start_parms)

  penalty <- get_qfuncmm_penalty(y, X, locs, covfun_name)
  pen <- penalty$pen
  dpen <- penalty$dpen
  ddpen <- penalty$ddpen

  # get an ordering and reorder everything
  if (reorder) {
    if (!silent) cat("Reordering...")
    ord <- order_qfuncmm(locs)
    if (!silent) cat("Done \n")
  } else {
    ord <- 1:n
  }
  yord <- y[ord]
  locsord <- locs[ord, , drop = FALSE]
  Xord <- as.matrix(X[ord, , drop = FALSE])

  # get neighbor array if not provided
  if (is.null(NNarray)) {
    if (!silent) cat("Finding nearest neighbors...")
    NNarray <- find_ordered_nn_region_aware(
      locsord, max(m_seq), st_scale, region_loc_id
    )
    if (!silent) cat("Done \n")
  }

  # refine the estimates using m in m_seq
  for (i in seq_along(m_seq)) {
    m <- m_seq[i]
    if (group) {
      NNlist <- group_obs(NNarray[, 1:(m + 1)])
      likfun <- function(logparms) {
        lp <- rep(NA, length(start_parms))
        lp[active] <- logparms
        lp[!active] <- invlink_startparms[!active]

        likobj <- vecchia_grouped_profbeta_loglik_grad_info(
          link(lp), covfun_name, yord, Xord, locsord, NNlist, stage1_parms
        )
        likobj$loglik <- -likobj$loglik - pen(link(lp))
        likobj$grad <- -c(likobj$grad) * dlink(lp) -
          dpen(link(lp)) * dlink(lp)
        likobj$info <- likobj$info * outer(dlink(lp), dlink(lp)) -
          ddpen(link(lp)) * outer(dlink(lp), dlink(lp))
        likobj$grad <- likobj$grad[active]
        likobj$info <- likobj$info[active, active]
        return(likobj)
      }
    } else {
      likfun <- function(logparms) {
        lp <- rep(NA, length(start_parms))
        lp[active] <- logparms
        lp[!active] <- invlink_startparms[!active]

        likobj <- vecchia_profbeta_loglik_grad_info(
          link(lp), covfun_name, yord, Xord, locsord, NNarray[, 1:(m + 1)], stage1_parms
        )
        likobj$loglik <- -likobj$loglik - pen(link(lp))
        likobj$grad <- -c(likobj$grad) * dlink(lp) -
          dpen(link(lp)) * dlink(lp)
        likobj$info <- likobj$info * outer(dlink(lp), dlink(lp)) -
          ddpen(link(lp)) * outer(dlink(lp), dlink(lp))
        likobj$grad <- likobj$grad[active]
        likobj$info <- likobj$info[active, active]
        return(likobj)
      }
    }
    fit <- fisher_scoring(likfun, invlink(start_parms)[active],
      link,
      silent = silent, convtol = convtol, max_iter = max_iter
    )
    invlink_startparms[active] <- fit$logparms
    # start_parms[active] <- fit$covparms
    start_parms <- link(invlink_startparms)
    fit$loglik <- -fit$loglik - pen(start_parms)
    invlink_startparms <- invlink(start_parms)
  }

  # return fit and information used for predictions
  fit$covfun_name <- covfun_name
  # fit$covparms <- start_parms
  lp <- rep(NA, length(start_parms))
  lp[active] <- fit$logparms
  lp[!active] <- invlink_startparms[!active]
  fit$covparms <- link(lp)
  fit$y <- y
  fit$locs <- locs
  fit$X <- X
  names(fit$covparms) <- c("rho", "k_eta1", "k_eta2", "tau_eta", "nugget_eta")
  names(fit$logparms) <- names(fit$covparms)
  class(fit) <- "GpGp_fit"
  return(fit)
}


# get default starting values
get_qfuncmm_start_parms <- function(start_parms) {
  start_parms
}


# rho in (-Inf, Inf) to [0,1], the original scale
invsigmoid <- function(x, lower, upper) {
  lower + (upper - lower) * stats::plogis(x)
}

d_invsigmoid <- function(x, lower, upper) {
  (upper - lower) * stats::plogis(x) * (1 - stats::plogis(x))
}

# rho in [0,1] to (-Inf, Inf), the optimization scale
sigmoid <- function(x, lower, upper) {
  stats::qlogis((x - lower) / (upper - lower))
}


# get the qfuncmm link function
# parameter list: (rho, k_eta1, k_eta2, tau_eta, nugget_eta)
get_qfuncmm_linkfun <- function() {
  link <- \(x) c(invsigmoid(x[1], -1, 1), exp(x[2:5])) # link puts optimization parameters to the real scale
  dlink <- \(x) c(d_invsigmoid(x[1], -1, 1), exp(x[2:5]))
  invlink <- \(x) c(sigmoid(x[1], -1, 1), log(x[2:5])) # invlink puts parameters on the optimization scale

  return(list(
    link = link, dlink = dlink, invlink = invlink
  ))
}

get_qfuncmm_penalty <- function(y, X, locs, covfun_name) {
  # pen <- function(x) 0.0
  # dpen <- function(x) rep(0,length(x))
  # ddpen <- function(x) matrix(0,length(x),length(x))
  pen_nug <- function(x, j) {
    pen_loglo(exp(x[j]), .01, log(0.01))
  }
  dpen_nug <- function(x, j) {
    dpen <- rep(0, length(x))
    dpen[j] <- dpen_loglo(exp(x[j]), .01, log(0.01))
    return(dpen)
  }
  ddpen_nug <- function(x, j) {
    ddpen <- matrix(0, length(x), length(x))
    ddpen[j, j] <- ddpen_loglo(exp(x[j]), .01, log(0.01))
    return(ddpen)
  }
  pen <- \(x) {
    rho_logit <- stats::qlogis((x[1] + 1) / 2)
    prho <- -matrixStats::logSumExp(c(0, rho_logit^2 / 3 - 6))
    pketa1 <- pen_nug(x, 2)
    pketa2 <- pen_nug(x, 3)
    ptaueta <- pen_nug(x, 4)
    pnug <- pen_nug(x, 5)
    prho + pketa1 + pketa2 + ptaueta + pnug
  }
  dpen <- function(x) {
    rho <- x[1]
    rho_logit <- stats::qlogis((rho + 1) / 2)
    ex <- exp(rho_logit^2 / 3 - 6)
    dpenrho <- 2 * ex * rho_logit / (3 * (1 + ex))
    dpenrho <- dpenrho * (-2 / (rho^2 + 1)) # chain rule logit'(rho)
    dpenrho <- c(-dpenrho, rep(0, length(x) - 1))
    dpenketa1 <- dpen_nug(x, 2)
    dpenketa2 <- dpen_nug(x, 3)
    dpentaueta <- dpen_nug(x, 4)
    dpennug <- dpen_nug(x, 5)
    dpenrho + dpenketa1 + dpenketa2 + dpentaueta + dpennug
  }
  ddpen <- function(x) {
    rho <- x[1]

    log_term <- log(-(1 + rho) / (rho - 1))
    exp_term <- exp(log_term^2 / 3)
    e6 <- exp(6)
    numerator <- 8 * exp_term * (3 * (e6 + exp_term) + 3 * (e6 + exp_term) * rho * log_term + 2 * e6 * log_term^2)
    denominator <- 9 * (e6 + exp_term)^2 * (1 - rho)^2 * (1 + rho)^2

    ddrho <- matrix(0, length(x), length(x))
    ddrho[1, 1] <- -numerator / denominator
    ddketa1 <- ddpen_nug(x, 2)
    ddketa2 <- ddpen_nug(x, 3)
    ddtaueta <- ddpen_nug(x, 4)
    ddnug <- ddpen_nug(x, 5)
    ddrho + ddketa1 + ddketa2 + ddtaueta + ddnug
  }
  return(list(pen = pen, dpen = dpen, ddpen = ddpen))
}
