# Helper function to compute RBF kernel
rbf <- function(tau, time_sqrd_mat) {
  exp(-tau^2 / 2 * time_sqrd_mat)
}

# Objective function for temporal covariance optimization
temporal_obj <- function(theta, data_cross, time_sqrd, xi_tim, tau_gamma) {
  chi <- theta[1] # chi_tim parameter
  tau <- theta[2] # tau_eta parameter

  # Model temporal covariance
  model_cov <- chi * rbf(tau, time_sqrd) + xi_tim * rbf(tau_gamma, time_sqrd)

  # Sum squared differences between empirical and model covariances
  sum(sapply(data_cross, \(slice) {
    diffmx <- slice - model_cov
    diag(diffmx) <- 0L # Ignore diagonal elements
    sum(diffmx^2)
  }))
}

# Initialize parameters for a single region
init_region <- function(region_info) {
  s1 <- unlist(region_info$stage1) # Stage 1 estimates
  data <- region_info$data_std # Standardized data
  n_time <- nrow(data)
  n_voxel <- ncol(data)

  # Get Stage 1 estimates
  xi_tot <- s1["k_gamma"] + s1["nugget_gamma"]
  xi_tim <- s1["k_gamma"]
  tau_gamma <- s1["tau_gamma"]
  sigma2_ep <- 1

  if (region_info$cov_setting == "noisy") {
    xi_tot <- s1["sigma2_ep"] * (xi_tot + 1)
    xi_tim <- s1["sigma2_ep"] * xi_tim
    sigma2_ep <- s1["sigma2_ep"]
  }

  # Compute total empirical variance
  mu_hat <- mean(data) # Overall mean
  chi_tot <- max(1e-2, sum((data - mu_hat)^2) / (n_time * n_voxel - 1) - xi_tot)

  # Compute temporal covariances for optimization
  time_sqrd <- outer(1:n_time, 1:n_time, `-`)^2
  data_cross <- lapply(1:n_voxel, \(l) {
    # Center the data for this voxel
    data_l <- data[, l] - mu_hat
    # Compute empirical covariances
    tcrossprod(data_l)
  })

  # Optimize temporal parameters
  optim_result <- stats::optim(
    c(1, 1), # Initial values for chi_tim and tau_eta
    temporal_obj,
    data_cross = data_cross,
    time_sqrd = time_sqrd,
    xi_tim = xi_tim,
    tau_gamma = tau_gamma,
    method = "L-BFGS-B",
    lower = c(0.01, 0.01)
  )

  chi_tim <- optim_result$par[1]
  tau_eta <- optim_result$par[2]

  # Compute k_eta based on setting
  k_eta <- if (region_info$cov_setting == "noisy") {
    chi_tot / (sigma2_ep + 1e-2)
  } else {
    chi_tot
  }

  # Compute nugget_eta using the ratio of total to temporal variance
  nugget_eta <- max(1e-2, chi_tot / chi_tim - 1)

  result <- c(k_eta = k_eta, tau_eta = tau_eta, nugget_eta = nugget_eta)
  names(result) <- c("k_eta", "tau_eta", "nugget_eta")
  return(result)
}

# Main stage 2 initialization function
stage2_init <- function(r1_info, r2_info) {
  # Get initial estimates for each region
  r1k <- init_region(r1_info)
  r2k <- init_region(r2_info)

  # Use correlation of EBLUEs for rho
  rho_eblue <- cor(r1_info$eblue, r2_info$eblue)

  # Average the tau parameters from the two regions
  tau_eta <- mean(c(r1k["tau_eta"], r2k["tau_eta"]))

  # Take maximum of nugget parameters
  nugget_eta <- max(r1k["nugget_eta"], r2k["nugget_eta"])

  # Return parameter vector
  init <- c(
    rho = rho_eblue,
    k_eta1 = r1k["k_eta"],
    k_eta2 = r2k["k_eta"],
    tau_eta = tau_eta,
    nugget_eta = nugget_eta
  )
  names(init) <- c("rho", "k_eta1", "k_eta2", "tau_eta", "nugget_eta")

  return(init)
}
