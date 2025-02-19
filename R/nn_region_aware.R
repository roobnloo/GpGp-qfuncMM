interleave <- function(x, y) {
  if (length(x) < length(y)) {
    x <- c(x, rep(NA, length(y) - length(x)))
  } else if (length(y) < length(x)) {
    y <- c(y, rep(NA, length(x) - length(y)))
  }
  xy <- as.numeric(rbind(x, y))
  xy[!is.na(xy)]
}

find_ordered_nn_region_aware <- function(locs, m, st_scale, region_col) {
  # convert locs to matrix (works on vectors and data frames)
  locs <- as.matrix(locs)
  regionid <- locs[, region_col]
  locs <- locs[, -region_col, drop = FALSE]

  # number of locations
  n <- nrow(locs)
  m <- min(m, n - 1)
  mult <- 2

  # FNN::get.knnx has strange behavior for exact matches
  # so add a small amount of noise to each location
  ee <- min(apply(locs, 2, stats::sd))
  locs <- locs + matrix(ee * 1e-4 * stats::rnorm(n * ncol(locs)), n, ncol(locs))


  d <- ncol(locs) - 1
  locs[, 1:d] <- locs[, 1:d] / st_scale[1]
  locs[, d + 1] <- locs[, d + 1] / st_scale[2]

  # to store the nearest neighbor indices
  NNarray <- matrix(NA, n, m + 1)

  # to the first mult*m+1 by brutce force
  maxval <- min(mult * m + 1, n)
  NNarray[1:maxval, ] <- find_ordered_nn_brute(locs[1:maxval, , drop = FALSE], m)

  query_inds <- min(maxval + 1, n):n
  data_inds <- 1:n

  msearch <- m

  while (length(query_inds) > 0) {
    msearch <- min(max(query_inds), 2 * msearch)
    data_inds <- 1:min(max(query_inds), n)
    data_inds_region <- cbind(seq_len(length(data_inds)), regionid[data_inds])
    data_inds_r0 <- data_inds[regionid[data_inds] == 0]
    data_inds_r1 <- data_inds[regionid[data_inds] == 1]

    msearch0 <- min(msearch, length(data_inds_r0))
    msearch1 <- min(msearch, length(data_inds_r1))
    nn0 <- FNN::get.knnx(locs[data_inds_r0, , drop = FALSE], locs[query_inds, , drop = FALSE], msearch0)$nn.index
    nn1 <- FNN::get.knnx(locs[data_inds_r1, , drop = FALSE], locs[query_inds, , drop = FALSE], msearch1)$nn.index

    NN <- matrix(Inf, length(query_inds), msearch)
    for (i in seq_along(query_inds)) {
      nn0i <- data_inds_region[data_inds_region[, 2] == 0, ][nn0[i, ]]
      nn1i <- data_inds_region[data_inds_region[, 2] == 1, ][nn1[i, ]]
      nni <- NULL
      if (nn0i[1] == query_inds[i]) {
        nni <- interleave(nn0i, nn1i)
      } else {
        nni <- interleave(nn1i, nn0i)
      }
      NN[i, ] <- nni[1:msearch]
    }

    # NN <- FNN::get.knnx(locs[data_inds, , drop = FALSE], locs[query_inds, , drop = FALSE], msearch)$nn.index
    less_than_k <- t(sapply(1:nrow(NN), function(k) NN[k, ] <= query_inds[k]))
    sum_less_than_k <- apply(less_than_k, 1, sum)
    ind_less_than_k <- which(sum_less_than_k >= m + 1)

    NN_m <- t(sapply(ind_less_than_k, function(k) NN[k, ][less_than_k[k, ]][1:(m + 1)]))

    NNarray[query_inds[ind_less_than_k], ] <- NN_m

    query_inds <- query_inds[-ind_less_than_k]
  }

  return(NNarray)
}

#' Time-Block Neighbor Selection
#'
#' Selects neighbors by first prioritizing spatial neighbors within the same time point
#' and region, then adding temporal neighbors from nearby time points.
#'
#' @inheritParams find_ordered_nn_region_aware
#' @return Matrix of neighbor indices, where first column is point index
#' @export
find_ordered_nn_time_block <- function(locs, m, st_scale, region_col = 5L) {
  # Convert locs to matrix for consistent handling
  locs <- as.matrix(locs)

  # Get dimensions
  n <- nrow(locs)
  m <- min(m, n - 1) # Cannot have more neighbors than previous points
  mult <- 2 # Multiplier for initial kNN search

  # Add small noise to prevent exact matches causing FNN issues
  ee <- min(apply(locs, 2, stats::sd))

  # Scale coordinates
  d <- ncol(locs) - 2 # Subtract time and region columns
  locs[, -region_col] <- locs[, -region_col] + matrix(ee * 1e-4 * stats::rnorm(n * (d + 1)), n, d + 1)
  locs[, 1:d] <- locs[, 1:d] / st_scale[1]
  locs[, d + 1] <- locs[, d + 1] / st_scale[2]

  # Initialize neighbor array
  NNarray <- matrix(NA, n, m + 1)

  # Handle first maxval points using brute force
  maxval <- min(mult * m + 1, n)
  NNarray[1:maxval, ] <- find_ordered_nn_brute(locs[1:maxval, -region_col], m)

  # Process remaining points
  query_inds <- min(maxval + 1, n):n
  msearch <- m

  while (length(query_inds) > 0) {
    msearch <- min(max(query_inds), 2 * msearch)
    data_inds <- 1:min(max(query_inds), n)

    # For each query point, we'll:
    # 1. First find spatial neighbors from same time point and region
    # 2. Then find temporal neighbors from other times

    # Get spatial neighbors first
    spatial_matches <- matrix(Inf, length(query_inds), msearch)

    for (i in seq_along(query_inds)) {
      current_point <- query_inds[i]
      current_time <- locs[current_point, d + 1]
      current_region <- locs[current_point, region_col]

      # Find points from same time and region
      same_time_region <- which(
        abs(locs[1:(current_point), d + 1] - current_time) <= 1e-4 &
          locs[1:(current_point), region_col] == current_region
      )

      if (length(same_time_region) > 0) {
        # Get spatial neighbors from same time/region
        nn_spatial <- FNN::get.knnx(
          locs[same_time_region, 1:d, drop = FALSE],
          locs[current_point, 1:d, drop = FALSE],
          k = min(msearch, length(same_time_region))
        )
        spatial_matches[i, 1:ncol(nn_spatial$nn.index)] <-
          same_time_region[nn_spatial$nn.index]
      }
    }

    # Get temporal neighbors
    temporal_matches <- FNN::get.knnx(
      locs[data_inds, ],
      locs[query_inds, ],
      k = msearch
    )$nn.index

    # Combine spatial and temporal neighbors
    # Prefer spatial neighbors when available
    NN <- matrix(Inf, length(query_inds), msearch)
    for (i in seq_along(query_inds)) {
      spatial_valid <- spatial_matches[i, ] < Inf
      n_spatial <- sum(spatial_valid)

      if (n_spatial > 0) {
        # Take available spatial neighbors
        NN[i, 1:n_spatial] <- spatial_matches[i, spatial_valid]

        # Fill remaining spots with temporal neighbors
        remaining_spots <- (n_spatial + 1):msearch
        if (length(remaining_spots) > 0) {
          temp_neighbors <- setdiff(
            temporal_matches[i, ],
            spatial_matches[i, spatial_valid]
          )
          NN[i, remaining_spots] <- temp_neighbors[seq_along(remaining_spots)]
        }
      } else {
        # If no spatial neighbors, use all temporal
        NN[i, ] <- temporal_matches[i, ]
      }
    }

    # Check which points have enough valid neighbors
    less_than_k <- t(sapply(seq_len(nrow(NN)), function(k) NN[k, ] <= query_inds[k]))
    sum_less_than_k <- apply(less_than_k, 1, sum)
    ind_less_than_k <- which(sum_less_than_k >= m + 1)

    if (length(ind_less_than_k) > 0) {
      NN_m <- t(sapply(
        ind_less_than_k,
        function(k) NN[k, ][less_than_k[k, ]][1:(m + 1)]
      ))
      NNarray[query_inds[ind_less_than_k], ] <- NN_m
      query_inds <- query_inds[-ind_less_than_k]
    }
  }

  return(NNarray)
}


#' Pure Random Neighbor Selection
#'
#' Selects neighbors completely randomly from available previous points,
#' with no consideration of regions or space-time structure.
#' This serves as a baseline for comparing against more structured approaches.
#'
#' @inheritParams find_ordered_nn_region_aware
#' @return Matrix of neighbor indices, where first column is point index
#' @export
find_ordered_nn_pure_random <- function(locs, m, st_scale, region_col = 5L) {
  # Get dimensions
  n <- nrow(locs)
  m <- min(m, n - 1)
  mult <- 2

  # Initialize neighbor array
  NNarray <- matrix(NA, n, m + 1)

  # Handle first maxval points using brute force
  maxval <- min(mult * m + 1, n)
  NNarray[1:maxval, ] <- find_ordered_nn_brute(locs[1:maxval, ], m)

  # For remaining points, just randomly sample from available previous points
  for (i in (maxval + 1):n) {
    # Current point is always first
    NNarray[i, 1] <- i

    # Randomly sample m points from all available previous points
    available_points <- 1:(i - 1)
    n_needed <- min(m, length(available_points))

    NNarray[i, 2:(n_needed + 1)] <- sample(available_points, n_needed)
  }

  return(NNarray)
}

#' Stratified Temporal Neighbor Selection
#'
#' Selects neighbors by first stratifying available points by temporal distance,
#' then sampling within strata. This ensures coverage of both nearby and distant
#' time points while maintaining spatial correlation through sampling weights.
#'
#' @inheritParams find_ordered_nn_region_aware
#' @return Matrix of neighbor indices, where first column is point index
#' @export
find_ordered_nn_temporal_strata <- function(locs, m, st_scale, region_col = 5L) {
  # Convert locs to matrix for consistent handling
  locs <- as.matrix(locs)

  # Get dimensions
  n <- nrow(locs)
  m <- min(m, n - 1)
  mult <- 2

  # Add small noise to prevent exact matches causing FNN issues
  d <- ncol(locs) - 2 # Subtract time and region columns
  time_col <- d + 1
  ee <- min(apply(locs, 2, stats::sd))
  locs[, -region_col] <- locs[, -region_col] + matrix(ee * 1e-4 * stats::rnorm(n * time_col), n, time_col)

  # Initialize neighbor array
  NNarray <- matrix(NA, n, m + 1)

  # Handle first maxval points using brute force
  maxval <- min(mult * m + 1, n)
  NNarray[1:maxval, ] <- find_ordered_nn_brute(locs[1:maxval, ], m)

  # For remaining points, use stratified sampling

  for (i in (maxval + 1):n) {
    # Current point is always first
    NNarray[i, 1] <- i

    current_time <- locs[i, time_col]
    available_points <- 1:(i - 1)

    # Compute temporal distances
    time_diffs <- abs(locs[available_points, time_col] - current_time)

    # Create temporal strata
    # More neighbors from nearby times, fewer from distant times
    strata <- list(
      near = which(time_diffs <= 2),
      medium = which(time_diffs > 2 & time_diffs <= 5),
      far = which(time_diffs > 5)
    )

    # Allocate neighbors to strata (50% near, 30% medium, 20% far)
    n_needed <- min(m, length(available_points))
    allocations <- c(
      near = ceiling(0.5 * n_needed),
      medium = ceiling(0.3 * n_needed),
      far = n_needed - ceiling(0.5 * n_needed) - ceiling(0.3 * n_needed)
    )

    # Sample from each stratum
    selected <- integer(0)
    for (stratum in names(strata)) {
      stratum_points <- available_points[strata[[stratum]]]
      if (length(stratum_points) > 0) {
        n_take <- min(allocations[[stratum]], length(stratum_points))
        selected <- c(selected, sample(stratum_points, n_take))
      }
    }

    # If we still need more points, sample from anywhere
    if (length(selected) < n_needed) {
      remaining <- setdiff(available_points, selected)
      additional <- sample(remaining, n_needed - length(selected))
      selected <- c(selected, additional)
    }

    # Randomly permute final selection
    selected <- sample(selected)
    NNarray[i, 2:(length(selected) + 1)] <- selected
  }

  return(NNarray)
}

#' Adaptive Distance Neighbor Selection
#'
#' Selects neighbors using an adaptive distance metric that changes based on
#' the position in the ordering. Early points focus more on temporal distance
#' to establish cross-region correlation, while later points balance spatial
#' and temporal distances.
#'
#' @inheritParams find_ordered_nn_region_aware
#' @return Matrix of neighbor indices, where first column is point index
#' @export
find_ordered_nn_adaptive <- function(locs, m, st_scale, region_col = 5L) {
  # Convert locs to matrix for consistent handling
  locs <- as.matrix(locs)

  # Get dimensions
  n <- nrow(locs)
  m <- min(m, n - 1)
  mult <- 2

  # Add small noise to prevent exact matches causing FNN issues
  ee <- min(apply(locs, 2, stats::sd))
  locs[, -region_col] <- locs[, -region_col] +
    matrix(ee * 1e-4 * stats::rnorm(n * (ncol(locs) - 1)), n, ncol(locs) - 1)

  # Initialize neighbor array
  NNarray <- matrix(NA, n, m + 1)

  # Handle first maxval points using brute force
  maxval <- min(mult * m + 1, n)
  NNarray[1:maxval, ] <- find_ordered_nn_brute(locs[1:maxval, ], m)

  # Initialize search parameters for remaining points
  query_inds <- min(maxval + 1, n):n
  msearch <- m

  while (length(query_inds) > 0) {
    msearch <- min(max(query_inds), 2 * msearch)
    data_inds <- 1:min(max(query_inds), n)

    # Create scaled locations with adaptive weights
    scaled_locs <- locs[, -region_col, drop = FALSE] # Remove region column
    d <- ncol(scaled_locs) - 1 # Last column is time

    # Scale spatial and temporal coordinates differently for each query point
    for (i in seq_along(query_inds)) {
      progress <- query_inds[i] / n # Position in ordering (0 to 1)

      # Early in ordering: emphasize temporal distance
      # Later in ordering: more balanced spatial-temporal distance
      spatial_weight <- progress # Increases from 0 to 1
      temporal_weight <- 1 - 0.5 * progress # Decreases from 1 to 0.5

      scaled_locs[, 1:d] <- locs[, 1:d] / st_scale[1] * spatial_weight
      scaled_locs[, d + 1] <- locs[, d + 1] / st_scale[2] * temporal_weight
    }

    # Find nearest neighbors using adaptive distances
    NN <- FNN::get.knnx(scaled_locs[data_inds, ],
      scaled_locs[query_inds, ],
      k = msearch
    )$nn.index

    # Process valid neighbors
    less_than_k <- t(sapply(seq_len(nrow(NN)), function(k) NN[k, ] <= query_inds[k]))
    sum_less_than_k <- apply(less_than_k, 1, sum)
    ind_less_than_k <- which(sum_less_than_k >= m + 1)

    if (length(ind_less_than_k) > 0) {
      NN_m <- t(sapply(
        ind_less_than_k,
        function(k) NN[k, ][less_than_k[k, ]][1:(m + 1)]
      ))
      NNarray[query_inds[ind_less_than_k], ] <- NN_m
      query_inds <- query_inds[-ind_less_than_k]
    }
  }

  return(NNarray)
}
