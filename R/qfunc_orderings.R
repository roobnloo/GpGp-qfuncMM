#' Ordering for Qfunc mixed model
#'
#' First the time dimension is ordered middle-out. For each time point, the voxels are
#' ordered using maximum minimum distance ordering independently for each region.
#' Finally, these within-region orderings are interleaved.
#'
#' @inheritParams order_dist_to_point
#' @param space_cols column indices of spatial dimensions of \code{locs}
#' @param time_col column index of time dimension of \code{locs}
#' @param region_col column index of region dimension of \code{locs}
#' @return A vector of indices giving the ordering.
#'
#' @export
order_qfuncmm <- function(locs, space_cols = seq(1L, 3L), time_col = 4L, region_col = 5L) {
  m <- max(locs[, time_col])
  n <- nrow(locs)
  mid_time <- ceiling(m / 2)
  left_time <- seq(mid_time, 1, by = -1)
  right_time <- seq(mid_time + 1, m, by = 1)
  if (m %% 2 == 1) {
    right_time <- c(right_time, NA)
  }
  time_order <- as.numeric(rbind(left_time, right_time))
  time_order <- time_order[!is.na(time_order)]

  ord <- integer(n)
  index <- 1
  for (t in time_order) {
    # Get indices for the current time
    time_indices <- which(locs[, time_col] == t)

    # Separate indices by region
    region0_indices <- time_indices[locs[time_indices, region_col] == 0L]
    region1_indices <- time_indices[locs[time_indices, region_col] == 1L]

    # Order voxels within each region using order_maxmin
    region0_order <- order_maxmin(locs[region0_indices, space_cols], space_time = FALSE)
    region1_order <- order_maxmin(locs[region1_indices, space_cols], space_time = FALSE)

    # Interleave the indices
    max_length <- max(length(region0_order), length(region1_order))
    interleaved_order <- matrix(NA, 2, max_length)
    interleaved_order[1, seq_along(region0_order)] <- region0_indices[region0_order]
    interleaved_order[2, seq_along(region1_order)] <- region1_indices[region1_order]
    interleaved_order <- as.numeric(interleaved_order)
    interleaved_order <- interleaved_order[!is.na(interleaved_order)]

    # Add to the final order
    ord[index:(index + length(interleaved_order) - 1)] <- interleaved_order
    index <- index + length(interleaved_order)
  }

  return(ord)
}

#' @inheritParams order_qfuncmm
#'
#' @export
order_region_middleout <- function(locs, space_cols = seq(1L, 3L), time_col = 4L, region_col = 5L) {
  m <- max(locs[, time_col])
  n <- nrow(locs)
  mid_time <- ceiling(m / 2)
  left_time <- seq(mid_time, 1, by = -1)
  right_time <- seq(mid_time + 1, m, by = 1)
  if (m %% 2 == 1) {
    right_time <- c(right_time, NA)
  }
  time_order <- as.numeric(rbind(left_time, right_time))
  time_order <- time_order[!is.na(time_order)]

  ord <- integer(n)
  index <- 1
  for (t in time_order) {
    # Get indices for the current time
    time_indices <- which(locs[, time_col] == t)

    # Separate indices by region
    region0_indices <- time_indices[locs[time_indices, region_col] == 0L]
    region1_indices <- time_indices[locs[time_indices, region_col] == 1L]

    # Order voxels within each region using order_maxmin
    region0_order <- order_middleout(locs[region0_indices, space_cols])
    region1_order <- order_middleout(locs[region1_indices, space_cols])

    # Interleave the indices
    max_length <- max(length(region0_order), length(region1_order))
    interleaved_order <- matrix(NA, 2, max_length)
    interleaved_order[1, seq_along(region0_order)] <- region0_indices[region0_order]
    interleaved_order[2, seq_along(region1_order)] <- region1_indices[region1_order]
    interleaved_order <- as.numeric(interleaved_order)
    interleaved_order <- interleaved_order[!is.na(interleaved_order)]

    # Add to the final order
    ord[index:(index + length(interleaved_order) - 1)] <- interleaved_order
    index <- index + length(interleaved_order)
  }

  return(ord)
}

#' @inheritParams order_qfuncmm
#'
#' @export
completely_random_order <- function(locs, space_cols = seq(1L, 3L), time_col = 4L, region_col = 5L) {
  sample(nrow(locs))
}

#' Spatial Block Configuration
#'
#' Orders points in blocks to maintain spatial correlation while allowing
#' temporal correlation through alternation
#'
#' @inheritParams order_qfuncmm
#' @param block_size Size of spatial blocks to alternate between regions
#' @return A vector of indices giving the ordering
#' @export
order_spatial_blocks <- function(locs, space_cols = seq(1L, 3L),
                                 time_col = 4L, region_col = 5L,
                                 block_size = 10) {
  # Get dimensions
  m <- max(locs[, time_col])
  n <- nrow(locs)

  # Order time points middle out
  mid_time <- ceiling(m / 2)
  left_time <- seq(mid_time, 1, by = -1)
  right_time <- seq(mid_time + 1, m, by = 1)
  if (m %% 2 == 1) {
    right_time <- c(right_time, NA)
  }
  time_order <- as.numeric(rbind(left_time, right_time))
  time_order <- time_order[!is.na(time_order)]

  ord <- integer(n)
  index <- 1

  for (t in time_order) {
    # Get indices for current time
    time_indices <- which(locs[, time_col] == t)

    # Separate by region
    region0_indices <- time_indices[locs[time_indices, region_col] == 0L]
    region1_indices <- time_indices[locs[time_indices, region_col] == 1L]

    # Order spatially within each region
    region0_spatial <- order_maxmin(locs[region0_indices, space_cols])
    region1_spatial <- order_maxmin(locs[region1_indices, space_cols])

    # Break into blocks and interleave
    n0 <- length(region0_spatial)
    n1 <- length(region1_spatial)

    block_starts0 <- seq(1, n0, by = block_size)
    block_starts1 <- seq(1, n1, by = block_size)

    # Interleave blocks
    for (i in seq_len(max(length(block_starts0), length(block_starts1)))) {
      if (i <= length(block_starts0)) {
        end0 <- min(block_starts0[i] + block_size - 1, n0)
        block0 <- region0_indices[region0_spatial[block_starts0[i]:end0]]
        ord[index:(index + length(block0) - 1)] <- block0
        index <- index + length(block0)
      }

      if (i <= length(block_starts1)) {
        end1 <- min(block_starts1[i] + block_size - 1, n1)
        block1 <- region1_indices[region1_spatial[block_starts1[i]:end1]]
        ord[index:(index + length(block1) - 1)] <- block1
        index <- index + length(block1)
      }
    }
  }

  return(ord)
}

#' Time-Stratified Configuration
#'
#' Groups time points into strata and maintains spatial coherence within strata
#'
#' @inheritParams order_qfuncmm
#' @param stratum_size Number of time points per stratum
#' @return A vector of indices giving the ordering
#' @export
order_time_strata <- function(locs, space_cols = seq(1L, 3L),
                              time_col = 4L, region_col = 5L,
                              stratum_size = 3) {
  # Get dimensions
  m <- max(locs[, time_col])
  n <- nrow(locs)

  # Create strata from middle outward
  mid_time <- ceiling(m / 2)
  strata <- list()
  stratum_idx <- 1

  left <- mid_time
  right <- mid_time + 1

  while (left > 0 || right <= m) {
    # Build current stratum
    stratum <- numeric()
    for (i in 1:stratum_size) {
      if (left > 0) {
        stratum <- c(stratum, left)
        left <- left - 1
      }
      if (right <= m) {
        stratum <- c(stratum, right)
        right <- right + 1
      }
    }
    strata[[stratum_idx]] <- stratum
    stratum_idx <- stratum_idx + 1
  }

  # Order within each stratum
  ord <- integer(n)
  index <- 1

  for (stratum in strata) {
    # Get all points in current stratum
    stratum_indices <- which(locs[, time_col] %in% stratum)

    # Order region 1 then region 2
    for (r in c(0L, 1L)) {
      region_indices <- stratum_indices[locs[stratum_indices, region_col] == r]
      if (length(region_indices) > 0) {
        # Order spatially within region
        spatial_order <- order_maxmin(locs[region_indices, space_cols])
        ord[index:(index + length(region_indices) - 1)] <-
          region_indices[spatial_order]
        index <- index + length(region_indices)
      }
    }
  }

  return(ord)
}

#' Progressive Resolution Configuration
#'
#' Orders points by starting with a coarse grid and progressively filling in
#' @inheritParams order_qfuncmm
#' @param stride Initial stride for coarse grid
#' @return A vector of indices giving the ordering
#' @export
order_progressive <- function(locs, space_cols = seq(1L, 3L),
                              time_col = 4L, region_col = 5L,
                              stride = 3) {
  m <- max(locs[, time_col])
  n <- nrow(locs)

  # Order time points middle out
  mid_time <- ceiling(m / 2)
  left_time <- seq(mid_time, 1, by = -1)
  right_time <- seq(mid_time + 1, m, by = 1)
  time_order <- as.numeric(rbind(left_time, right_time))
  time_order <- time_order[!is.na(time_order)]

  ord <- integer(n)
  index <- 1

  # Process each resolution level
  for (current_stride in c(stride, 2, 1)) {
    for (t in time_order) {
      # Get indices for current time
      time_indices <- which(locs[, time_col] == t)

      # Skip points already ordered
      remaining_indices <- setdiff(time_indices, ord[1:(index - 1)])
      if (length(remaining_indices) == 0) next

      # Take every stride-th point
      for (r in c(0L, 1L)) {
        region_indices <- remaining_indices[locs[remaining_indices, region_col] == r]
        if (length(region_indices) == 0) next

        # Get unique spatial coordinates
        spatial_coords <- unique(locs[region_indices, space_cols])
        stride_indices <- seq(1, nrow(spatial_coords), by = current_stride)
        stride_coords <- spatial_coords[stride_indices, , drop = FALSE]

        # Match back to full indices
        for (i in seq_len(nrow(stride_coords))) {
          coord_match <- which(
            locs[region_indices, space_cols[1]] == stride_coords[i, 1] &
              locs[region_indices, space_cols[2]] == stride_coords[i, 2] &
              locs[region_indices, space_cols[3]] == stride_coords[i, 3]
          )
          if (length(coord_match) > 0) {
            ord[index] <- region_indices[coord_match]
            index <- index + 1
          }
        }
      }
    }
  }

  return(ord)
}

#' Distance-Based Hybrid Configuration
#'
#' Uses MMD for first portion then forces region alternation
#' @inheritParams order_qfuncmm
#' @param mmd_fraction Fraction of points to order using pure MMD
#' @param spatial_scale Scale factor for spatial coordinates
#' @return A vector of indices giving the ordering
#' @export
order_hybrid <- function(locs, space_cols = seq(1L, 3L),
                         time_col = 4L, region_col = 5L,
                         mmd_fraction = 0.7,
                         spatial_scale = 15) {
  n <- nrow(locs)
  n_mmd <- floor(n * mmd_fraction)

  # Scale coordinates for MMD portion
  scaled_locs <- locs
  scaled_locs[, space_cols] <- scaled_locs[, space_cols] / spatial_scale

  # Do MMD on scaled coordinates for first portion
  all_coords <- scaled_locs[, c(space_cols, time_col)]
  mmd_order <- order_maxmin(all_coords)
  ord <- mmd_order[1:n_mmd]

  # For remaining points, alternate regions
  remaining <- setdiff(1:n, ord)
  remaining_locs <- locs[remaining, ]

  region0 <- remaining[remaining_locs[, region_col] == 0L]
  region1 <- remaining[remaining_locs[, region_col] == 1L]

  # Order spatially within regions
  if (length(region0) > 0) {
    region0 <- region0[order_maxmin(locs[region0, space_cols])]
  }
  if (length(region1) > 0) {
    region1 <- region1[order_maxmin(locs[region1, space_cols])]
  }

  # Interleave remaining points
  i0 <- 1
  i1 <- 1
  while (i0 <= length(region0) || i1 <= length(region1)) {
    if (i0 <= length(region0)) {
      ord <- c(ord, region0[i0])
      i0 <- i0 + 1
    }
    if (i1 <= length(region1)) {
      ord <- c(ord, region1[i1])
      i1 <- i1 + 1
    }
  }

  return(ord)
}
