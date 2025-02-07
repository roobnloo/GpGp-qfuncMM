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
