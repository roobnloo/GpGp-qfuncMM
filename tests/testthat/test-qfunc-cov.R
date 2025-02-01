test_that("QFunc Covariance Matrix", {
  locs <- matrix(
    c(
      21, 21, 10, 1, 0,
      21, 23, 10, 1, 0,
      15, 4, 11, 1, 1,
      16, 3, 11, 2, 1
    ),
    ncol = 5,
    byrow = TRUE
  )

  s1 <- matrix(
    c(
      1, 0.1, 0.3, 2, 1,
      1.1, 0.1, 0.25, 1, 0.9
    ),
    ncol = 5
  )

  s2 <- c(0.6, 1, 2, 0.5, 1.1)

  covmx <- round(test_qfunc_cov(s1, s2, locs), 2)

  expect_c <- matrix(
    c(
      4.40, 3.36, 1.60, 0.67,
      3.36, 4.40, 1.60, 0.67,
      1.60, 1.60, 5.91, 1.47,
      0.67, 0.67, 1.47, 5.91
    ),
    ncol = 4, byrow = TRUE
  )
  expect_equal(covmx, expect_c)

  dcov <- test_d_qfunc_cov(s1, s2, locs)
  expect_equal(dim(dcov), c(4, 4, 5))
})
