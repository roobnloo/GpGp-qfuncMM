test_that("Exclude dimensions in ordering and neighbors", {
  n1 <- 15
  n2 <- 15
  n <- n1 * n2
  locs <- as.matrix(expand.grid(1:n1, 1:n2))

  set.seed(123)
  ord <- order_maxmin(locs)
  set.seed(123)
  ord_exclude <- order_maxmin(locs, exclude_dims = 2)

  expect_true(sum(abs(ord - ord_exclude)) > 0)

  locsord <- locs[ord, ]
  m <- 20

  set.seed(456)
  nnarray <- find_ordered_nn(locsord, m = m)
  set.seed(456)
  nnarray_exclude_null <- find_ordered_nn(locsord, m = m, exclude_dims = NULL)
  expect_equal(nnarray[n, ], nnarray_exclude_null[n, ])
  set.seed(456)
  nnarray_exclude <- find_ordered_nn(locsord, m = m, exclude_dims = 1)
  expect_true(sum(abs(nnarray[n, ] - nnarray_exclude[n, ])) > 0)
})
