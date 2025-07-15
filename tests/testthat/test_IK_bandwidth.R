test_that("IK bandwidth works as expected", {
  n = 1000; threshold = 2.4
  X = runif(n, -7, 12)
  W = as.numeric(X >= threshold)
  Y = (1 + 2*W)*(1 + X^2) + 1 / (1 + exp(X)) + rnorm(n, sd = .5)


  kernel = "triangular"
  expect_equal(IK_bandwidth(Y, X, threshold, kernel = kernel),
              IK_bandwidth(-Y, X, threshold, kernel = kernel))
  expect_equal(IK_bandwidth(Y, X, threshold, kernel = kernel),
              IK_bandwidth(Y, -X, -threshold, kernel = kernel))
  expect_equal(IK_bandwidth(Y, X, threshold, kernel = kernel),
              IK_bandwidth(Y + 42, X, threshold, kernel = kernel))

  kernel = "uniform"
  expect_equal(IK_bandwidth(Y, X, threshold, kernel = kernel),
              IK_bandwidth(-Y, X, threshold, kernel = kernel))
  expect_equal(IK_bandwidth(Y, X, threshold, kernel = kernel),
              IK_bandwidth(Y, -X, -threshold, kernel = kernel))
  expect_equal(IK_bandwidth(Y, X, threshold, kernel = kernel),
              IK_bandwidth(Y + 42, X, threshold, kernel = kernel))
})
