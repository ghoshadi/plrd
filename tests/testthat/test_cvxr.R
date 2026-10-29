test_that("cvxr primal agrees with the exact quadprog dual (no spline)", {
  set.seed(42)
  n = 1000; threshold = 0
  X = runif(n, -1, 1)
  W = as.numeric(X >= threshold)
  Y = (1 + 2*W)*(1 + X^2) + 1 / (1 + exp(X)) + rnorm(n, sd = .5)
  Y = Y / diff(range(Y))
  for (dc in c(FALSE, TRUE)) {
    args = list(Y = Y, X = X, threshold = threshold, Lipschitz.constant = 0.5,
                sigma.sq = 0.01, diff.curvatures = dc, num.bucket = 100)
    q = do.call(plrd:::plrd.optim, c(args, use.spline = FALSE, optimizer = "quadprog"))
    cv = do.call(plrd:::plrd.optim, c(args, optimizer = "cvxr"))
    expect_equal(cv$tau.hat, q$tau.hat, tolerance = 1e-3)
    expect_equal(cv$max.bias, q$max.bias, tolerance = 1e-3)
    expect_equal(cv$gamma, q$gamma, tolerance = 1e-2)
    expect_equal(cv$gamma.fun.0, q$gamma.fun.0, tolerance = 1e-2)
    expect_equal(cv$gamma.fun.1, q$gamma.fun.1, tolerance = 1e-2)
  }
})

test_that("cvxr weights satisfy the moment constraints", {
  set.seed(7)
  n = 1000; threshold = 0
  X = runif(n, -1, 1)
  W = as.numeric(X >= threshold)
  Y = (1 + 2*W)*(1 + X^2) + 1 / (1 + exp(X)) + rnorm(n, sd = .5)
  out = plrd(Y, X, threshold)
  # The moment constraints in X hold at the bin centers, not at the raw X,
  # so only the normalizations are exact at the observation level.
  g = out$gamma
  expect_equal(sum(g[W == 1]), 1, tolerance = 1e-6)
  expect_equal(sum(g[W == 0]), -1, tolerance = 1e-6)
})
