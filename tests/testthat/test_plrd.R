
test_that("plrd basic example runs", {
  library(plrd)
  # Simple example of regression discontinuity design 
  set.seed(42)
  n = 1000; threshold = 0
  X = runif(n, -1, 1)
  W = as.numeric(X >= threshold)
  Y = (1 + 2*W)*(1 + X^2) + 1 / (1 + exp(X)) + rnorm(n, sd = .5)
  out = plrd(Y, X, threshold)
  print(out)
  summary(out)
  plot(out)
})


test_that("plrd has not changed", {
  set.seed(42)
  n = 10000; threshold = 1.42
  X = runif(n, -3, 3)
  W = as.numeric(X >= threshold)
  Y = (1 + 2*W)*(1 + X^2) + 1 / (1 + exp(X)) + rnorm(n, sd = .5)
  out = plrd(Y, X, threshold)
  expect_equal(out$tau.hat, 6.0581497597, tol = 1e-6)
  expect_equal(out$sampling.se, 0.04943550781, tol = 1e-6)
})

test_that("plrd works as expected under translation and scaling", {
  n = 1000; threshold = 1.42
  X = runif(n, -3, 3)
  W = as.numeric(X >= threshold)
  Y = (1 + 2*W)*(1 + X^2) + 1 / (1 + exp(X)) + rnorm(n, sd = .5)

  out = plrd(Y, X, threshold)
  out.flipped.Y = plrd(-Y, X, threshold)
  out.flipped.X = plrd(Y, -X, -threshold)
  out.flipped.both = plrd(-Y, -X, -threshold)
  out.translated.Y = plrd(Y + 42, X, threshold)
  out.translated.X = plrd(Y, X + 42, threshold = threshold + 42)
  out.scaled.Y = plrd(10*Y, X, threshold)

  expect_equal(out$tau.hat, -1 * out.flipped.X$tau.hat, tolerance =  1e-6)
  expect_equal(out$tau.hat, -1 * out.flipped.Y$tau.hat, tolerance =  1e-6)
  expect_equal(out$tau.hat, out.flipped.both$tau.hat, tolerance =  1e-6)
  expect_equal(out$tau.hat, out.translated.Y$tau.hat, tolerance =  1e-6)
  expect_equal(out$tau.hat, out.translated.X$tau.hat, tolerance =  1e-6)
  expect_equal(out$tau.hat, out.scaled.Y$tau.hat/10, tolerance =  1e-6)
})

test_that("Treatment effect flips if running variable is reversed", {
  n = 1000; threshold = 0
  X = runif(n, -1, 1)
  W = as.numeric(X >= threshold)
  Y = (1 + 2*W)*(1 + X^2) + 1 / (1 + exp(X)) + rnorm(n, sd = .5)

  out = plrd(Y, X, threshold)
  out.flipped = plrd(Y, -X, threshold)

  expect_equal(out$tau.hat, -1 * out.flipped$tau.hat, tolerance =  1e-6)
})

test_that("Treatment effect flips if discrete running variable is reversed", {
  n = 1000; threshold = 0
  X = runif(n, -1, 1)
  # Discretize running var
  X = round(X,2)
  W = as.numeric(X >= threshold)
  Y = (1 + 2*W)*(1 + X^2) + 1 / (1 + exp(X)) + rnorm(n, sd = .5)
  
  out = plrd(Y, X, threshold)
  # Take care of units exactly at threshold, should be not treated when reversed
  X.rev = X
  X.rev[X.rev == threshold] = .Machine$double.eps
  X.rev = -X.rev
  out.flipped = plrd(Y, X.rev, threshold)
  expect_equal(out$tau.hat, -1 * out.flipped$tau.hat, tolerance = 1e-6)
})
