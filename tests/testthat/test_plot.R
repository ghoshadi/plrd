test_that("plotting runs across settings", {

  n <- 10000
  c <- 2
  set.seed(42)

  # setting 1: uniform X, shallow sigmoid at c
  X = runif(n, 0, 4); W = as.numeric(X >= c)
  Y = 0.5 - 0.2 * X + 1.5 / (1 + exp(-(X - c) / 0.3)) + 1.2 * W + rnorm(n, sd = 0.3)
  fit1 = plrd(Y, X, c)
  plot(fit1, type = "default")
  plot(fit1, type = "weights")
  plot(fit1, type = "combined")

  # setting 2: normal X, gentle slope + curvature and positive effect
  X = rnorm(n, mean = c, sd = 1); W = as.numeric(X >= c)
  Y = 0.5 - 0.25 * X - 0.04 * X^2 + 1.5 * W + rnorm(n, sd = 0.35)
  fit2 = plrd(Y, X, c)
  plot(fit2, type = "default")
  plot(fit2, type = "weights")
  plot(fit2, type = "combined")

  # setting 3: normal X, more  noise, threshold in tail of distribution
  X = rnorm(n, mean = c, sd = 1); W = as.numeric(X >= 4)
  Y = 1 - 0.4 * X^2 + 2.0 * W + rnorm(n, sd = 0.9)
  fit3 = plrd(Y, X, 4)
  plot(fit3, type = "default")
  plot(fit3, type = "weights")
  plot(fit3, type = "combined")

  # setting 4: normal X, small jump, curved baseline
  X = rnorm(n, mean = c, sd = 1); W = as.numeric(X >= 1)
  Y = 0.5 - 0.1 * X + 0.3 * X^2 - 0.08 * X^3 + 2.0 * W + rnorm(n, sd = 0.3)
  fit4 = plrd(Y, X, 1)
  plot(fit4, type = "default")
  plot(fit4, type = "weights")
  plot(fit4, type = "combined")

  # setting 5: normal X, different curvature above vs. below c
  X = rnorm(n, mean = c, sd = 1); W = as.numeric(X >= c)
  Y = 0.5 - 0.10 * X + 0.15 * X^2 * (1 - W) - 0.08 * X^2 * W + 2 * W + rnorm(n, sd = 0.3)
  fit5 = plrd(Y, X, c)
  plot(fit5, type = "default")
  plot(fit5, type = "weights")
  plot(fit5, type = "combined")

  # setting 6: normal X, different curvature above vs. below c, negative effect
  X = rnorm(n, mean = c, sd = 1); W = as.numeric(X >= 3)
  Y = 0.6 * X^2 -2  * W + rnorm(n, sd = 0.8) + (X >= 2) * (-0.3 * abs(X - c) * X^2)
  fit6 = plrd(Y, X, 3)
  plot(fit6, type = "default")
  plot(fit6, type = "weights")
  plot(fit6, type = "combined")

  # setting 7: lower noise
  X = rnorm(n, mean = c, sd = 1); W = as.numeric(X >= 3)
  Y = 0.6 * X^2 + 1.5 * W + rnorm(n, sd = 0.5) + (X >= 2) * (-0.4 * abs(X - c) * X^2)
  fit7 = plrd(Y, X, 3)
  plot(fit7, type = "default")
  plot(fit7, type = "weights")
  plot(fit7, type = "combined")

  expect_true(TRUE)
})
