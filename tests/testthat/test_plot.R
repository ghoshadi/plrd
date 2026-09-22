test_that("plot plrd rejects invalid plot types and effective support windows", {
  set.seed(42)
  n <- 5000
  threshold <- 2
  X <- rnorm(n, mean = threshold, sd = 1)
  Y <- 0.5 - 0.2 * X + 1.2 * as.numeric(X >= threshold) + rnorm(n, sd = 0.3)
  fit <- plrd(Y, X, threshold)

  # Test whether invalid test inputs are rejected
  expect_error(find_weight_window(fit, percentage.cumulative.weights = 1.5))
  expect_error(find_weight_window(fit, percentage.cumulative.weights = 0))
  expect_error(find_weight_window(fit, percentage.cumulative.weights = -0.2))

  # Test whether computing weight window errors on valid cumulative percentage
  expect_no_error(find_weight_window(fit, percentage.cumulative.weights = 0.99))
  expect_no_error(find_weight_window(fit, percentage.cumulative.weights = 1)) # Shows all weights

  expect_error(plot(fit, type = "incompatible.plot"), "select plot type")

  # Test extraction of coordinates for main and weight curves
  plotting.coordinates <- plot(fit, type = "combined")
  expect_no_error({plot(plotting.coordinates$main$x_coordinates, plotting.coordinates$main$y_coordinates,
       xlab = "Running Variable (X)",
       ylab = "Response (Y)")
       abline(v = plotting.coordinates$gamma$window_weights, lty = 2, col = "darkgreen", lwd = 2)})
  expect_equal(length(plotting.coordinates$main$x_coordinates),
               length(plotting.coordinates$main$y_coordinates))
  expect_no_error(plot(plotting.coordinates$gamma$x_coordinates, plotting.coordinates$gamma$y_coordinates,
                       xlab = "Running Variable (X)",
                       ylab = "Weight"
                       ))
  expect_equal(length(plotting.coordinates$gamma$x_coordinates),
               length(plotting.coordinates$gamma$y_coordinates))
})


test_that("plotting runs across settings and types without error", {
  set.seed(42)
  n <- 10000
  threshold <- 2

  # setting 1: uniform X, shallow sigmoid at threshold
  X <- runif(n, 0, 4); W <- as.numeric(X >= threshold)
  Y <- 0.5 - 0.2 * X + 1.5 / (1 + exp(-(X - threshold) / 0.3)) + 1.2 * W + rnorm(n, sd = 0.3)
  fit1 <- plrd(Y, X, threshold)
  for (t in c("default", "weights", "combined")) expect_no_error(plot(fit1, type = t))

  # setting 2: normal X, gentle slope + curvature and positive effect
  X <- rnorm(n, mean = threshold, sd = 1); W <- as.numeric(X >= threshold)
  Y <- 0.5 - 0.25 * X - 0.04 * X^2 + 1.5 * W + rnorm(n, sd = 0.35)
  fit2 <- plrd(Y, X, threshold)
  for (t in c("default", "weights", "combined")) expect_no_error(plot(fit2, type = t))

  # setting 3: normal X, more noise, threshold in tail of distribution
  X <- rnorm(n, mean = threshold, sd = 1); W <- as.numeric(X >= 4)
  Y <- 1 - 0.4 * X^2 + 2.0 * W + rnorm(n, sd = 0.9)
  fit3 <- plrd(Y, X, 4)
  for (t in c("default", "weights", "combined")) expect_no_error(plot(fit3, type = t))

  # setting 4: normal X, small jump, curved baseline
  X <- rnorm(n, mean = threshold, sd = 1); W <- as.numeric(X >= 1)
  Y <- 0.5 - 0.1 * X + 0.3 * X^2 - 0.08 * X^3 + 2.0 * W + rnorm(n, sd = 0.3)
  fit4 <- plrd(Y, X, 1)
  for (t in c("default", "weights", "combined")) expect_no_error(plot(fit4, type = t))

  # setting 5: normal X, different curvature above vs. below threshold
  X <- rnorm(n, mean = threshold, sd = 1); W <- as.numeric(X >= threshold)
  Y <- 0.5 - 0.10 * X + 0.15 * X^2 * (1 - W) - 0.08 * X^2 * W + 2 * W + rnorm(n, sd = 0.3)
  fit5 <- plrd(Y, X, threshold)
  for (t in c("default", "weights", "combined")) expect_no_error(plot(fit5, type = t))

  # setting 6: normal X, different curvature above vs. below threshold, negative effect
  X <- rnorm(n, mean = threshold, sd = 1); W <- as.numeric(X >= 3)
  Y <- 0.6 * X^2 - 2 * W + rnorm(n, sd = 0.8) + (X >= 2) * (-0.3 * abs(X - threshold) * X^2)
  fit6 <- plrd(Y, X, 3)
  for (t in c("default", "weights", "combined")) expect_no_error(plot(fit6, type = t))

  # setting 7: lower noise
  X <- rnorm(n, mean = threshold, sd = 1); W <- as.numeric(X >= 3)
  Y <- 0.6 * X^2 + 1.5 * W + rnorm(n, sd = 0.5) + (X >= 2) * (-0.4 * abs(X - threshold) * X^2)
  fit7 <- plrd(Y, X, 3)
  for (t in c("default", "weights", "combined")) expect_no_error(plot(fit7, type = t))
})
