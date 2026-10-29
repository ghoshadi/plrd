## Issue 2: the three plrd.optim() paths on the first cross-fitting fold of the design in
## tests/testthat/test_plrd.R (n = 1000, X ~ U(-3, 3), threshold 1.42), with B and
## sigma.sq as plrd() estimates them. Each path is scored only by the weights it returns
## (common.R): true objective, true worst-case bias against the reported one, and the
## largest constraint violation. Seeds 1-3 are the PR's table; the centered design
## (threshold 0, X ~ U(-1, 1)) is the control.
##
##   Rscript Experiments/cvxr/audit_paths.R            # seeds 1-3 off-center, 1-2 centered
source(file.path("Experiments", "cvxr", "common.R"))
versions()

audit_fold <- function(seed, c0, lo, hi) {
  captured <- NULL
  inner <- plrd.optim
  plrd.optim <<- function(...) { a <- list(...); if (is.null(captured)) captured <<- a; inner(...) }
  on.exit(plrd.optim <<- inner)
  set.seed(seed)
  X <- runif(1000, lo, hi); W <- as.numeric(X >= c0)
  Y <- (1 + 2*W)*(1 + X^2) + 1 / (1 + exp(X)) + rnorm(1000, sd = .5)
  invisible(suppressMessages(plrd(Y, X, c0)))
  plrd.optim <<- inner
  a <- captured
  cat(sprintf("\nseed %d, X ~ U(%g, %g), threshold %.2f | fold n = %d, B = %.3g, sigma.sq = %.3g, diff.curvatures = %s\n",
              seed, lo, hi, c0, length(a$Y), a$Lipschitz.constant, a$sigma.sq, a$diff.curvatures))
  paths <- list(`quadprog spline df=40` = list(spline.df = 40, optimizer = "quadprog"),
                `quadprog exact` = list(use.spline = FALSE, optimizer = "quadprog"),
                cvxr = list(optimizer = "cvxr"))
  for (nm in names(paths)) {
    r <- tryCatch(do.call(plrd.optim, c(list(a$Y, a$X, a$threshold, Lipschitz.constant = a$Lipschitz.constant,
                                             sigma.sq = a$sigma.sq, diff.curvatures = a$diff.curvatures,
                                             num.bucket = a$num.bucket), paths[[nm]])),
                  error = function(e) e)
    if (inherits(r, "error")) { cat(sprintf("  %-22s ERROR: %s\n", nm, conditionMessage(r))); next }
    s <- score_weights(a, r$gamma)
    cat(sprintf("  %-22s true obj %.6e | true bias %.6e (reported %.6e) | max|constraint| %.1e | tau %.6f\n",
                nm, s[["obj"]], s[["bias"]], r$max.bias, s[["cons"]], r$tau.hat))
  }
}

for (s in 1:3) audit_fold(s, 1.42, -3, 3)
for (s in 1:2) audit_fold(s, 0, -1, 1)
