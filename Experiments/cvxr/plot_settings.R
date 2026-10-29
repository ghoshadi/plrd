## The seven designs of tests/testthat/test_plot.R (same seed and draws), run through
## plrd() with optimizer = "cvxr". Every fold problem is captured and re-solved with PIQP
## and ECOS (packages piqp, ECOSolveR); each solution is scored against plrd.optim()'s own
## problem (common.R), relative to the best of the three. Shows that the package's
## CLARABEL solutions are optimal (gap ~1e-9) and that no fold reports "inaccurate".
## With CLARABEL_DEFAULTS=1 the package's tol_gap_rel = 1e-7 is dropped: 3 of the 14
## folds then report optimal_inaccurate, at the same solution (see tolerance.R).
##
##   Rscript Experiments/cvxr/plot_settings.R
##   CLARABEL_DEFAULTS=1 Rscript Experiments/cvxr/plot_settings.R
source(file.path("Experiments", "cvxr", "common.R"))
suppressMessages(library(CVXR))
versions()

primal <- primal_with_solver(defaults = Sys.getenv("CLARABEL_DEFAULTS") == "1")
captured <- list()
plrd.primal.cvxr <- function(...) {
  args <- list(...)
  out <- withCallingHandlers(primal(...), warning = function(w) {
    args$warned <<- conditionMessage(w); invokeRestart("muffleWarning")
  })
  args$out <- out
  captured[[length(captured) + 1]] <<- args
  out
}
problem_args <- function(a) a[setdiff(names(a), c("out", "warned"))]

n <- 10000
designs <- list(
  s1 = function() { X <- runif(n, 0, 4); W <- as.numeric(X >= 2)
    list(Y = 0.5 - 0.2 * X + 1.5 / (1 + exp(-(X - 2) / 0.3)) + 1.2 * W + rnorm(n, sd = 0.3), X = X, c = 2) },
  s2 = function() { X <- rnorm(n, 2, 1); W <- as.numeric(X >= 2)
    list(Y = 0.5 - 0.25 * X - 0.04 * X^2 + 1.5 * W + rnorm(n, sd = 0.35), X = X, c = 2) },
  s3 = function() { X <- rnorm(n, 2, 1); W <- as.numeric(X >= 4)
    list(Y = 1 - 0.4 * X^2 + 2.0 * W + rnorm(n, sd = 0.9), X = X, c = 4) },
  s4 = function() { X <- rnorm(n, 2, 1); W <- as.numeric(X >= 1)
    list(Y = 0.5 - 0.1 * X + 0.3 * X^2 - 0.08 * X^3 + 2.0 * W + rnorm(n, sd = 0.3), X = X, c = 1) },
  s5 = function() { X <- rnorm(n, 2, 1); W <- as.numeric(X >= 2)
    list(Y = 0.5 - 0.10 * X + 0.15 * X^2 * (1 - W) - 0.08 * X^2 * W + 2 * W + rnorm(n, sd = 0.3), X = X, c = 2) },
  s6 = function() { X <- rnorm(n, 2, 1); W <- as.numeric(X >= 3)
    list(Y = 0.6 * X^2 - 2 * W + rnorm(n, sd = 0.8) + (X >= 2) * (-0.3 * abs(X - 2) * X^2), X = X, c = 3) },
  s7 = function() { X <- rnorm(n, 2, 1); W <- as.numeric(X >= 3)
    list(Y = 0.6 * X^2 + 1.5 * W + rnorm(n, sd = 0.5) + (X >= 2) * (-0.4 * abs(X - 2) * X^2), X = X, c = 3) })

set.seed(42)    # one seed, designs drawn in order, as in test_plot.R
rows <- list()
for (nm in names(designs)) {
  d <- designs[[nm]]()
  before <- length(captured)
  invisible(suppressMessages(plrd(d$Y, d$X, d$c)))
  for (k in seq(before + 1, length(captured))) {
    a <- captured[[k]]
    s <- rbind(CLARABEL = score_bins(a, a$out),
               PIQP = score_bins(a, do.call(primal, c(problem_args(a), solver = "PIQP"))),
               ECOS = score_bins(a, do.call(primal, c(problem_args(a), solver = "ECOS"))))
    best <- min(s[, "obj"])
    rows[[length(rows) + 1]] <- data.frame(
      design = nm, fold = k - before, diff.curvatures = isTRUE(a$diff.curvatures),
      clarabel_status = if (is.null(a$warned)) "optimal" else "inaccurate",
      gap_clarabel = s["CLARABEL", "obj"] / best - 1, gap_piqp = s["PIQP", "obj"] / best - 1,
      gap_ecos = s["ECOS", "obj"] / best - 1,
      bias_err_clarabel = s["CLARABEL", "reported"] / s["CLARABEL", "bias"] - 1)
  }
}
res <- do.call(rbind, rows)
print(format(res, digits = 3), row.names = FALSE)
cat(sprintf("\nCLARABEL: %d of %d folds optimal_inaccurate; max gap %.1e; max |bias error| %.1e\n",
            sum(res$clarabel_status == "inaccurate"), nrow(res), max(res$gap_clarabel),
            max(abs(res$bias_err_clarabel))))
cat("gap = relative excess of the true objective over the best of CLARABEL, PIQP, ECOS\n")
