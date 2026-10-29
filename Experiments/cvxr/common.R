## Shared helpers for the scripts in Experiments/cvxr/. Run every script from the
## repository root (Rscript Experiments/cvxr/<script>.R); they source the package's own
## R/ files, so they test this checkout without installing it.
##
## Every solution is scored by the same yardstick: the problem plrd.optim() itself sets
## up (its grid, its third-difference operator, its centering conditions). No solver's
## self-report is trusted.

if (!file.exists("DESCRIPTION") || read.dcf("DESCRIPTION", "Package")[1] != "plrd")
  stop("Run from the root of the plrd repository: Rscript Experiments/cvxr/<script>.R")
for (f in list.files("R", full.names = TRUE)) source(f)
data_dir <- file.path("Experiments", "Datasets")

## Output directory for scripts that save results: CVXR_EXPT_OUT, else a temporary one.
out_dir <- function() {
  d <- Sys.getenv("CVXR_EXPT_OUT", file.path(tempdir(), "plrd_cvxr_expts"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

versions <- function() {
  cat("R", format(getRversion()), "| CVXR", format(packageVersion("CVXR")),
      "| clarabel", format(packageVersion("clarabel")), "\n")
}

## The bias operator on plrd.optim()'s grid: h^3 * (M^{-T})[-(1:3), ], M = rbind(C, D3).
bias_operator <- function(N, h, zeroed.idx) {
  D3 <- Matrix::bandSparse(n = N - 3, m = N, k = 0:3,
          diag = list(rep(1, N), rep(-3, N), rep(3, N - 1), rep(-1, N - 2)))
  C <- Matrix::sparseMatrix(dims = c(3, N), i = 1:3, j = zeroed.idx, x = 1)
  h^3 * solve(t(as.matrix(rbind(C, D3))))[-(1:3), , drop = FALSE]
}

## plrd.optim()'s grid, rebuilt from its arguments exactly as plrd.optim() builds it.
optim_grid <- function(X, threshold, num.bucket) {
  Xm <- matrix(X, ncol = 1)
  h <- (max(Xm) - min(Xm)) / num.bucket
  breaks <- seq(min(Xm) - h/2, max(Xm) + h, by = h)
  xx <- breaks[-1] - h/2; N <- length(xx); xc <- xx - threshold
  zi <- max(which(xc < 0)) + c(0, 1, 2)
  list(N = N, xc = xc, bin = as.numeric(cut(Xm[, 1], breaks = breaks)),
       A = bias_operator(N, h, zi))
}

## Score observation-level weights gamma (as plrd.optim() returns them) against the
## problem for these arguments: the TRUE worst-case bias B * ||A v||_1 (v = weight per
## bin; one norm per side when the curvatures differ), the largest violation of the
## constraints at the bin centers, and the fold objective sigma.sq * sum(gamma^2) + bias^2.
score_weights <- function(a, gamma, num.bucket = if (is.null(a$num.bucket)) 400 else a$num.bucket) {
  g <- optim_grid(a$X, a$threshold, num.bucket)
  W <- as.numeric(a$X >= a$threshold)
  per_bin <- function(x) { v <- as.numeric(tapply(x, factor(g$bin, levels = 1:g$N), sum)); v[is.na(v)] <- 0; v }
  v0 <- per_bin(gamma * (W == 0)); v1 <- per_bin(gamma * (W == 1))
  dc <- isTRUE(a$diff.curvatures)
  l1 <- if (dc) sum(abs(g$A %*% v0)) + sum(abs(g$A %*% v1)) else sum(abs(g$A %*% (v0 + v1)))
  xb <- g$xc[g$bin]
  cons <- c(sum(gamma[W == 1]) - 1, sum(gamma[W == 0]) + 1,
            sum(gamma * xb * W), sum(gamma * xb * (1 - W)),
            if (dc) c(sum(gamma * xb^2 * W), sum(gamma * xb^2 * (1 - W))) else sum(gamma * xb^2))
  bias <- a$Lipschitz.constant * l1
  c(bias = bias, cons = max(abs(cons)), obj = a$sigma.sq * sum(gamma^2) + bias^2)
}

## Score per-bin weights (as plrd.primal.cvxr() returns them) against the same problem.
score_bins <- function(p, sol) {
  A <- bias_operator(length(p$xx.centered), p$bin.width, p$zeroed.idx)
  xc <- as.numeric(p$xx.centered)
  v0 <- p$counts.0 * sol$gamma.0; v1 <- p$counts.1 * sol$gamma.1
  dc <- isTRUE(p$diff.curvatures)
  l1 <- if (dc) sum(abs(A %*% v0)) + sum(abs(A %*% v1)) else sum(abs(A %*% (v0 + v1)))
  cons <- c(sum(v1) - 1, sum(v0) + 1, sum(v1 * xc), sum(v0 * xc),
            if (dc) c(sum(v1 * xc^2), sum(v0 * xc^2)) else sum((v0 + v1) * xc^2))
  bias <- p$Lipschitz.constant * l1
  c(bias = bias, reported = p$Lipschitz.constant * sol$t.hat, cons = max(abs(cons)),
    obj = p$sigma.sq * (sum(p$counts.0 * sol$gamma.0^2) + sum(p$counts.1 * sol$gamma.1^2)) + bias^2)
}

## plrd.primal.cvxr() with the solver and solver settings as arguments; the body is
## otherwise the package's own. defaults = TRUE drops the package's tol_gap_rel, i.e.
## uses CLARABEL's default tolerances.
primal_with_solver <- function(defaults = FALSE) {
  edit <- function(e) {
    if (is.character(e) && length(e) == 1 && identical(e, "CLARABEL")) return(quote(solver))
    if (is.call(e)) {
      if (identical(e[[1]], quote(CVXR::psolve))) {
        if (defaults) e$tol_gap_rel <- NULL
        e <- as.call(c(as.list(e), quote(...)))
      }
      return(as.call(lapply(as.list(e), edit)))
    }
    e
  }
  f <- plrd.primal.cvxr
  body(f) <- edit(body(f))
  formals(f) <- c(formals(f), alist(solver = "CLARABEL", ... = ))
  f
}

## The four empirical datasets of Table 2, prepared as in Experiments/Datasets/
## Results_for_real_datasets.R, plus the README example.
table2_datasets <- function() {
  d <- list()
  sen <- read.csv(file.path(data_dir, "senate.csv"))
  d$`Cattaneo et al. (2015)` <- list(Y = sen$y, X = sen$x, c = 0, paper = "5.830 ± 2.127")
  hs <- read.csv(file.path(data_dir, "headstart.csv"))
  ok <- !is.na(hs$mort_age59_related_postHS) & !is.na(hs$povrate60)
  d$`Ludwig & Miller (2007)` <- list(Y = hs$mort_age59_related_postHS[ok], X = hs$povrate60[ok],
                                     c = 59.1984, paper = "-1.249 ± 1.507")
  mm <- read.csv(file.path(data_dir, "mats_math.csv"))
  d$`Matsudaira (2008) math` <- list(Y = mm$y, X = -mm$x, c = 0, paper = "0.073 ± 0.028")
  jl <- read.csv(file.path(data_dir, "jl_math.csv"))
  d$`Jacob & Lefgren (2004) math` <- list(Y = jl$y, X = -jl$x, c = 0, paper = "0.168 ± 0.039")
  d$`README example` <- c(readme_example(), paper = "-")
  d
}

## The README / ?plrd example (n = 1000, X ~ U(-1, 1), threshold 0, seed 42).
readme_example <- function() {
  set.seed(42)
  X <- runif(1000, -1, 1); W <- as.numeric(X >= 0)
  list(Y = (1 + 2*W)*(1 + X^2) + 1 / (1 + exp(X)) + rnorm(1000, sd = .5), X = X, c = 0)
}
