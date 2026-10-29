## CLARABEL (the package's solver) vs PIQP, and the effect of the rescaling in R/cvxr.R.
## Captures the 120 fold problems of 60 plrd() runs (3 designs x seeds 1-20 x 2 folds:
## X ~ U(-3, 3) with threshold 1.42; X ~ U(-1, 1) with threshold 0; the first with X
## rounded to 0.01), solves each with both solvers in the scaled formulation (the
## package's) and unscaled, scores every solution against plrd.optim()'s own problem
## (common.R), and times 6 of them with bench::mark. Each solver runs in its own R
## process so the timings do not share a session. Needs packages piqp and bench.
## About 5 minutes.
##
##   Rscript Experiments/cvxr/solver_comparison.R
source(file.path("Experiments", "cvxr", "common.R"))
suppressPackageStartupMessages(library(CVXR))
solvers <- c("CLARABEL", "PIQP")
solver <- Sys.getenv("SOLVER")

if (solver == "") {
  ## Driver: one worker process per solver, then the summary.
  out <- out_dir()
  for (s in solvers) {
    status <- system2(file.path(R.home("bin"), "Rscript"), "Experiments/cvxr/solver_comparison.R",
                      env = c(paste0("SOLVER=", s), paste0("CVXR_EXPT_OUT=", out)))
    if (status != 0) stop("worker for ", s, " failed")
  }
  versions(); cat("piqp", format(packageVersion("piqp")), "\n")
  acc <- do.call(rbind, lapply(solvers, function(s) readRDS(file.path(out, paste0("acc_", s, ".rds")))))
  best <- tapply(acc$obj, acc$id, min, na.rm = TRUE)
  acc$gap <- (acc$obj - best[as.character(acc$id)]) / best[as.character(acc$id)]
  cat("\n== Accuracy on 120 fold problems; gap = relative excess over the best true objective found ==\n")
  for (sc in c(FALSE, TRUE)) for (des in c(unique(acc$design), "ALL")) {
    cat(sprintf("\n[%s | %s]\n", if (sc) "scaled (the package)" else "unscaled", des))
    for (s in solvers) {
      a <- acc[acc$solver == s & acc$scaled == sc & (des == "ALL" | acc$design == des), ]
      cat(sprintf("  %-8s optimal %3d  inaccurate %3d  other %2d | gap median %.1e max %.1e | max cons %.1e | max |bias error| %.1e\n",
                  s, sum(a$status == "optimal"), sum(a$status == "optimal_inaccurate"),
                  sum(!a$status %in% c("optimal", "optimal_inaccurate")),
                  median(a$gap, na.rm = TRUE), max(a$gap, na.rm = TRUE),
                  max(a$cons, na.rm = TRUE), max(abs(a$bias_err), na.rm = TRUE)))
    }
  }
  tm <- do.call(rbind, lapply(solvers, function(s) readRDS(file.path(out, paste0("time_", s, ".rds")))))
  tm$min_ms <- tm$min * 1e3; tm$median_ms <- tm$median * 1e3
  cat("\n== Time per fold solve (bench::mark, 10 iterations; includes building the operator and the model) ==\n")
  print(aggregate(cbind(min_ms, median_ms) ~ solver + scaled, tm, median), digits = 3, row.names = FALSE)
  cat("(median over the 6 timed problems)\n")
  quit(save = "no")
}

## ---- Worker: one solver ----
problems <- list()
inner <- plrd.primal.cvxr
plrd.primal.cvxr <- function(...) {
  problems[[length(problems) + 1]] <<- c(list(...), design = design, seed = s)
  suppressWarnings(inner(...))
}
designs <- list(
  shifted = function() { X <- runif(1000, -3, 3); list(X = X, c = 1.42) },
  centered = function() { X <- runif(1000, -1, 1); list(X = X, c = 0) },
  shifted_discrete = function() { X <- round(runif(1000, -3, 3), 2); list(X = X, c = 1.42) })
for (design in names(designs)) for (s in 1:20) {
  set.seed(s); d <- designs[[design]]()
  W <- as.numeric(d$X >= d$c)
  Y <- (1 + 2*W)*(1 + d$X^2) + 1 / (1 + exp(d$X)) + rnorm(1000, sd = .5)
  invisible(suppressMessages(plrd(Y, d$X, d$c)))
}
plrd.primal.cvxr <- inner

## R/cvxr.R's program with the solver and the scaling as switches. CLARABEL gets the
## package's tol_gap_rel = 1e-7; PIQP its defaults.
primal <- function(p, solver, scaled) {
  N <- length(p$xx.centered); xc <- as.numeric(p$xx.centered)
  Aop <- bias_operator(N, p$bin.width, p$zeroed.idx)
  r0 <- p$realized.idx.0; r1 <- p$realized.idx.1
  w0 <- p$counts.0[r0]; w1 <- p$counts.1[r1]; x0 <- xc[r0]; x1 <- xc[r1]
  A0 <- Aop[, r0, drop = FALSE] %*% diag(w0, length(w0))
  A1 <- Aop[, r1, drop = FALSE] %*% diag(w1, length(w1))
  B <- p$Lipschitz.constant
  n.obs <- if (scaled) sum(w0) + sum(w1) else 1
  op.scale <- if (scaled) 1 / max(abs(A0), abs(A1)) else 1
  A0 <- op.scale * A0; A1 <- op.scale * A1
  u0 <- Variable(length(w0)); u1 <- Variable(length(w1)); t <- Variable(1)
  bias <- if (p$diff.curvatures) p_norm(A0 %*% u0, 1) + p_norm(A1 %*% u1, 1)
          else p_norm(A0 %*% u0 + A1 %*% u1, 1)
  cons <- list(sum(w1 * u1) == n.obs, sum(w0 * u0) == -n.obs,
               sum(w1 * x1 * u1) == 0, sum(w0 * x0 * u0) == 0,
               (if (scaled) B / op.scale else 1) * bias <= t)
  cons <- c(cons, if (p$diff.curvatures) list(sum(w1 * x1^2 * u1) == 0, sum(w0 * x0^2 * u0) == 0)
                  else list(sum(w1 * x1^2 * u1) + sum(w0 * x0^2 * u0) == 0))
  obj <- p$sigma.sq * (sum_squares(sqrt(w0) * u0) + sum_squares(sqrt(w1) * u1)) +
         (if (scaled) 1 else B^2) * square(t)
  prob <- Problem(Minimize(obj), cons)
  settings <- if (solver == "CLARABEL") list(tol_gap_rel = 1e-7) else list()
  tryCatch(suppressWarnings(do.call(psolve, c(list(prob, solver = solver), settings))),
           error = function(e) NULL)
  if (is.null(value(u0))) return(list(status = status(prob)))
  g0 <- rep(0, N); g1 <- rep(0, N)
  g0[r0] <- as.numeric(value(u0)) / n.obs; g1[r1] <- as.numeric(value(u1)) / n.obs
  list(status = status(prob), gamma.0 = g0, gamma.1 = g1,
       t.hat = as.numeric(value(t)) / (if (scaled) B * n.obs else 1))
}

## Faithfulness: the scaled CLARABEL formulation here must reproduce the package exactly.
if (solver == "CLARABEL") {
  p1 <- problems[[1]]
  pk <- do.call(plrd.primal.cvxr, p1[setdiff(names(p1), c("design", "seed"))])
  mr <- primal(p1, "CLARABEL", TRUE)
  stopifnot(identical(pk$gamma.0, mr$gamma.0), identical(pk$gamma.1, mr$gamma.1), identical(pk$t.hat, mr$t.hat))
}

rows <- list()
for (i in seq_along(problems)) for (sc in c(FALSE, TRUE)) {
  p <- problems[[i]]; r <- primal(p, solver, sc)
  sco <- if (is.null(r$gamma.0)) c(obj = NA, cons = NA, bias = NA, reported = NA) else score_bins(p, r)
  rows[[length(rows) + 1]] <- data.frame(id = i, design = p$design, seed = p$seed, solver = solver,
    scaled = sc, status = r$status, obj = sco[["obj"]], cons = sco[["cons"]],
    bias_err = sco[["reported"]] / sco[["bias"]] - 1)
}
saveRDS(do.call(rbind, rows), file.path(out_dir(), paste0("acc_", solver, ".rds")))

## Timing: the first fold of seeds 1 and 2 in each design, scaled and unscaled.
pick <- which(!duplicated(paste(sapply(problems, `[[`, "design"), sapply(problems, `[[`, "seed"))) &
              sapply(problems, `[[`, "seed") %in% 1:2)
invisible(primal(problems[[1]], solver, TRUE))   # session warm-up, outside the timings
tm <- list()
for (i in pick) for (sc in c(FALSE, TRUE)) {
  b <- bench::mark(primal(problems[[i]], solver, sc), iterations = 10, check = FALSE, filter_gc = FALSE)
  tm[[length(tm) + 1]] <- data.frame(id = i, design = problems[[i]]$design, solver = solver, scaled = sc,
                                     min = as.numeric(b$min), median = as.numeric(b$median))
}
saveRDS(do.call(rbind, tm), file.path(out_dir(), paste0("time_", solver, ".rds")))
