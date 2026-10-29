## Table 2 end to end through plrd(), three runs per dataset differing ONLY in the
## optimizer:
##   spline = the shipped default (optimizer = "quadprog", use.spline = TRUE, spline.df = 40)
##   exact  = the exact dual   (optimizer = "quadprog", use.spline = FALSE)
##   cvxr   = optimizer = "cvxr"
## Folds, B-hat, sigma-hat and the standard error come from the same code in all three.
## Each fold's weights are scored against plrd.optim()'s own problem (common.R), giving
## the TRUE worst-case bias, the half-width it implies, and the fold objective.
## Behind the PR's Table 2 comparison, Issue 1's table, Issue 2's numbers on Matsudaira
## and Jacob & Lefgren, and the timings. About 1 minute (the exact path is slow).
##
##   Rscript Experiments/cvxr/table2.R
source(file.path("Experiments", "cvxr", "common.R"))
options(warn = 1)
versions()

## Run plrd() while recording every fold's plrd.optim() arguments and result.
run <- function(d, ...) {
  folds <- list()
  inner <- plrd.optim
  plrd.optim <<- function(...) { a <- list(...); r <- inner(...); folds[[length(folds) + 1]] <<- list(a = a, r = r); r }
  on.exit(plrd.optim <<- inner)
  t0 <- Sys.time()
  out <- tryCatch(suppressMessages(plrd(d$Y, d$X, d$c, ...)), error = function(e) e)
  secs <- as.numeric(Sys.time() - t0, units = "secs")
  if (inherits(out, "error")) return(list(error = conditionMessage(out)))
  sc <- sapply(folds, function(f) c(score_weights(f$a, f$r$gamma), reported = f$r$max.bias))
  ## plrd()'s max.bias is the fold average times scale.Y; the true bias scales the same way
  true_bias <- out$max.bias * sum(sc["bias", ]) / sum(sc["reported", ])
  list(out = out, true_bias = true_bias, cons = max(sc["cons", ]), secs = secs, obj = sc["obj", ],
       honest_pm = get.plusminus(true_bias, out$sampling.se, out$alpha))
}

datasets <- table2_datasets()
for (nm in names(datasets)) {
  d <- datasets[[nm]]
  cat(sprintf("\n=== %s (n = %d)   paper Table 2: %s ===\n", nm, length(d$Y), d$paper))
  runs <- list(spline = run(d, optimizer = "quadprog"),
               exact  = run(d, optimizer = "quadprog", use.spline = FALSE),
               cvxr   = run(d, optimizer = "cvxr"))
  for (k in names(runs)) {
    r <- runs[[k]]
    if (!is.null(r$error)) { cat(sprintf("  %-6s ERROR: %s\n", k, r$error)); next }
    o <- r$out
    cat(sprintf("  %-6s %9.4f ± %.4f | bias reported %.6g true %.6g (%+.2f%%) | honest ± %.4f (%+.2f%%) | max|cons| %.0e | %5.1fs\n",
                k, o$tau.hat, o$tau.plusminus, o$max.bias, r$true_bias,
                100 * (o$max.bias - r$true_bias) / r$true_bias, r$honest_pm,
                100 * (o$tau.plusminus - r$honest_pm) / r$honest_pm, r$cons, r$secs))
  }
  ok <- names(runs)[sapply(runs, function(r) is.null(r$error))]
  best <- apply(sapply(runs[ok], `[[`, "obj"), 1, min)
  cat("  fold objective, relative excess over the best of the three:\n")
  for (k in ok) cat(sprintf("    %-6s fold 1 %+.1e   fold 2 %+.1e\n", k,
                            (runs[[k]]$obj[1] - best[1]) / best[1], (runs[[k]]$obj[2] - best[2]) / best[2]))
}
cat("\n'true' = worst-case bias of the returned weights on plrd.optim()'s own grid;",
    "'honest ±' = the half-width that bias implies, same standard error.\n")
