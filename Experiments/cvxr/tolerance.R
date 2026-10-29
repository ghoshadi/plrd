## Why R/cvxr.R sets CLARABEL's tol_gap_rel = 1e-7 (default 1e-8), and that it is set
## only to avoid a spurious "Solution may be inaccurate" warning. Runs plot_settings.R
## under CLARABEL's default tolerances, takes the folds that come back
## optimal_inaccurate, prints each solve's iteration log (it stalls, step 0, at a relative
## gap of 2-3e-8 with the residuals already below 1e-8), then re-solves under candidate
## settings: only tol_gap_rel clears the status, and the true objective is the same to
## 13 digits under every setting.
##
##   Rscript Experiments/cvxr/tolerance.R
Sys.setenv(CLARABEL_DEFAULTS = "1")
source(file.path("Experiments", "cvxr", "plot_settings.R"))
print(clarabel::clarabel_control()[c("max_iter", "tol_gap_abs", "tol_gap_rel", "tol_feas", "tol_ktratio")])

flagged <- Filter(function(a) !is.null(a$warned), captured)
cat("\nfolds reported optimal_inaccurate:", length(flagged), "\n")
for (k in seq_along(flagged)) {
  cat(sprintf("\n######## flagged fold %d: CLARABEL log, default settings ########\n", k))
  invisible(suppressWarnings(do.call(primal, modifyList(problem_args(flagged[[k]]), list(verbose = TRUE)))))
}

settings <- list(default = list(), max_iter_500 = list(max_iter = 500L),
                 tol_gap_rel_5e8 = list(tol_gap_rel = 5e-8), tol_gap_rel_1e7 = list(tol_gap_rel = 1e-7),
                 tol_feas_1e7 = list(tol_feas = 1e-7), tol_ktratio_1e5 = list(tol_ktratio = 1e-5))
cat("\n")
for (k in seq_along(flagged)) for (nm in names(settings)) {
  st <- "optimal"
  r <- withCallingHandlers(do.call(primal, c(problem_args(flagged[[k]]), settings[[nm]])),
         warning = function(w) { st <<- "inaccurate"; invokeRestart("muffleWarning") })
  cat(sprintf("fold %d  %-16s %-10s true objective %.12e\n", k, nm, st, score_bins(flagged[[k]], r)[["obj"]]))
}
