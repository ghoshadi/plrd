## How much does the optimizer matter, against the method's own sampling noise? On the
## senate data (Cattaneo et al. 2015), plrd() is run with seeds 1-20 (the seed sets the
## cross-fitting folds, after plrd()'s own fold-balancing loop), under the shipped default
## and under optimizer = "cvxr". Reported: the sd of tau-hat across seeds, the mean
## half-width, and the per-seed optimizer effect |tau(cvxr) - tau(spline)|. About 10 s.
##
##   Rscript Experiments/cvxr/fold_noise.R
source(file.path("Experiments", "cvxr", "common.R"))
versions()

sen <- read.csv(file.path(data_dir, "senate.csv"))
seeds <- 1:20
res <- do.call(rbind, lapply(seeds, function(s) {
  sp <- suppressMessages(plrd(sen$y, sen$x, 0, seed = s, optimizer = "quadprog"))
  cv <- suppressMessages(plrd(sen$y, sen$x, 0, seed = s, optimizer = "cvxr"))
  data.frame(seed = s, tau_spline = sp$tau.hat, tau_cvxr = cv$tau.hat,
             hw_spline = sp$tau.plusminus, hw_cvxr = cv$tau.plusminus)
}))
print(res, digits = 5, row.names = FALSE)

cat(sprintf("\nacross %d seeds:\n", length(seeds)))
cat(sprintf("  sd(tau-hat):       spline %.4f   cvxr %.4f\n", sd(res$tau_spline), sd(res$tau_cvxr)))
cat(sprintf("  mean half-width:   spline %.4f   cvxr %.4f\n", mean(res$hw_spline), mean(res$hw_cvxr)))
eff <- abs(res$tau_cvxr - res$tau_spline)
cat(sprintf("  optimizer effect |tau(cvxr) - tau(spline)|: median %.4f, max %.4f\n", median(eff), max(eff)))
d42 <- with(list(sp = suppressMessages(plrd(sen$y, sen$x, 0, optimizer = "quadprog")),
                 cv = suppressMessages(plrd(sen$y, sen$x, 0, optimizer = "cvxr"))), cv$tau.hat - sp$tau.hat)
cat(sprintf("  at the default seed (42): %.4f\n", d42))
