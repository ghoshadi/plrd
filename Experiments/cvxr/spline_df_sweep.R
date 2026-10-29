## Issue 1: the error in the reported worst-case bias as a function of spline.df, on the
## README example, at plrd.optim()'s default grid (num.bucket = 400). The inputs are
## scaled and B, sigma.sq estimated exactly as plrd() does before it calls plrd.optim().
## The error is relative to the exact problem's optimal worst-case bias, from
## use.spline = FALSE (quadprog) and from optimizer = "cvxr"; on this centered design
## the two agree. It is non-monotone in spline.df and changes sign.
##
##   Rscript Experiments/cvxr/spline_df_sweep.R
source(file.path("Experiments", "cvxr", "common.R"))
versions()

d <- readme_example()
th <- d$c; W <- as.numeric(d$X >= th)
scale.Y <- diff(range(d$Y)); Y.s <- d$Y / scale.Y
mw <- max(abs(d$X - th)); X.s <- th + (d$X - th) / mw
args <- list(Y = Y.s, X = X.s, threshold = th,
             Lipschitz.constant = get.Lipschitz.constant(Y.s, X.s, th, FALSE),
             sigma.sq = summary(stats::lm(Y.s ~ X.s * W))$sigma^2,
             diff.curvatures = FALSE, num.bucket = 400)

ex <- do.call(plrd.optim, c(args, use.spline = FALSE, optimizer = "quadprog"))
cv <- do.call(plrd.optim, c(args, optimizer = "cvxr"))
cat(sprintf("%-22s max.bias %.10f  tau %.10f\n", "quadprog exact", ex$max.bias, ex$tau.hat))
cat(sprintf("%-22s max.bias %.10f  tau %.10f\n", "cvxr", cv$max.bias, cv$tau.hat))
cat("\nrelative error of the spline path's reported max.bias vs the exact optimum:\n")
for (df in c(10, 20, 40, 80, 160)) {
  s <- do.call(plrd.optim, c(args, use.spline = TRUE, spline.df = df, optimizer = "quadprog"))
  cat(sprintf("  spline.df = %-4d max.bias %.10f | vs quadprog exact %+7.3f%% | vs cvxr %+7.3f%%\n",
              df, s$max.bias, 100 * (s$max.bias / ex$max.bias - 1), 100 * (s$max.bias / cv$max.bias - 1)))
}
