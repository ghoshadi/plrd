if (!requireNamespace("RDHonest", quietly = TRUE)) install.packages("RDHonest")
if (!requireNamespace("rdrobust", quietly = TRUE)) install.packages("rdrobust")
if (!requireNamespace("glmnet", quietly = TRUE)) install.packages("glmnet")
if (!requireNamespace("splines", quietly = TRUE)) install.packages("splines")
if (!requireNamespace("parallel", quietly = TRUE)) install.packages("parallel")
if (!requireNamespace("pbmcapply", quietly = TRUE)) install.packages("pbmcapply")
if (!requireNamespace("R.matlab", quietly = TRUE)) install.packages("R.matlab")
if (!requireNamespace("withr", quietly = TRUE)) install.packages("withr")

library(RDHonest)   # version 1.0.0
library(rdrobust)   # version 2.2
library(glmnet)     # version 4.1.9
library(splines)    # version 4.4.2
library(parallel)   # version 4.4.2
library(pbmcapply)  # version 1.5.1
library(R.matlab)   # version 3.7.0
library(plrd)       # our plrd package
library(withr)      # version 3.0.2

mred = "#CC3311"; mgreen = "#009E73"

# num_cores <- detectCores() - 1
# n = 500
# num_replications <- 1e4

#===========================================
#           Create your own
#===========================================

mu <- function(x){
  w = as.numeric(x>0)
  return( (1 + 2*w)*(1 +x^2) + 1 / (1 + exp(x)) )
}
truth <-  mu(1e-9) - mu(-1e-9)

single_experiment <- function(i) {
  set.seed(i)
  x = runif(n, -1, 1); y = sapply(x,  mu) + rnorm(n, sd = 0.5)

  plrd.out <- plrd(y, x, cutoff = 0) # PLRD estimates and CI's
  rdrob = ests.from.ci(rdrobust(y, x, c = 0, vce = 'hc1')$ci) # rdrobust estimates and CI's
  rdh = RDHonest(y ~ x, cutoff = 0)$coefficients # RDHonest estimates and CI's

  w1 = 2*rdrob[,2] # widths of rdrobust CI's (conventional, bias corrected, robust)
  w2 = rdh$conf.high-rdh$conf.low # width of RDHonest CI
  w3 <- 2*plrd.out$tau.plusminus # width of PLRD CI

  c1 = apply(rdrob, 1, function(x) as.numeric(abs(x[1]-truth) < x[2])) # coverages of rdrobust CI's (conventional, bias corrected, robust)
  c2 = as.numeric((rdh$conf.low<truth)*(rdh$conf.high>truth)) # coverage of RDHonest CI
  c3 <- as.numeric(abs(plrd.out$tau.hat - truth) <= w3/2) # coverage of PLRD CI

  return(c(c1, c2, c3, w1, w2, w3))
}

# set.seed(123)
# x = runif(n, -1, 1); y = sapply(x,  mu) + rnorm(n, sd = 0.5)
# par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))
# plot(x, y, pch = 16, cex = 0.7, xlim = c(-1,1), col = c(mred, mgreen)[1+as.numeric(x>0)])
# abline(v = 0, col = "black", lty = 1, lwd = 1.5)
# xx = seq(min(x), max(x), length = n)
# yy = sapply(xx, mu)
# points(xx, yy, type = "l", lwd = 2.5)
# plot(our <- plrd(y, x, cutoff = 0))

# results.list <- pbmclapply(1:num_replications, single_experiment, mc.cores = num_cores)
# results <- do.call(rbind, results.list)
#
# print(matrix(colMeans(results), ncol = 2,
#              dimnames = list(c("rdrobust conv", "rdrobust biascorr", "rdrobust robust",
#                                "rdhonest", "plrd"),
#                              c("coverage", "avg width")
#              )
# ), digits = 3
# )

#===========================================
#           Pure noise example
#===========================================

single_experiment_pure_noise <- function(i, n) {
  truth <- 0

  set.seed(i)
  x = runif(n, -1, 1); y = rnorm(n) # uniform(-1,1) running variable & pure noise response

  plrd.out <- plrd(y, x, cutoff = 0) # PLRD estimates and CI's
  rdrob = ests.from.ci(rdrobust(y, x, c = 0)$ci) # rdrobust estimates and CI's
  rdh = RDHonest(y ~ x, cutoff = 0)$coefficients # RDHonest estimates and CI's

  w1 = 2*rdrob[,2] # widths of rdrobust CI's (conventional, bias corrected, robust)
  w2 = rdh$conf.high-rdh$conf.low # width of RDHonest CI
  w3 <- 2*plrd.out$tau.plusminus # width of PLRD CI

  c1 = apply(rdrob, 1, function(x) as.numeric(abs(x[1]-truth) < x[2])) # coverages of rdrobust CI's (conventional, bias corrected, robust)
  c2 = as.numeric((rdh$conf.low<truth)*(rdh$conf.high>truth)) # coverage of RDHonest CI
  c3 <- as.numeric(abs(plrd.out$tau.hat - truth) <= w3/2) # coverage of PLRD CI

  return(c(c1, c2, c3, w1, w2, w3))
}

# results <- pbmclapply(1:num_replications, single_experiment_pure_noise, mc.cores = num_cores)
# results.pure.noise <- do.call(rbind, results)
#
# print(matrix(colMeans(results.pure.noise), ncol = 2,
#              dimnames = list(c("rdrobust conv", "rdrobust biascorr", "rdrobust robust",
#                                "rdhonest", "plrd"),
#                              c("coverage", "avg width")
#                              )
#              ), digits = 3
# )
#
# write.csv(results.pure.noise, "expt_pure_noise.csv", row.names = FALSE)


sigma.eps = 0.1295 # this noise level is from Calonico et al. (2014)

#===========================================
#           Setting 1 from CCT
#===========================================

mu_1 <- function(x){
  if(x < 0) return(0.48 + 1.27*x + 7.18*x^2 + 20.21*x^3 + 21.54*x^4 + 7.33*x^5)
  else return(0.52 + 0.84*x -3*x^2+ 7.99*x^3 - 9.01*x^4 + 3.56*x^5)
}

single_experiment_CCT1 <- function(i, n = 500) {
  truth <-  mu_1(1e-9) - mu_1(-1e-9)

  set.seed(i)
  x = 2*rbeta(n, 2, 4) - 1; y = sapply(x,  mu_1) + rnorm(n, sd = sigma.eps)

  plrd.out <- plrd(y, x, cutoff = 0) # PLRD estimates and CI's
  rdrob = ests.from.ci(rdrobust(y, x, c = 0, vce = 'hc1')$ci) # rdrobust estimates and CI's
  rdh = RDHonest(y ~ x, cutoff = 0)$coefficients # RDHonest estimates and CI's

  w1 = 2*rdrob[,2] # widths of rdrobust CI's (conventional, bias corrected, robust)
  w2 = rdh$conf.high-rdh$conf.low # width of RDHonest CI
  w3 <- 2*plrd.out$tau.plusminus # width of PLRD CI

  c1 = apply(rdrob, 1, function(x) as.numeric(abs(x[1]-truth) < x[2])) # coverages of rdrobust CI's (conventional, bias corrected, robust)
  c2 = as.numeric((rdh$conf.low<truth)*(rdh$conf.high>truth)) # coverage of RDHonest CI
  c3 <- as.numeric(abs(plrd.out$tau.hat - truth) <= w3/2) # coverage of PLRD CI

  return(c(c1, c2, c3, w1, w2, w3))
}

results <- pbmclapply(1:num_replications, function(i) single_experiment_CCT1(i,n), mc.cores = num_cores)
results.CCT1 <- do.call(rbind, results)

print(matrix(colMeans(results.CCT1), ncol = 2,
             dimnames = list(c("rdrobust conv", "rdrobust biascorr", "rdrobust robust",
                               "rdhonest", "plrd"),
                             c("coverage", "avg width")
             )
), digits = 3
)

# write.csv(results.CCT1, "expt_CCT1.csv", row.names = FALSE)

#===========================================
#           Setting 2 from CCT
#===========================================

mu_2 <- function(x){
  if(x < 0) return(3.71 + 2.30*x + 3.28*x^2 + 1.45*x^3 + 0.23*x^4 + 0.03*x^5)
  else return(0.26 + 18.49*x - 54.81*x^2 + 74.30*x^3 - 45.02*x^4 + 9.83*x^5)
}

single_experiment_CCT2 <- function(i, n = 500) {
  truth <-  mu_2(1e-9) - mu_2(-1e-9)

  set.seed(i)
  x = 2*rbeta(n, 2, 4) - 1; y = sapply(x,  mu_2) + rnorm(n, sd = sigma.eps)

  plrd.out <- plrd(y, x, cutoff = 0) # PLRD estimates and CI's
  rdrob = ests.from.ci(rdrobust(y, x, c = 0, vce = 'hc1')$ci) # rdrobust estimates and CI's
  rdh = RDHonest(y ~ x, cutoff = 0)$coefficients # RDHonest estimates and CI's
  plrd.out2 <- plrd(y, x, cutoff = 0, diff.curvatures = TRUE) # PLRD with curvature term

  w1 = 2*rdrob[,2] # widths of rdrobust CI's (conventional, bias corrected, robust)
  w2 = rdh$conf.high-rdh$conf.low # width of RDHonest CI
  w3 <- 2*plrd.out$tau.plusminus # width of PLRD CI
  w4 <- 2*plrd.out2$tau.plusminus # width of PLRD w/ diff curv CI

  c1 = apply(rdrob, 1, function(x) as.numeric(abs(x[1]-truth) < x[2])) # coverages of rdrobust CI's (conventional, bias corrected, robust)
  c2 = as.numeric((rdh$conf.low<truth)*(rdh$conf.high>truth)) # coverage of RDHonest CI
  c3 <- as.numeric(abs(plrd.out$tau.hat - truth) <= w3/2) # coverage of PLRD CI
  c4 <- as.numeric(abs(plrd.out2$tau.hat - truth) <= w4/2) # coverage of PLRD w/ diff curv CI

  return(c(c1, c2, c3, c4, w1, w2, w3, w4))
}

# results <- pbmclapply(1:num_replications, function(i) single_experiment_CCT2(i,n), mc.cores = num_cores)
# results.CCT2 <- do.call(rbind, results)
#
# print(matrix(colMeans(results.CCT2), ncol = 2,
#              dimnames = list(c("rdrobust conv", "rdrobust biascorr", "rdrobust robust",
#                                "rdhonest", "plrd", "plrd w/ diff curv"),
#                              c("coverage", "avg width")
#              )
# ), digits = 3
# )

# write.csv(results.CCT2, "expt_CCT2.csv", row.names = FALSE)


#===========================================
#           Setting 3 from CCT
#===========================================

mu_3 <- function(x){
  if(x < 0) return(0.48 + 1.27*x - 0.5*7.18*x^2 +
                     0.7*20.21*x^3 + 1.1*21.54*x^4 + 1.5*7.33*x^5)
  else return(0.52 + 0.84*x -0.1*3*x^2 - 0.3*7.99*x^3 - 0.1*9.01*x^4 + 3.56*x^5)
}

single_experiment_CCT3 <- function(i, n = 500) {
  truth <-  mu_3(1e-9) - mu_3(-1e-9)

  set.seed(i)
  x = 2*rbeta(n, 2, 4) - 1; y = sapply(x,  mu_3) + rnorm(n, sd = sigma.eps)

  plrd.out <- plrd(y, x, cutoff = 0) # PLRD estimates and CI's
  rdrob = ests.from.ci(rdrobust(y, x, c = 0, vce = 'hc1')$ci) # rdrobust estimates and CI's
  rdh = RDHonest(y ~ x, cutoff = 0)$coefficients # RDHonest estimates and CI's

  w1 = 2*rdrob[,2] # widths of rdrobust CI's (conventional, bias corrected, robust)
  w2 = rdh$conf.high-rdh$conf.low # width of RDHonest CI
  w3 <- 2*plrd.out$tau.plusminus # width of PLRD CI

  c1 = apply(rdrob, 1, function(x) as.numeric(abs(x[1]-truth) < x[2])) # coverages of rdrobust CI's (conventional, bias corrected, robust)
  c2 = as.numeric((rdh$conf.low<truth)*(rdh$conf.high>truth)) # coverage of RDHonest CI
  c3 <- as.numeric(abs(plrd.out$tau.hat - truth) <= w3/2) # coverage of PLRD CI

  return(c(c1, c2, c3, w1, w2, w3))
}

# results <- pbmclapply(1:num_replications, function(i) single_experiment_CCT3(i,n), mc.cores = num_cores)
# results.CCT3 <- do.call(rbind, results)
#
# print(matrix(colMeans(results.CCT3), ncol = 2,
#              dimnames = list(c("rdrobust conv", "rdrobust biascorr", "rdrobust robust",
#                                "rdhonest", "plrd"),
#                              c("coverage", "avg width")
#              )
# ), digits = 3
# )
#
# write.csv(results.CCT3, "expt_CCT3.csv", row.names = FALSE)

#=========================================
#         Setting 4 from IK 2012
#=========================================

mu_4 <- function(x) return((3+as.numeric(x>0))*x^2)

single_experiment_IK <- function(i, n = 500) {
  truth <-  mu_4(1e-9) - mu_4(-1e-9)

  set.seed(i)
  x = 2*rbeta(n, 2, 4) - 1; y = sapply(x,  mu_4) + rnorm(n, sd = sigma.eps)

  plrd.out <- plrd(y, x, cutoff = 0) # PLRD estimates and CI's
  rdrob = ests.from.ci(rdrobust(y, x, c = 0, vce = 'hc1')$ci) # rdrobust estimates and CI's
  rdh = RDHonest(y ~ x, cutoff = 0)$coefficients # RDHonest estimates and CI's

  w1 = 2*rdrob[,2] # widths of rdrobust CI's (conventional, bias corrected, robust)
  w2 = rdh$conf.high-rdh$conf.low # width of RDHonest CI
  w3 <- 2*plrd.out$tau.plusminus # width of PLRD CI

  c1 = apply(rdrob, 1, function(x) as.numeric(abs(x[1]-truth) < x[2])) # coverages of rdrobust CI's (conventional, bias corrected, robust)
  c2 = as.numeric((rdh$conf.low<truth)*(rdh$conf.high>truth)) # coverage of RDHonest CI
  c3 <- as.numeric(abs(plrd.out$tau.hat - truth) <= w3/2) # coverage of PLRD CI

  return(c(c1, c2, c3, w1, w2, w3))
}

# results <- pbmclapply(1:num_replications, function(i) single_experiment_IK(i,n), mc.cores = num_cores)
# results.IK <- do.call(rbind, results)
#
# print(matrix(colMeans(results.IK), ncol = 2,
#              dimnames = list(c("rdrobust conv", "rdrobust biascorr", "rdrobust robust",
#                                "rdhonest", "plrd"),
#                              c("coverage", "avg width")
#              )
# ), digits = 3
# )
#
# write.csv(results.IK, "expt_IK.csv", row.names = FALSE)

#------------------------------------------------------------------------------------------------
#
# print(matrix(colMeans(read.csv("expt_pure_noise.csv")), ncol = 2,
#              dimnames = list(c("rdrobust conv", "rdrobust biascorr", "rdrobust robust",
#                                "rdhonest", "plrd"),
#                              c("coverage", "avg width")
#              )
# ), digits = 3
# )
#
# print(matrix(colMeans(read.csv("expt_CCT1.csv")), ncol = 2,
#              dimnames = list(c("rdrobust conv", "rdrobust biascorr", "rdrobust robust",
#                                "rdhonest", "plrd"),
#                              c("coverage", "avg width")
#              )
# ), digits = 3
# )
#
# print(matrix(colMeans(read.csv("expt_CCT2.csv")), ncol = 2,
#              dimnames = list(c("rdrobust conv", "rdrobust biascorr", "rdrobust robust",
#                                "rdhonest", "plrd", "plrd w/ diff curv"),
#                              c("coverage", "avg width")
#              )
# ), digits = 3
# )
#
# print(matrix(colMeans(read.csv("expt_CCT3.csv")), ncol = 2,
#              dimnames = list(c("rdrobust conv", "rdrobust biascorr", "rdrobust robust",
#                                "rdhonest", "plrd"),
#                              c("coverage", "avg width")
#              )
# ), digits = 3
# )
#
# print(matrix(colMeans(read.csv("expt_IK.csv")), ncol = 2,
#              dimnames = list(c("rdrobust conv", "rdrobust biascorr", "rdrobust robust",
#                                "rdhonest", "plrd"),
#                              c("coverage", "avg width")
#              )
# ), digits = 3
# )


#------------------------------------------------------------------------------------------------


run_experiments <- function(num_replications = 100, n = 500,
                            num_cores = detectCores() - 1,
                            seed = 42){

  mred = "#CC3311"; mgreen = "#009E73"

  par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))
  withr::with_seed(seed, {x = runif(n, -1, 1); y = rnorm(n)})
  plot(x, y, pch = 16, cex = 0.7, xlim = c(-1,1),
       col = c(mred, mgreen)[1+as.numeric(x>0)], main = "Setting 0")
  abline(v = 0, col = "black", lty = 1, lwd = 1.5)
  abline(h = 0, col = "black", lty = 1, lwd = 2.5)

  cat("Running simulation setting 0, where X is uniform and Y is pure noise\n")
  results <- pbmclapply(1:num_replications,
                        function(i) single_experiment_pure_noise(i,n),
                        mc.cores = num_cores)
  results.pure.noise <- do.call(rbind, results)
  print(matrix(colMeans(results.pure.noise), ncol = 2,
               dimnames = list(c("rdrobust conv", "rdrobust biascorr", "rdrobust robust",
                                 "rdhonest", "plrd"),
                               c("coverage", "avg width")
               )
  ), digits = 3
  )
  cat("\n")

  x = 2*rbeta(n, 2, 4) - 1; y = sapply(x, mu_1) + rnorm(n, sd = sigma.eps)
  plot(x, y, pch = 16, cex = 0.7, xlim = c(-1,1),
       col = c(mred, mgreen)[1+as.numeric(x>0)], main = "Setting 1")
  abline(v = 0, col = "black", lty = 1, lwd = 1.5)
  xx = seq(min(x), max(x), length = 500)
  yy = sapply(xx, mu_1)
  points(xx, yy, type = "l", lwd = 2.5)

  cat("Running simulation setting 1 from Calonico et al. (2014)\n")
  results <- pbmclapply(1:num_replications,
                        function(i) single_experiment_CCT1(i, n),
                        mc.cores = num_cores)
  results.CCT1 <- do.call(rbind, results)
  print(matrix(colMeans(results.CCT1), ncol = 2,
               dimnames = list(c("rdrobust conv", "rdrobust biascorr", "rdrobust robust",
                                 "rdhonest", "plrd"),
                               c("coverage", "avg width")
               )
  ), digits = 3
  )
  cat("\n")

  x = 2*rbeta(n, 2, 4) - 1; y = sapply(x, mu_2) + rnorm(n, sd = sigma.eps)
  plot(x, y, pch = 16, cex = 0.7, xlim = c(-1,1),
       col = c(mred, mgreen)[1+as.numeric(x>0)], main = "Setting 2")
  abline(v = 0, col = "black", lty = 1, lwd = 1.5)
  xx = seq(min(x), max(x), length = 500)
  yy = sapply(xx, mu_2)
  points(xx, yy, type = "l", lwd = 2.5)

  cat("Running simulation setting 2 from Calonico et al. (2014)\n")
  results <- pbmclapply(1:num_replications,
                        function(i) single_experiment_CCT2(i,n),
                        mc.cores = num_cores)
  results.CCT2 <- do.call(rbind, results)
  print(matrix(colMeans(results.CCT2), ncol = 2,
               dimnames = list(c("rdrobust conv", "rdrobust biascorr", "rdrobust robust",
                                 "rdhonest", "plrd", "plrd w/ diff curvature"),
                               c("coverage", "avg width")
               )
  ), digits = 3
  )
  cat("\n")

  x = 2*rbeta(n, 2, 4) - 1; y = sapply(x, mu_3) + rnorm(n, sd = sigma.eps)
  plot(x, y, pch = 16, cex = 0.7, xlim = c(-1,1),
       col = c(mred, mgreen)[1+as.numeric(x>0)], main = "Setting 3")
  abline(v = 0, col = "black", lty = 1, lwd = 1.5)
  xx = seq(min(x), max(x), length = 500)
  yy = sapply(xx, mu_3)
  points(xx, yy, type = "l", lwd = 2.5)

  cat("Running simulation setting 3 from Calonico et al. (2014)\n")
  results <- pbmclapply(1:num_replications,
                        function(i) single_experiment_CCT3(i, n),
                        mc.cores = num_cores)
  results.CCT3 <- do.call(rbind, results)
  print(matrix(colMeans(results.CCT3), ncol = 2,
               dimnames = list(c("rdrobust conv", "rdrobust biascorr", "rdrobust robust",
                                 "rdhonest", "plrd"),
                               c("coverage", "avg width")
               )
  ), digits = 3
  )
  cat("\n")

  x = 2*rbeta(n, 2, 4) - 1; y = sapply(x, mu_4) + rnorm(n, sd = sigma.eps)
  plot(x, y, pch = 16, cex = 0.7, xlim = c(-1,1),
       col = c(mred, mgreen)[1+as.numeric(x>0)], main = "Setting 4")
  abline(v = 0, col = "black", lty = 1, lwd = 1.5)
  xx = seq(min(x), max(x), length = 500)
  yy = sapply(xx, mu_4)
  points(xx, yy, type = "l", lwd = 2.5)

  cat("Running simulation setting 4: an experiment from Imbens and Kalyanamaran (2012)\n")
  results <- pbmclapply(1:num_replications,
                        function(i) single_experiment_IK(i,n),
                        mc.cores = num_cores)
  results.IK <- do.call(rbind, results)
  print(matrix(colMeans(results.IK), ncol = 2,
               dimnames = list(c("rdrobust conv", "rdrobust biascorr", "rdrobust robust",
                                 "rdhonest", "plrd"),
                               c("coverage", "avg width")
               )
  ), digits = 3
  )
}


run_experiments(100)
