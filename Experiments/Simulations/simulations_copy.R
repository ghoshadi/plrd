if (!requireNamespace("RDHonest", quietly = TRUE)) install.packages("RDHonest")
if (!requireNamespace("rdrobust", quietly = TRUE)) install.packages("rdrobust")
if (!requireNamespace("glmnet", quietly = TRUE)) install.packages("glmnet")
if (!requireNamespace("splines", quietly = TRUE)) install.packages("splines")
if (!requireNamespace("parallel", quietly = TRUE)) install.packages("parallel")
if (!requireNamespace("R.matlab", quietly = TRUE)) install.packages("R.matlab")
library(RDHonest) # version 1.0.0
library(rdrobust) # version 2.2
library(glmnet) # version 4.1.8
library(splines) # version 4.4.2
library(parallel) # version 4.4.2
library(R.matlab) # version 3.7.0
library(plrd) # our package
mred = "#CC3311"; mgreen = "#009E73"
splrd.outce('ploptrdd.R')

#===========================================
#            Create your own
#===========================================

num_cores <- detectCores() - 1
n = 500
num_iterations <- 1000
mu <- function(x){
  return(1 + 2*as.numeric(x>0) + 3*x^2)
}
truth <-  mu(1e-9) - mu(-1e-9)

process_file <- function(i) {
  x = runif(n, -1, 1); y = sapply(x, mu) + rnorm(n, sd = 0.5)
  plrd.out <- plrd(y, x, 0, verbose = FALSE) # PLRD estimates and CI's
  rdrob = ests.from.ci(rdrobust(y, x, c = 0, vce = 'hc1')$ci) # rdrobust estimates and CI's
  rdh = RDHonest(y ~ x)$coefficients # RDHonest estimates and CI's

  w1 = 2*rdrob[,2] # widths of rdrobust CI's (conventional, bias corrected, robust)
  w2 = rdh$conf.high-rdh$conf.low # width of RDHonest CI
  w3 <- 2*plrd.out$tau.plusminus # width of PLRD CI

  c1 = apply(rdrob, 1, function(x) as.numeric(abs(x[1]-truth) < x[2])) # coverages of rdrobust CI's (conventional, bias corrected, robust)
  c2 = as.numeric((rdh$conf.low<truth)*(rdh$conf.high>truth)) # coverage of RDHonest CI
  c3 <- as.numeric(abs(plrd.out$tau.hat - truth) <= w3/2) # coverage of PLRD CI
  return(c(c1, c2, c3, w1, w2, w3))
}

set.seed(123)
x = runif(n, -1, 1); y = sapply(x, mu) + rnorm(n, sd = 0.5)
plot(x, y, pch = 16, cex = 0.7, xlim = c(-1,1), col = c(mred, mgreen)[1+as.numeric(x>0)])
abline(v = 0, col = "black", lty = 1, lwd = 1.5)
xx = seq(min(x), max(x), length = n)
yy = sapply(xx, mu)
points(xx, yy, type = "l", lwd = 2.5)
plot(plrd.out <- plrd(y, x, 0))

system.time({
  results <- do.call(rbind, mclapply(1:num_iterations, process_file, mc.cores = num_cores))
})

print(matrix(colMeans(results), ncol = 2,
             dimnames = list(c("rdrobust conv", "rdrobust biascorr", "rdrobust robust",
                               "rdhonest", "plrd"),
                             c("coverage", "avg width")
             )
), digits = 3
)

#-----------------------------------------
#          Pure noise example
#-----------------------------------------

n = 1000
truth <- 0
num_cores <- detectCores() - 1

process_file <- function(i) {
  x = runif(n, -1, 1); y = rnorm(n) # uniform(-1,1) running variable & pure noise response
  plrd.out <- plrd(y, x, 0, diff.curvatures = T, verbose = FALSE) # PLRD estimates and CI's
  ploptrdd.out <- ploptrdd(y, x, 0, verbose = FALSE)
  rdrob = ests.from.ci(rdrobust(y, x, c = 0, vce = 'hc1')$ci) # rdrobust estimates and CI's
  rdrob2 = ests.from.ci(rdrobust(y, x, c = 0, p = 2, vce = 'hc1')$ci)
  w1 = 2*rdrob[1,2]*2.181/1.96 # Column (a)
  w2 = 2*rdrob2[1,2]*2.113/1.96 # Column (b)
  w3 <- 2*plrd.out$tau.plusminus # Column (c)
  w4 <- 2*ploptrdd.out$tau.plusminus # Column (d)

  c1 <- as.numeric(abs(rdrob[1,1] - truth) <= w1/2)
  c2 <- as.numeric(abs(rdrob[1,1] - truth) <= w2/2)
  c3 <- as.numeric(abs(plrd.out$tau.hat - truth) <= w3/2) # coverage of PLRD CI w/ diff curv
  c4 <- as.numeric(abs(ploptrdd.out$tau.hat - truth) <= w4/2) # coverage of PLOPTRDD CI
  return(c(c1, c2, c3, c4, w1, w2, w3, w4))
}

set.seed(123)
system.time({
  results <- mclapply(1:num_iterations, process_file, mc.cores = num_cores)
  results.pure.noise <- do.call(rbind, results)
})

print(matrix(colMeans(results.pure.noise), ncol = 2,
             dimnames = list(paste0("Alt.(",letters[1:4],")"),
                             c("coverage", "avg width")
             )
), digits = 3
)
write.csv(results.pure.noise, "expt_pure_noise_appendix.csv", row.names = FALSE)



#-----------------------------------------
#           Setting 1 from CCT
#-----------------------------------------

n <- 500
mu_1 <- function(x){
  if(x < 0) return(0.48 + 1.27*x + 7.18*x^2 + 20.21*x^3 + 21.54*x^4 + 7.33*x^5)
  else return(0.52 + 0.84*x -3*x^2+ 7.99*x^3 - 9.01*x^4 + 3.56*x^5)
}
truth <-  mu_1(1e-9) - mu_1(-1e-9)
#
process_file <- function(i) {
  x = 2*rbeta(n, 2, 4) - 1; y = sapply(x, mu_1) + rnorm(n, sd = 0.1295)
  plrd.out <- plrd(y, x, 0, diff.curvatures = T, verbose = FALSE) # PLRD estimates and CI's
  ploptrdd.out <- ploptrdd(y, x, 0, verbose = FALSE)
  rdrob = ests.from.ci(rdrobust(y, x, c = 0, vce = 'hc1')$ci) # rdrobust estimates and CI's
  rdrob2 = ests.from.ci(rdrobust(y, x, c = 0, p = 2, vce = 'hc1')$ci)
  w1 = 2*rdrob[1,2]*2.181/1.96 # Column (a)
  w2 = 2*rdrob2[1,2]*2.113/1.96 # Column (b)
  w3 <- 2*plrd.out$tau.plusminus # Column (c)
  w4 <- 2*ploptrdd.out$tau.plusminus # Column (d)

  c1 <- as.numeric(abs(rdrob[1,1] - truth) <= w1/2)
  c2 <- as.numeric(abs(rdrob[1,1] - truth) <= w2/2)
  c3 <- as.numeric(abs(plrd.out$tau.hat - truth) <= w3/2) # coverage of PLRD CI w/ diff curv
  c4 <- as.numeric(abs(ploptrdd.out$tau.hat - truth) <= w4/2) # coverage of PLOPTRDD CI
  return(c(c1, c2, c3, c4, w1, w2, w3, w4))
}

set.seed(123)
system.time({
  results <- mclapply(1:num_iterations, process_file, mc.cores = num_cores)
  results.CCT1 <- do.call(rbind, results)
})

print(matrix(colMeans(results.CCT1), ncol = 2,
             dimnames = list(paste0("Alt.(",letters[1:4],")"),
                             c("coverage", "avg width")
             )
), digits = 3
)
write.csv(results.CCT1, "expt_CCT1_appendix.csv", row.names = FALSE)

#-----------------------------------------
#           Setting 2 from CCT
#-----------------------------------------

n <- 500
mu_2 <- function(x){
  if(x < 0) return(3.71 + 2.30*x + 3.28*x^2 + 1.45*x^3 + 0.23*x^4 + 0.03*x^5)
  else return(0.26 + 18.49*x - 54.81*x^2 + 74.30*x^3 - 45.02*x^4 + 9.83*x^5)
}
truth <-  mu_2(1e-9) - mu_2(-1e-9)

process_file <- function(i) {
  x = 2*rbeta(n, 2, 4) - 1; y = sapply(x, mu_2) + rnorm(n, sd = 0.1295)
  plrd.out <- plrd(y, x, 0, diff.curvatures = T, verbose = FALSE) # PLRD estimates and CI's
  ploptrdd.out <- ploptrdd(y, x, 0, verbose = FALSE)
  rdrob = ests.from.ci(rdrobust(y, x, c = 0, vce = 'hc1')$ci) # rdrobust estimates and CI's
  rdrob2 = ests.from.ci(rdrobust(y, x, c = 0, p = 2, vce = 'hc1')$ci)
  w1 = 2*rdrob[1,2]*2.181/1.96 # Column (a)
  w2 = 2*rdrob2[1,2]*2.113/1.96 # Column (b)
  w3 <- 2*plrd.out$tau.plusminus # Column (c)
  w4 <- 2*ploptrdd.out$tau.plusminus # Column (d)

  c1 <- as.numeric(abs(rdrob[1,1] - truth) <= w1/2)
  c2 <- as.numeric(abs(rdrob[1,1] - truth) <= w2/2)
  c3 <- as.numeric(abs(plrd.out$tau.hat - truth) <= w3/2) # coverage of PLRD CI w/ diff curv
  c4 <- as.numeric(abs(ploptrdd.out$tau.hat - truth) <= w4/2) # coverage of PLOPTRDD CI
  return(c(c1, c2, c3, c4, w1, w2, w3, w4))
}

set.seed(123)
system.time({
  results <- mclapply(1:num_iterations, process_file, mc.cores = num_cores)
  results.CCT2 <- do.call(rbind, results)
})

print(matrix(colMeans(results.CCT2), ncol = 2,
             dimnames = list(paste0("Alt.(",letters[1:4],")"),
                             c("coverage", "avg width")
             )
), digits = 3
)
write.csv(results.CCT2, "expt_CCT2_appendix.csv", row.names = FALSE)


#-----------------------------------------
#           Setting 3 from CCT
#-----------------------------------------

n <- 500
mu_3 <- function(x){
  if(x < 0) return(0.48 + 1.27*x - 0.5*7.18*x^2 +
                     0.7*20.21*x^3 + 1.1*21.54*x^4 + 1.5*7.33*x^5)
  else return(0.52 + 0.84*x -0.1*3*x^2 - 0.3*7.99*x^3 - 0.1*9.01*x^4 + 3.56*x^5)
}
truth <-  mu_3(1e-9) - mu_3(-1e-9)

process_file <- function(i) {
  x = 2*rbeta(n, 2, 4) - 1; y = sapply(x, mu_3) + rnorm(n, sd = 0.1295)
  plrd.out <- plrd(y, x, 0, diff.curvatures = T, verbose = FALSE) # PLRD estimates and CI's
  ploptrdd.out <- ploptrdd(y, x, 0, verbose = FALSE)
  rdrob = ests.from.ci(rdrobust(y, x, c = 0, vce = 'hc1')$ci) # rdrobust estimates and CI's
  rdrob2 = ests.from.ci(rdrobust(y, x, c = 0, p = 2, vce = 'hc1')$ci)
  w1 = 2*rdrob[1,2]*2.181/1.96 # Column (a)
  w2 = 2*rdrob2[1,2]*2.113/1.96 # Column (b)
  w3 <- 2*plrd.out$tau.plusminus # Column (c)
  w4 <- 2*ploptrdd.out$tau.plusminus # Column (d)

  c1 <- as.numeric(abs(rdrob[1,1] - truth) <= w1/2)
  c2 <- as.numeric(abs(rdrob[1,1] - truth) <= w2/2)
  c3 <- as.numeric(abs(plrd.out$tau.hat - truth) <= w3/2) # coverage of PLRD CI w/ diff curv
  c4 <- as.numeric(abs(ploptrdd.out$tau.hat - truth) <= w4/2) # coverage of PLOPTRDD CI
  return(c(c1, c2, c3, c4, w1, w2, w3, w4))
}

set.seed(123)
system.time({
  results <- mclapply(1:num_iterations, process_file, mc.cores = num_cores)
  results.CCT3 <- do.call(rbind, results)
})

print(matrix(colMeans(results.CCT3), ncol = 2,
             dimnames = list(paste0("Alt.(",letters[1:4],")"),
                             c("coverage", "avg width")
             )
), digits = 3
)
write.csv(results.CCT3, "expt_CCT3_appendix.csv", row.names = FALSE)


#-----------------------------------------
#           Setting 4 from IK
#-----------------------------------------

n <- 500
mu_4 <- function(x) return((3+as.numeric(x>0))*x^2)
truth <-  mu_4(1e-9) - mu_4(-1e-9)

process_file <- function(i) {
  x = 2*rbeta(n, 2, 4) - 1; y = sapply(x, mu_4) + rnorm(n, sd = 0.1295)
  plrd.out <- plrd(y, x, 0, diff.curvatures = T, verbose = FALSE) # PLRD estimates and CI's
  ploptrdd.out <- ploptrdd(y, x, 0, verbose = FALSE)
  rdrob = ests.from.ci(rdrobust(y, x, c = 0, vce = 'hc1')$ci) # rdrobust estimates and CI's
  rdrob2 = ests.from.ci(rdrobust(y, x, c = 0, p = 2, vce = 'hc1')$ci)
  w1 = 2*rdrob[1,2]*2.181/1.96 # Column (a)
  w2 = 2*rdrob2[1,2]*2.113/1.96 # Column (b)
  w3 <- 2*plrd.out$tau.plusminus # Column (c)
  w4 <- 2*ploptrdd.out$tau.plusminus # Column (d)

  c1 <- as.numeric(abs(rdrob[1,1] - truth) <= w1/2)
  c2 <- as.numeric(abs(rdrob[1,1] - truth) <= w2/2)
  c3 <- as.numeric(abs(plrd.out$tau.hat - truth) <= w3/2) # coverage of PLRD CI w/ diff curv
  c4 <- as.numeric(abs(ploptrdd.out$tau.hat - truth) <= w4/2) # coverage of PLOPTRDD CI
  return(c(c1, c2, c3, c4, w1, w2, w3, w4))
}

set.seed(123)
system.time({
  results <- mclapply(1:num_iterations, process_file, mc.cores = num_cores)
  results.IK <- do.call(rbind, results)
})

print(matrix(colMeans(results.IK), ncol = 2,
             dimnames = list(paste0("Alt.(",letters[1:4],")"),
                             c("coverage", "avg width")
             )
), digits = 3
)
write.csv(results.IK, "expt_IK_appendix.csv", row.names = FALSE)

print(t(matrix(colMeans(read.csv("expt_CCT2_appendix.csv")), ncol = 2,
             dimnames = list(paste0("Alt.(",letters[1:4],")"),
                             c("coverage", "avg width")
             )
)), digits = 3
)

#------------------------------------------------------------------------------------------------

print(matrix(colMeans(read.csv("expt_pure_noise.csv")), ncol = 2,
             dimnames = list(c("rdrobust conv", "rdrobust biascorr", "rdrobust robust",
                               "rdhonest", "plrd"),
                             c("coverage", "avg width")
             )
), digits = 3
)

print(matrix(colMeans(read.csv("expt_CCT1.csv")), ncol = 2,
             dimnames = list(c("rdrobust conv", "rdrobust biascorr", "rdrobust robust",
                               "rdhonest", "plrd"),
                             c("coverage", "avg width")
             )
), digits = 3
)

print(matrix(colMeans(read.csv("expt_CCT2.csv")), ncol = 2,
             dimnames = list(c("rdrobust conv", "rdrobust biascorr", "rdrobust robust",
                               "rdhonest", "plrd", "plrd w/ diff curv"),
                             c("coverage", "avg width")
             )
), digits = 3
)

print(matrix(colMeans(read.csv("expt_CCT3.csv")), ncol = 2,
             dimnames = list(c("rdrobust conv", "rdrobust biascorr", "rdrobust robust",
                               "rdhonest", "plrd"),
                             c("coverage", "avg width")
             )
), digits = 3
)

print(matrix(colMeans(read.csv("expt_IK.csv")), ncol = 2,
             dimnames = list(c("rdrobust conv", "rdrobust biascorr", "rdrobust robust",
                               "rdhonest", "plrd"),
                             c("coverage", "avg width")
             )
), digits = 3
)

#------------------------------------------------------------------------------------------------

## Plotting one realizations for each of the above simulation settings

mred = "#CC3311"; mgreen = "#009E73"

par(mar=c(4.1,4.1,1,1))
set.seed(123)
n = 500
x = 2*rbeta(n, 2, 4) - 1; y = sapply(x, mu_1) + rnorm(n, sd = 0.1295)
plot(x, y, pch = 16, cex = 0.7, xlim = c(-1,1),
     col = c(mred, mgreen)[1+as.numeric(x>0)])
abline(v = 0, col = "black", lty = 1, lwd = 1.5)
xx = seq(min(x), max(x), length = 500)
yy = sapply(xx, mu_1)
points(xx, yy, type = "l", lwd = 2.5)


set.seed(123)
n = 500
x = 2*rbeta(n, 2, 4) - 1; y = sapply(x, mu_2) + rnorm(n, sd = 0.1295)
plot(x, y, pch = 16, cex = 0.7, xlim = c(-1,1),
     col = c(mred, mgreen)[1+as.numeric(x>0)])
abline(v = 0, col = "black", lty = 1, lwd = 1.5)
xx = seq(min(x), max(x), length = 500)
yy = sapply(xx, mu_2)
points(xx, yy, type = "l", lwd = 2.5)


set.seed(123)
n = 500
x = 2*rbeta(n, 2, 4) - 1; y = sapply(x, mu_3) + rnorm(n, sd = 0.1295)
plot(x, y, pch = 16, cex = 0.7, xlim = c(-1,1),
     col = c(mred, mgreen)[1+as.numeric(x>0)])
abline(v = 0, col = "black", lty = 1, lwd = 1.5)
xx = seq(min(x), max(x), length = 500)
yy = sapply(xx, mu_3)
points(xx, yy, type = "l", lwd = 2.5)


set.seed(123)
n = 500
x = 2*rbeta(n, 2, 4) - 1; y = sapply(x, mu_4) + rnorm(n, sd = 0.1295)
plot(x, y, pch = 16, cex = 0.7, xlim = c(-1,1),
     col = c(mred, mgreen)[1+as.numeric(x>0)])
abline(v = 0, col = "black", lty = 1, lwd = 1.5)
xx = seq(min(x), max(x), length = 500)
yy = sapply(xx, mu_4)
points(xx, yy, type = "l", lwd = 2.5)
