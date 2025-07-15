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
source('ploptrdd.R') # temporary function implementing OPTRDD under partially linear model



# Meyersson data
#
# dfa <- read.csv("./fake_Meyersson/dfa.csv")
# dfb <- read.csv("./fake_Meyersson/dfb.csv")
# print(truth <- mean(dfa$y) - mean(dfb$y))
#
# num_iterations = 1000
# num_cores <- detectCores() - 1
#
# process_file <- function(i) {
#   df <- read.csv(paste0("./fake_Meyersson/df_fake_", i, ".csv"))
#   plrd.out <- plrd(df$y, df$x, verbose = FALSE) # PLRD estimates and CI's.
#   plrd.out2<- plrd(df$y, df$x, diff.curvatures = TRUE, verbose = FALSE)
#   ploptrdd.out <- ploptrdd(df$y, df$x, verbose = FALSE)
#   rdrob = ests.from.ci(rdrobust(df$y, df$x, c = 0)$ci) # rdrobust estimates and CI's
#   rdh = RDHonest(df$y ~ df$x)$coefficients # RDHonest estimates and CI's
#
#   w1 = 2*rdrob[,2] # widths of rdrobust CI's (conventional, bias corrected, robust)
#   w2 = rdh$conf.high-rdh$conf.low # width of RDHonest CI
#   w3 <- 2*plrd.out$tau.plusminus # width of PLRD CI
#   w5 = 2*rdrob[1,2]*2.181/1.96 # width of AK20 CI with 2.181 instead of 1.96
#
#   rdrob2 = ests.from.ci(rdrobust(df$y, df$x, c = 0, p = 2)$ci)
#   w4 = 2*rdrob2[1,2]*2.113/1.96 # width of AK20 CI with 2.113 instead of 1.96
#
#   w6 = 2*plrd.out2$tau.plusminus # width of OPTRDD with our smoothness
#   w7 = 2*ploptrdd.out$tau.plusminus # width of OPTRDD under partially linear model
#
#   c1 = apply(rdrob, 1, function(x) as.numeric(abs(x[1]-truth) <= x[2])) # coverages of rdrobust CI's (conventional, bias corrected, robust)
#   c2 = as.numeric((rdh$conf.low<=truth)*(rdh$conf.high>=truth)) # coverage of RDHonest CI
#   c3 <- as.numeric(abs(plrd.out$tau.hat - truth) <= w3/2) # coverage of PLRD CI
#   c4 <- as.numeric(abs(rdrob2[1,1] - truth) <= w4/2) # coverage of AK20 CI with 2.113 instead of 1.96
#   c5 <- as.numeric(abs(rdrob[1,1] - truth) <= w5/2) # coverage of AK20 CI with 2.181 instead of 1.96
#   c6 <- as.numeric(abs(plrd.out2$tau.hat - truth) <= w6/2) # coverage of OPTRDD with our smoothness
#   c7 <- as.numeric(abs(ploptrdd.out$tau.hat - truth) <= w7/2) # coverage of OPTRDD under partially linear model
#   return(c(c1, c2, c3, c4, c5, c6, c7, w1, w2, w3, w4, w5, w6, w7))
# }
#
# set.seed(123)
#
# system.time({
#   results <- mclapply(1:num_iterations, process_file, mc.cores = num_cores)
#   results <- do.call(rbind, results)
# })
# write.csv(results, "expt_wgan_Meyersson.csv", row.names = FALSE)
print(matrix(colMeans(read.csv("expt_wgan_Meyersson.csv")), ncol = 2,
             dimnames = list(c("rdrobust conv", "rdrobust biascorr", "rdrobust robust",
                               "rdhonest", "plrd", "AK20 2.113", "AK20 2.181",
                               "smoother optrdd", "pl+optrdd"),
                             c("coverage", "avg width")
             )
), digits = 4
)

# Lee (2008) data
#
# num_cores <- detectCores() - 1
# num_iterations <- 1000
#
# dfa <- read.csv("./fake_Lee08/dfa.csv")
# dfb <- read.csv("./fake_Lee08/dfb.csv")
# print(truth <- mean(dfa$y) - mean(dfb$y))
#
# process_file <- function(i) {
#   df <- read.csv(paste0("./fake_Lee08/df_fake_", i, ".csv"))
#   plrd.out <- plrd(df$y, df$x, verbose = FALSE) # PLRD estimates and CI's.
#   plrd.out2<- plrd(df$y, df$x, diff.curvatures = TRUE, verbose = FALSE) # OPTRDD with our smoothness
#   ploptrdd.out <- ploptrdd(df$y, df$x, verbose = FALSE) # OPTRDD under partially linear model
#   rdrob = ests.from.ci(rdrobust(df$y, df$x, c = 0)$ci) # rdrobust estimates and CI's
#   rdh = RDHonest(df$y ~ df$x)$coefficients # RDHonest estimates and CI's
#   rdrob2 = ests.from.ci(rdrobust(df$y, df$x, c = 0, p = 2)$ci)
#
#   w1 = 2*rdrob[,2] # widths of rdrobust CI's (conventional, bias corrected, robust)
#   w2 = rdh$conf.high-rdh$conf.low # width of RDHonest CI
#   w3 <- 2*plrd.out$tau.plusminus # width of PLRD CI
#   w4 = 2*rdrob2[1,2]*2.113/1.96 # width of AK20 CI with 2.113 instead of 1.96
#   w5 = 2*rdrob[1,2]*2.181/1.96 # width of AK20 CI with 2.181 instead of 1.96
#   w6 = 2*plrd.out2$tau.plusminus # width of OPTRDD with our smoothness
#   w7 = 2*ploptrdd.out$tau.plusminus # width of OPTRDD under partially linear model
#
#   c1 = apply(rdrob, 1, function(x) as.numeric(abs(x[1]-truth) <= x[2])) # coverages of rdrobust CI's (conventional, bias corrected, robust)
#   c2 = as.numeric((rdh$conf.low<=truth)*(rdh$conf.high>=truth)) # coverage of RDHonest CI
#   c3 <- as.numeric(abs(plrd.out$tau.hat - truth) <= w3/2) # coverage of PLRD CI
#   c4 <- as.numeric(abs(rdrob2[1,1] - truth) <= w4/2) # coverage of AK20 CI with 2.113 instead of 1.96
#   c5 <- as.numeric(abs(rdrob[1,1] - truth) <= w5/2) # coverage of AK20 CI with 2.181 instead of 1.96
#   c6 <- as.numeric(abs(plrd.out2$tau.hat - truth) <= w6/2) # coverage of OPTRDD with our smoothness
#   c7 <- as.numeric(abs(ploptrdd.out$tau.hat - truth) <= w7/2) # coverage of OPTRDD under partially linear model
#   return(c(c1, c2, c3, c4, c5, c6, c7, w1, w2, w3, w4, w5, w6, w7))
# }
#
# set.seed(123)
#
# system.time({
#   results <- mclapply(1:num_iterations, process_file, mc.cores = num_cores)
#   results <- do.call(rbind, results)
# })
# write.csv(results, "expt_wgan_lee08.csv", row.names = FALSE)

print(matrix(colMeans(read.csv("expt_wgan_lee08.csv")), ncol = 2,
             dimnames = list(c("rdrobust conv", "rdrobust biascorr", "rdrobust robust",
                               "rdhonest", "plrd", "AK20 2.113", "AK20 2.181",
                               "smoother optrdd", "pl+optrdd"),
                             c("coverage", "avg width")
             )
), digits = 3
)

# Matsudaira Reading
#
# dfa <- read.csv("./fake_mats_read/dfa.csv")
# dfb <- read.csv("./fake_mats_read/dfb.csv")
# print(truth <- mean(dfa$y) - mean(dfb$y))
#
# process_file <- function(i) {
#   df <- read.csv(paste0("./fake_mats_read/df_fake_", i, ".csv"))
#   plrd.out <- plrd(df$y, df$x, max.window = 50, verbose = FALSE) # PLRD estimates and CI's.
#   plrd.out2<- plrd(df$y, df$x, max.window = 50, diff.curvatures = TRUE, verbose = FALSE)
#   ploptrdd.out <- ploptrdd(df$y, df$x, max.window = 50, verbose = FALSE)
#   rdrob = ests.from.ci(rdrobust(df$y, df$x, c = 0)$ci) # rdrobust estimates and CI's
#   rdh = RDHonest(df$y ~ df$x)$coefficients # RDHonest estimates and CI's
#   rdrob2 = ests.from.ci(rdrobust(df$y, df$x, c = 0, p = 2)$ci)
#
#   w1 = 2*rdrob[,2] # widths of rdrobust CI's (conventional, bias corrected, robust)
#   w2 = rdh$conf.high-rdh$conf.low # width of RDHonest CI
#   w3 <- 2*plrd.out$tau.plusminus # width of PLRD CI
#   w4 = 2*rdrob2[1,2]*2.113/1.96 # width of AK20 CI with 2.113 instead of 1.96
#   w5 = 2*rdrob[1,2]*2.181/1.96 # width of AK20 CI with 2.181 instead of 1.96
#   w6 = 2*plrd.out2$tau.plusminus # width of OPTRDD with our smoothness
#   w7 = 2*ploptrdd.out$tau.plusminus # width of OPTRDD under partially linear model
#
#   c1 = apply(rdrob, 1, function(x) as.numeric(abs(x[1]-truth) <= x[2])) # coverages of rdrobust CI's (conventional, bias corrected, robust)
#   c2 = as.numeric((rdh$conf.low<=truth)*(rdh$conf.high>=truth)) # coverage of RDHonest CI
#   c3 <- as.numeric(abs(plrd.out$tau.hat - truth) <= w3/2) # coverage of PLRD CI
#   c4 <- as.numeric(abs(rdrob2[1,1] - truth) <= w4/2) # coverage of AK20 CI with 2.113 instead of 1.96
#   c5 <- as.numeric(abs(rdrob[1,1] - truth) <= w5/2) # coverage of AK20 CI with 2.181 instead of 1.96
#   c6 <- as.numeric(abs(plrd.out2$tau.hat - truth) <= w6/2) # coverage of OPTRDD with our smoothness
#   c7 <- as.numeric(abs(ploptrdd.out$tau.hat - truth) <= w7/2) # coverage of OPTRDD under partially linear model
#   return(c(c1, c2, c3, c4, c5, c6, c7, w1, w2, w3, w4, w5, w6, w7))
# }
#
# set.seed(123)
#
# system.time({
#   results <- mclapply(1:num_iterations, process_file, mc.cores = num_cores)
#   results <- do.call(rbind, results)
# })
# write.csv(results, "expt_wgan_mats_read.csv", row.names = FALSE)
print(matrix(colMeans(read.csv("expt_wgan_mats_read.csv")), ncol = 2,
             dimnames = list(c("rdrobust conv", "rdrobust biascorr", "rdrobust robust",
                               "rdhonest", "plrd", "AK20 2.113", "AK20 2.181",
                               "smoother optrdd", "pl+optrdd"),
                             c("coverage", "avg width")
             )
), digits = 3
)

# Matsudaira Math
#
# dfa <- read.csv("./fake_mats_math/dfa.csv")
# dfb <- read.csv("./fake_mats_math/dfb.csv")
# print(truth <- mean(dfa$y) - mean(dfb$y))
#
# process_file <- function(i) {
#   df <- read.csv(paste0("./fake_mats_math/df_fake_", i, ".csv"))
#   plrd.out <- plrd(df$y, df$x, max.window = 50, verbose = FALSE) # PLRD estimates and CI's.
#   plrd.out2<- plrd(df$y, df$x, max.window = 50, diff.curvatures = TRUE, verbose = FALSE)
#   ploptrdd.out <- ploptrdd(df$y, df$x, max.window = 50, verbose = FALSE)
#   rdrob = ests.from.ci(rdrobust(df$y, df$x, c = 0)$ci) # rdrobust estimates and CI's
#   rdh = RDHonest(df$y ~ df$x)$coefficients # RDHonest estimates and CI's
#   rdrob2 = ests.from.ci(rdrobust(df$y, df$x, c = 0, p = 2)$ci)
#
#   w1 = 2*rdrob[,2] # widths of rdrobust CI's (conventional, bias corrected, robust)
#   w2 = rdh$conf.high-rdh$conf.low # width of RDHonest CI
#   w3 <- 2*plrd.out$tau.plusminus # width of PLRD CI
#   w4 = 2*rdrob2[1,2]*2.113/1.96 # width of AK20 CI with 2.113 instead of 1.96
#   w5 = 2*rdrob[1,2]*2.181/1.96 # width of AK20 CI with 2.181 instead of 1.96
#   w6 = 2*plrd.out2$tau.plusminus # width of OPTRDD with our smoothness
#   w7 = 2*ploptrdd.out$tau.plusminus # width of OPTRDD under partially linear model
#
#   c1 = apply(rdrob, 1, function(x) as.numeric(abs(x[1]-truth) <= x[2])) # coverages of rdrobust CI's (conventional, bias corrected, robust)
#   c2 = as.numeric((rdh$conf.low<=truth)*(rdh$conf.high>=truth)) # coverage of RDHonest CI
#   c3 <- as.numeric(abs(plrd.out$tau.hat - truth) <= w3/2) # coverage of PLRD CI
#   c4 <- as.numeric(abs(rdrob2[1,1] - truth) <= w4/2) # coverage of AK20 CI with 2.113 instead of 1.96
#   c5 <- as.numeric(abs(rdrob[1,1] - truth) <= w5/2) # coverage of AK20 CI with 2.181 instead of 1.96
#   c6 <- as.numeric(abs(plrd.out2$tau.hat - truth) <= w6/2) # coverage of OPTRDD with our smoothness
#   c7 <- as.numeric(abs(ploptrdd.out$tau.hat - truth) <= w7/2) # coverage of OPTRDD under partially linear model
#   return(c(c1, c2, c3, c4, c5, c6, c7, w1, w2, w3, w4, w5, w6, w7))
# }
#
# set.seed(123)
#
# system.time({
#   results <- mclapply(1:num_iterations, process_file, mc.cores = num_cores)
#   results <- do.call(rbind, results)
# })
# write.csv(results, "expt_wgan_mats_math.csv", row.names = FALSE)
print(matrix(colMeans(read.csv("expt_wgan_mats_math.csv")), ncol = 2,
             dimnames = list(c("rdrobust conv", "rdrobust biascorr", "rdrobust robust",
                               "rdhonest", "plrd", "AK20 2.113", "AK20 2.181",
                               "smoother optrdd", "pl+optrdd"),
                             c("coverage", "avg width")
             )
), digits = 3
)


# Jacob-Lefgren Reading
#
# dfa <- read.csv("./fake_jl_read/dfa.csv")
# dfb <- read.csv("./fake_jl_read/dfb.csv")
# print(truth <- mean(dfa$y) - mean(dfb$y))
#
# process_file <- function(i) {
#   df <- read.csv(paste0("./fake_jl_read/df_fake_", i, ".csv"))
#   plrd.out <- plrd(df$y, df$x, max.window = 2, verbose = FALSE) # PLRD estimates and CI's.
#   plrd.out2<- plrd(df$y, df$x, max.window = 2, diff.curvatures = TRUE, verbose = FALSE)
#   ploptrdd.out <- ploptrdd(df$y, df$x, max.window = 2, verbose = FALSE)
#   rdrob = ests.from.ci(rdrobust(df$y, df$x, c = 0)$ci) # rdrobust estimates and CI's
#   rdh = RDHonest(df$y ~ df$x)$coefficients # RDHonest estimates and CI's
#   rdrob2 = ests.from.ci(rdrobust(df$y, df$x, c = 0, p = 2)$ci)
#
#   w1 = 2*rdrob[,2] # widths of rdrobust CI's (conventional, bias corrected, robust)
#   w2 = rdh$conf.high-rdh$conf.low # width of RDHonest CI
#   w3 <- 2*plrd.out$tau.plusminus # width of PLRD CI
#   w4 = 2*rdrob2[1,2]*2.113/1.96 # width of AK20 CI with 2.113 instead of 1.96
#   w5 = 2*rdrob[1,2]*2.181/1.96 # width of AK20 CI with 2.181 instead of 1.96
#   w6 = 2*plrd.out2$tau.plusminus # width of OPTRDD with our smoothness
#   w7 = 2*ploptrdd.out$tau.plusminus # width of OPTRDD under partially linear model
#
#   c1 = apply(rdrob, 1, function(x) as.numeric(abs(x[1]-truth) <= x[2])) # coverages of rdrobust CI's (conventional, bias corrected, robust)
#   c2 = as.numeric((rdh$conf.low<=truth)*(rdh$conf.high>=truth)) # coverage of RDHonest CI
#   c3 <- as.numeric(abs(plrd.out$tau.hat - truth) <= w3/2) # coverage of PLRD CI
#   c4 <- as.numeric(abs(rdrob2[1,1] - truth) <= w4/2) # coverage of AK20 CI with 2.113 instead of 1.96
#   c5 <- as.numeric(abs(rdrob[1,1] - truth) <= w5/2) # coverage of AK20 CI with 2.181 instead of 1.96
#   c6 <- as.numeric(abs(plrd.out2$tau.hat - truth) <= w6/2) # coverage of OPTRDD with our smoothness
#   c7 <- as.numeric(abs(ploptrdd.out$tau.hat - truth) <= w7/2) # coverage of OPTRDD under partially linear model
#   return(c(c1, c2, c3, c4, c5, c6, c7, w1, w2, w3, w4, w5, w6, w7))
# }
#
# set.seed(123)
#
# system.time({
#   results <- mclapply(1:num_iterations, process_file, mc.cores = num_cores)
#   results <- do.call(rbind, results)
# })
# write.csv(results, "expt_wgan_jl_read.csv", row.names = FALSE)
print(matrix(colMeans(read.csv("expt_wgan_jl_read.csv")), ncol = 2,
             dimnames = list(c("rdrobust conv", "rdrobust biascorr", "rdrobust robust",
                               "rdhonest", "plrd", "AK20 2.113", "AK20 2.181",
                               "smoother optrdd", "pl+optrdd"),
                             c("coverage", "avg width")
             )
), digits = 3
)


# Jacob-Lefgren Math
#
# dfa <- read.csv("./fake_jl_math/dfa.csv")
# dfb <- read.csv("./fake_jl_math/dfb.csv")
# print(truth <- mean(dfa$y) - mean(dfb$y))
#
# process_file <- function(i) {
#   df <- read.csv(paste0("./fake_jl_math/df_fake_", i, ".csv"))
#   plrd.out <- plrd(df$y, df$x, max.window = 2, verbose = FALSE) # PLRD estimates and CI's.
#   plrd.out2<- plrd(df$y, df$x, max.window = 2, diff.curvatures = TRUE, verbose = FALSE)
#   ploptrdd.out <- ploptrdd(df$y, df$x, max.window = 2, verbose = FALSE)
#   rdrob = ests.from.ci(rdrobust(df$y, df$x, c = 0)$ci) # rdrobust estimates and CI's
#   rdh = RDHonest(df$y ~ df$x)$coefficients # RDHonest estimates and CI's
#   rdrob2 = ests.from.ci(rdrobust(df$y, df$x, c = 0, p = 2)$ci)
#
#   w1 = 2*rdrob[,2] # widths of rdrobust CI's (conventional, bias corrected, robust)
#   w2 = rdh$conf.high-rdh$conf.low # width of RDHonest CI
#   w3 <- 2*plrd.out$tau.plusminus # width of PLRD CI
#   w4 = 2*rdrob2[1,2]*2.113/1.96 # width of AK20 CI with 2.113 instead of 1.96
#   w5 = 2*rdrob[1,2]*2.181/1.96 # width of AK20 CI with 2.181 instead of 1.96
#   w6 = 2*plrd.out2$tau.plusminus # width of OPTRDD with our smoothness
#   w7 = 2*ploptrdd.out$tau.plusminus # width of OPTRDD under partially linear model
#
#   c1 = apply(rdrob, 1, function(x) as.numeric(abs(x[1]-truth) <= x[2])) # coverages of rdrobust CI's (conventional, bias corrected, robust)
#   c2 = as.numeric((rdh$conf.low<=truth)*(rdh$conf.high>=truth)) # coverage of RDHonest CI
#   c3 <- as.numeric(abs(plrd.out$tau.hat - truth) <= w3/2) # coverage of PLRD CI
#   c4 <- as.numeric(abs(rdrob2[1,1] - truth) <= w4/2) # coverage of AK20 CI with 2.113 instead of 1.96
#   c5 <- as.numeric(abs(rdrob[1,1] - truth) <= w5/2) # coverage of AK20 CI with 2.181 instead of 1.96
#   c6 <- as.numeric(abs(plrd.out2$tau.hat - truth) <= w6/2) # coverage of OPTRDD with our smoothness
#   c7 <- as.numeric(abs(ploptrdd.out$tau.hat - truth) <= w7/2) # coverage of OPTRDD under partially linear model
#   return(c(c1, c2, c3, c4, c5, c6, c7, w1, w2, w3, w4, w5, w6, w7))
# }
#
# set.seed(123)
#
# system.time({
#   results <- mclapply(1:num_iterations, process_file, mc.cores = num_cores)
#   results <- do.call(rbind, results)
# })
# write.csv(results, "expt_wgan_jl_math.csv", row.names = FALSE)
print(matrix(colMeans(read.csv("expt_wgan_jl_math.csv")), ncol = 2,
             dimnames = list(c("rdrobust conv", "rdrobust biascorr", "rdrobust robust",
                               "rdhonest", "plrd", "AK20 2.113", "AK20 2.181",
                               "smoother optrdd", "pl+optrdd"),
                             c("coverage", "avg width")
             )
), digits = 3
)

# Oreopoulos data
#
# dfa <- read.csv("./fake_Oreopoulos/dfa.csv")
# dfb <- read.csv("./fake_Oreopoulos/dfb.csv")
# print(truth <- mean(dfa$y) - mean(dfb$y))
#
# process_file <- function(i) {
#   df <- read.csv(paste0("./fake_Oreopoulos/df_fake_", i, ".csv"))
#   plrd.out <- plrd(df$y, df$x, max.window = 10, verbose = FALSE) # PLRD estimates and CI's.
#   plrd.out2<- plrd(df$y, df$x, max.window = 10, diff.curvatures = TRUE, verbose = FALSE)
#   ploptrdd.out <- ploptrdd(df$y, df$x, max.window = 10, verbose = FALSE)
#   rdrob = ests.from.ci(rdrobust(df$y, df$x, c = 0)$ci) # rdrobust estimates and CI's
#   rdh = RDHonest(df$y ~ df$x)$coefficients # RDHonest estimates and CI's
#   rdrob2 = ests.from.ci(rdrobust(df$y, df$x, c = 0, p = 2)$ci)
#
#   w1 = 2*rdrob[,2] # widths of rdrobust CI's (conventional, bias corrected, robust)
#   w2 = rdh$conf.high-rdh$conf.low # width of RDHonest CI
#   w3 <- 2*plrd.out$tau.plusminus # width of PLRD CI
#   w4 = 2*rdrob2[1,2]*2.113/1.96 # width of AK20 CI with 2.113 instead of 1.96
#   w5 = 2*rdrob[1,2]*2.181/1.96 # width of AK20 CI with 2.181 instead of 1.96
#   w6 = 2*plrd.out2$tau.plusminus # width of OPTRDD with our smoothness
#   w7 = 2*ploptrdd.out$tau.plusminus # width of OPTRDD under partially linear model
#
#   c1 = apply(rdrob, 1, function(x) as.numeric(abs(x[1]-truth) <= x[2])) # coverages of rdrobust CI's (conventional, bias corrected, robust)
#   c2 = as.numeric((rdh$conf.low<=truth)*(rdh$conf.high>=truth)) # coverage of RDHonest CI
#   c3 <- as.numeric(abs(plrd.out$tau.hat - truth) <= w3/2) # coverage of PLRD CI
#   c4 <- as.numeric(abs(rdrob2[1,1] - truth) <= w4/2) # coverage of AK20 CI with 2.113 instead of 1.96
#   c5 <- as.numeric(abs(rdrob[1,1] - truth) <= w5/2) # coverage of AK20 CI with 2.181 instead of 1.96
#   c6 <- as.numeric(abs(plrd.out2$tau.hat - truth) <= w6/2) # coverage of OPTRDD with our smoothness
#   c7 <- as.numeric(abs(ploptrdd.out$tau.hat - truth) <= w7/2) # coverage of OPTRDD under partially linear model
#   return(c(c1, c2, c3, c4, c5, c6, c7, w1, w2, w3, w4, w5, w6, w7))
# }
#
# set.seed(123)
#
# system.time({
#   results <- mclapply(1:num_iterations, process_file, mc.cores = num_cores)
#   results <- do.call(rbind, results)
# })
# write.csv(results, "expt_wgan_Oreopoulos.csv", row.names = FALSE)
print(matrix(colMeans(read.csv("expt_wgan_Oreopoulos.csv")), ncol = 2,
             dimnames = list(c("rdrobust conv", "rdrobust biascorr", "rdrobust robust",
                               "rdhonest", "plrd", "AK20 2.113", "AK20 2.181",
                               "smoother optrdd", "pl+optrdd"),
                             c("coverage", "avg width")
             )
), digits = 3
)

# Ludwig-Miller headstart data
#
# dfa <- read.csv("./fake_LudMil/dfa.csv")
# dfb <- read.csv("./fake_LudMil/dfb.csv")
# print(truth <- mean(dfa$y) - mean(dfb$y))
#
# process_file <- function(i) {
#   df <- read.csv(paste0("./fake_LudMil/df_fake_", i, ".csv"))
#   plrd.out <- plrd(df$y, df$x, verbose = FALSE) # PLRD estimates and CI's.
#   plrd.out2<- plrd(df$y, df$x, diff.curvatures = TRUE, verbose = FALSE)
#   ploptrdd.out <- ploptrdd(df$y, df$x, verbose = FALSE)
#   rdrob = ests.from.ci(rdrobust(df$y, df$x, c = 0)$ci) # rdrobust estimates and CI's
#   rdh = RDHonest(df$y ~ df$x)$coefficients # RDHonest estimates and CI's
#   rdrob2 = ests.from.ci(rdrobust(df$y, df$x, c = 0, p = 2)$ci)
#
#   w1 = 2*rdrob[,2] # widths of rdrobust CI's (conventional, bias corrected, robust)
#   w2 = rdh$conf.high-rdh$conf.low # width of RDHonest CI
#   w3 <- 2*plrd.out$tau.plusminus # width of PLRD CI
#   w4 = 2*rdrob2[1,2]*2.113/1.96 # width of AK20 CI with 2.113 instead of 1.96
#   w5 = 2*rdrob[1,2]*2.181/1.96 # width of AK20 CI with 2.181 instead of 1.96
#   w6 = 2*plrd.out2$tau.plusminus # width of OPTRDD with our smoothness
#   w7 = 2*ploptrdd.out$tau.plusminus # width of OPTRDD under partially linear model
#
#   c1 = apply(rdrob, 1, function(x) as.numeric(abs(x[1]-truth) <= x[2])) # coverages of rdrobust CI's (conventional, bias corrected, robust)
#   c2 = as.numeric((rdh$conf.low<=truth)*(rdh$conf.high>=truth)) # coverage of RDHonest CI
#   c3 <- as.numeric(abs(plrd.out$tau.hat - truth) <= w3/2) # coverage of PLRD CI
#   c4 <- as.numeric(abs(rdrob2[1,1] - truth) <= w4/2) # coverage of AK20 CI with 2.113 instead of 1.96
#   c5 <- as.numeric(abs(rdrob[1,1] - truth) <= w5/2) # coverage of AK20 CI with 2.181 instead of 1.96
#   c6 <- as.numeric(abs(plrd.out2$tau.hat - truth) <= w6/2) # coverage of OPTRDD with our smoothness
#   c7 <- as.numeric(abs(ploptrdd.out$tau.hat - truth) <= w7/2) # coverage of OPTRDD under partially linear model
#   return(c(c1, c2, c3, c4, c5, c6, c7, w1, w2, w3, w4, w5, w6, w7))
# }
#
# set.seed(123)
#
# system.time({
#   results <- mclapply(1:num_iterations, process_file, mc.cores = num_cores)
#   results <- do.call(rbind, results)
# })
# write.csv(results, "expt_wgan_LudMil.csv", row.names = FALSE)
print(matrix(colMeans(read.csv("expt_wgan_LudMil.csv")), ncol = 2,
             dimnames = list(c("rdrobust conv", "rdrobust biascorr", "rdrobust robust",
                               "rdhonest", "plrd", "AK20 2.113", "AK20 2.181",
                               "smoother optrdd", "pl+optrdd"),
                             c("coverage", "avg width")
             )
), digits = 4
)
