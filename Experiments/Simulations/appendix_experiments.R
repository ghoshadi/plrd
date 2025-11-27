rm(list = ls())

if (!requireNamespace("RDHonest", quietly = TRUE)) install.packages("RDHonest")
if (!requireNamespace("rdrobust", quietly = TRUE)) install.packages("rdrobust")
if (!requireNamespace("glmnet", quietly = TRUE)) install.packages("glmnet")
if (!requireNamespace("splines", quietly = TRUE)) install.packages("splines")
if (!requireNamespace("parallel", quietly = TRUE)) install.packages("parallel")
if (!requireNamespace("pbmcapply", quietly = TRUE)) install.packages("pbmcapply")
if (!requireNamespace("R.matlab", quietly = TRUE)) install.packages("R.matlab")

library(RDHonest)
library(rdrobust)
library(glmnet)
library(splines)
library(parallel)
library(pbmcapply)
library(R.matlab)
library(plrd)

source("ploptrdd.R")

sigma.eps <- 0.1295

ests.from.ci <- function(ci){
  if (!is.null(dim(ci))) {
    ci.pm <- t(apply(ci, 1, function(x) c(mean(x), diff(x) / 2, x)))
    colnames(ci.pm) <- c("tau hat", "half width", "CI lower", "CI upper")
    ci.pm
  } else {
    ci.pm <- c(mean(ci), diff(ci) / 2, ci)
    names(ci.pm) <- c("tau hat", "half width", "CI lower", "CI upper")
    ci.pm
  }
}

mu_1 <- function(x){
  if (x < 0) return(0.48 + 1.27 * x + 7.18 * x^2 + 20.21 * x^3 + 21.54 * x^4 + 7.33 * x^5)
  0.52 + 0.84 * x - 3 * x^2 + 7.99 * x^3 - 9.01 * x^4 + 3.56 * x^5
}

mu_2 <- function(x){
  if (x < 0) return(3.71 + 2.30 * x + 3.28 * x^2 + 1.45 * x^3 + 0.23 * x^4 + 0.03 * x^5)
  0.26 + 18.49 * x - 54.81 * x^2 + 74.30 * x^3 - 45.02 * x^4 + 9.83 * x^5
}

mu_3 <- function(x){
  if (x < 0) return(0.48 + 1.27 * x - 0.5 * 7.18 * x^2 +
                      0.7 * 20.21 * x^3 + 1.1 * 21.54 * x^4 + 1.5 * 7.33 * x^5)
  0.52 + 0.84 * x - 0.1 * 3 * x^2 - 0.3 * 7.99 * x^3 - 0.1 * 9.01 * x^4 + 3.56 * x^5
}

mu_4 <- function(x){
  (3 + as.numeric(x > 0)) * x^2
}

single_experiment_pure_noise_appendix <- function(i, n = 500){
  truth <- 0
  set.seed(i)
  x <- runif(n, -1, 1)
  y <- rnorm(n)
  
  plrd.out <- plrd(y, x, 0, diff.curvatures = TRUE, verbose = FALSE)
  ploptrdd.out <- ploptrdd(y, x, 0, verbose = FALSE)
  rdrob <- ests.from.ci(rdrobust(y, x, c = 0, vce = "hc3")$ci)
  rdrob2 <- ests.from.ci(rdrobust(y, x, c = 0, p = 2, vce = "hc3")$ci)
  
  w1 <- 2 * rdrob[1, 2] * 2.181 / 1.96
  w2 <- 2 * rdrob2[1, 2] * 2.113 / 1.96
  w3 <- 2 * plrd.out$tau.plusminus
  w4 <- 2 * ploptrdd.out$tau.plusminus
  
  c1 <- as.numeric(abs(rdrob[1, 1] - truth) <= w1 / 2)
  c2 <- as.numeric(abs(rdrob[1, 1] - truth) <= w2 / 2)
  c3 <- as.numeric(abs(plrd.out$tau.hat - truth) <= w3 / 2)
  c4 <- as.numeric(abs(ploptrdd.out$tau.hat - truth) <= w4 / 2)
  
  c(c1, c2, c3, c4, w1, w2, w3, w4)
}

single_experiment_CCT1_appendix <- function(i, n = 500){
  truth <- mu_1(1e-9) - mu_1(-1e-9)
  set.seed(i)
  x <- 2 * rbeta(n, 2, 4) - 1
  y <- sapply(x, mu_1) + rnorm(n, sd = sigma.eps)
  
  plrd.out <- plrd(y, x, 0, diff.curvatures = TRUE, verbose = FALSE)
  ploptrdd.out <- ploptrdd(y, x, 0, verbose = FALSE)
  rdrob <- ests.from.ci(rdrobust(y, x, c = 0, vce = "hc3")$ci)
  rdrob2 <- ests.from.ci(rdrobust(y, x, c = 0, p = 2, vce = "hc3")$ci)
  
  w1 <- 2 * rdrob[1, 2] * 2.181 / 1.96
  w2 <- 2 * rdrob2[1, 2] * 2.113 / 1.96
  w3 <- 2 * plrd.out$tau.plusminus
  w4 <- 2 * ploptrdd.out$tau.plusminus
  
  c1 <- as.numeric(abs(rdrob[1, 1] - truth) <= w1 / 2)
  c2 <- as.numeric(abs(rdrob[1, 1] - truth) <= w2 / 2)
  c3 <- as.numeric(abs(plrd.out$tau.hat - truth) <= w3 / 2)
  c4 <- as.numeric(abs(ploptrdd.out$tau.hat - truth) <= w4 / 2)
  
  c(c1, c2, c3, c4, w1, w2, w3, w4)
}

single_experiment_CCT2_appendix <- function(i, n = 500){
  truth <- mu_2(1e-9) - mu_2(-1e-9)
  set.seed(i)
  x <- 2 * rbeta(n, 2, 4) - 1
  y <- sapply(x, mu_2) + rnorm(n, sd = sigma.eps)
  
  plrd.out <- plrd(y, x, 0, diff.curvatures = TRUE, verbose = FALSE)
  ploptrdd.out <- ploptrdd(y, x, 0, verbose = FALSE)
  rdrob <- ests.from.ci(rdrobust(y, x, c = 0, vce = "hc3")$ci)
  rdrob2 <- ests.from.ci(rdrobust(y, x, c = 0, p = 2, vce = "hc3")$ci)
  
  w1 <- 2 * rdrob[1, 2] * 2.181 / 1.96
  w2 <- 2 * rdrob2[1, 2] * 2.113 / 1.96
  w3 <- 2 * plrd.out$tau.plusminus
  w4 <- 2 * ploptrdd.out$tau.plusminus
  
  c1 <- as.numeric(abs(rdrob[1, 1] - truth) <= w1 / 2)
  c2 <- as.numeric(abs(rdrob[1, 1] - truth) <= w2 / 2)
  c3 <- as.numeric(abs(plrd.out$tau.hat - truth) <= w3 / 2)
  c4 <- as.numeric(abs(ploptrdd.out$tau.hat - truth) <= w4 / 2)
  
  c(c1, c2, c3, c4, w1, w2, w3, w4)
}

single_experiment_CCT3_appendix <- function(i, n = 500){
  truth <- mu_3(1e-9) - mu_3(-1e-9)
  set.seed(i)
  x <- 2 * rbeta(n, 2, 4) - 1
  y <- sapply(x, mu_3) + rnorm(n, sd = sigma.eps)
  
  plrd.out <- plrd(y, x, 0, diff.curvatures = TRUE, verbose = FALSE)
  ploptrdd.out <- ploptrdd(y, x, 0, verbose = FALSE)
  rdrob <- ests.from.ci(rdrobust(y, x, c = 0, vce = "hc3")$ci)
  rdrob2 <- ests.from.ci(rdrobust(y, x, c = 0, p = 2, vce = "hc3")$ci)
  
  w1 <- 2 * rdrob[1, 2] * 2.181 / 1.96
  w2 <- 2 * rdrob2[1, 2] * 2.113 / 1.96
  w3 <- 2 * plrd.out$tau.plusminus
  w4 <- 2 * ploptrdd.out$tau.plusminus
  
  c1 <- as.numeric(abs(rdrob[1, 1] - truth) <= w1 / 2)
  c2 <- as.numeric(abs(rdrob[1, 1] - truth) <= w2 / 2)
  c3 <- as.numeric(abs(plrd.out$tau.hat - truth) <= w3 / 2)
  c4 <- as.numeric(abs(ploptrdd.out$tau.hat - truth) <= w4 / 2)
  
  c(c1, c2, c3, c4, w1, w2, w3, w4)
}

single_experiment_IK_appendix <- function(i, n = 500){
  truth <- mu_4(1e-9) - mu_4(-1e-9)
  set.seed(i)
  x <- 2 * rbeta(n, 2, 4) - 1
  y <- sapply(x, mu_4) + rnorm(n, sd = sigma.eps)
  
  plrd.out <- plrd(y, x, 0, diff.curvatures = TRUE, verbose = FALSE)
  ploptrdd.out <- ploptrdd(y, x, 0, verbose = FALSE)
  rdrob <- ests.from.ci(rdrobust(y, x, c = 0, vce = "hc3")$ci)
  rdrob2 <- ests.from.ci(rdrobust(y, x, c = 0, p = 2, vce = "hc3")$ci)
  
  w1 <- 2 * rdrob[1, 2] * 2.181 / 1.96
  w2 <- 2 * rdrob2[1, 2] * 2.113 / 1.96
  w3 <- 2 * plrd.out$tau.plusminus
  w4 <- 2 * ploptrdd.out$tau.plusminus
  
  c1 <- as.numeric(abs(rdrob[1, 1] - truth) <= w1 / 2)
  c2 <- as.numeric(abs(rdrob[1, 1] - truth) <= w2 / 2)
  c3 <- as.numeric(abs(plrd.out$tau.hat - truth) <= w3 / 2)
  c4 <- as.numeric(abs(ploptrdd.out$tau.hat - truth) <= w4 / 2)
  
  c(c1, c2, c3, c4, w1, w2, w3, w4)
}

run_appendix_setting <- function(setting,
                                 num_replications = 1000,
                                 n = 500,
                                 num_cores = parallel::detectCores() - 1){
  
  if (setting == "pure_noise") {
    cat("Running appendix pure-noise setting\n")
    results <- pbmclapply(1:num_replications,
                          function(i) single_experiment_pure_noise_appendix(i, n),
                          mc.cores = num_cores)
    results <- do.call(rbind, results)
    print(matrix(colMeans(results), ncol = 2,
                 dimnames = list(paste0("Alt.(", letters[1:4], ")"),
                                 c("coverage", "avg width"))),
          digits = 3)
    write.csv(results, "expt_pure_noise_appendix.csv", row.names = FALSE)
    
  } else if (setting == "CCT1") {
    cat("Running appendix CCT1 setting\n")
    results <- pbmclapply(1:num_replications,
                          function(i) single_experiment_CCT1_appendix(i, n),
                          mc.cores = num_cores)
    results <- do.call(rbind, results)
    print(matrix(colMeans(results), ncol = 2,
                 dimnames = list(paste0("Alt.(", letters[1:4], ")"),
                                 c("coverage", "avg width"))),
          digits = 3)
    write.csv(results, "expt_CCT1_appendix.csv", row.names = FALSE)
    
  } else if (setting == "CCT2") {
    cat("Running appendix CCT2 setting\n")
    results <- pbmclapply(1:num_replications,
                          function(i) single_experiment_CCT2_appendix(i, n),
                          mc.cores = num_cores)
    results <- do.call(rbind, results)
    print(matrix(colMeans(results), ncol = 2,
                 dimnames = list(paste0("Alt.(", letters[1:4], ")"),
                                 c("coverage", "avg width"))),
          digits = 3)
    write.csv(results, "expt_CCT2_appendix.csv", row.names = FALSE)
    
  } else if (setting == "CCT3") {
    cat("Running appendix CCT3 setting\n")
    results <- pbmclapply(1:num_replications,
                          function(i) single_experiment_CCT3_appendix(i, n),
                          mc.cores = num_cores)
    results <- do.call(rbind, results)
    print(matrix(colMeans(results), ncol = 2,
                 dimnames = list(paste0("Alt.(", letters[1:4], ")"),
                                 c("coverage", "avg width"))),
          digits = 3)
    write.csv(results, "expt_CCT3_appendix.csv", row.names = FALSE)
    
  } else if (setting == "IK") {
    cat("Running appendix IK setting\n")
    results <- pbmclapply(1:num_replications,
                          function(i) single_experiment_IK_appendix(i, n),
                          mc.cores = num_cores)
    results <- do.call(rbind, results)
    print(matrix(colMeans(results), ncol = 2,
                 dimnames = list(paste0("Alt.(", letters[1:4], ")"),
                                 c("coverage", "avg width"))),
          digits = 3)
    write.csv(results, "expt_IK_appendix.csv", row.names = FALSE)
    
  } else {
    stop("Unknown setting: use one of 'pure_noise', 'CCT1', 'CCT2', 'CCT3', 'IK'")
  }
}

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    stop("Usage: Rscript appendix_experiments.R <setting> <num_replications> [n]")
  }
  setting <- args[1]
  num_replications <- as.integer(args[2])
  n <- if (length(args) >= 3) as.integer(args[3]) else 500L
  
  n_cores_env <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
  num_cores <- if (!is.na(n_cores_env) && n_cores_env > 0) n_cores_env else parallel::detectCores() - 1L
  
  run_appendix_setting(setting,
                       num_replications = num_replications,
                       n = n,
                       num_cores = num_cores)
}
