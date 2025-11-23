if (!requireNamespace("RDHonest", quietly = TRUE)) install.packages("RDHonest")
if (!requireNamespace("rdrobust", quietly = TRUE)) install.packages("rdrobust")
if (!requireNamespace("glmnet", quietly = TRUE)) install.packages("glmnet")
if (!requireNamespace("splines", quietly = TRUE)) install.packages("splines")
if (!requireNamespace("parallel", quietly = TRUE)) install.packages("parallel")
if (!requireNamespace("R.matlab", quietly = TRUE)) install.packages("R.matlab")
if (!requireNamespace("pbmcapply", quietly = TRUE)) install.packages("pbmcapply")
if (!requireNamespace("arrow", quietly = TRUE)) install.packages("arrow")
if (!requireNamespace("plrd", quietly = TRUE)) remotes::install_github("ghoshadi/plrd")

# If having trouble loading the package arrow, try nanoparquet

library(RDHonest) # version 1.0.1
library(rdrobust) # version 3.0.0
library(glmnet) # version 4.1.10
library(splines) # version 4.5.1
library(parallel) # version 4.5.1
library(R.matlab) # version 3.7.0
library(pbmcapply) # version 1.5.1
library(arrow) # version 22.0.0
library(plrd) # our package, version 0.0.1

num_cores = parallel::detectCores()

# helper function to compute estimates and standard errors from rdrobust summary
ests.from.ci <- function(ci){
  if(!is.null(dim(ci))){
    ci.pm <- t(apply(ci, 1, function(x) c(mean(x), diff(x)/2, x)))
    colnames(ci.pm) = c("tau hat", "half width", "CI lower", "CI upper")
    return(ci.pm)
  } else{
    ci.pm <- c(mean(ci), diff(ci)/2, ci)
    names(ci.pm) = c("tau hat", "half width", "CI lower", "CI upper")
    return(ci.pm)
  }
}

num_reps= 20
data.names = c("Lee", "senate", "LudMill",  "Meyersson",
               "Mats_read", "Mats_math", "JL_read", "JL_math", "Oreopoulos")

results_list = list()
for(d in seq_along(data.names)){
  print(gt <- read.csv("ground_truths.csv")[d, ])
  truth = gt[2]; dname = gt[1]
  set.seed(42)
  results <- do.call(rbind,
                     pbmclapply(1:num_reps,
                                function(i) {
                                  df <- read_parquet(paste0("./fake_datasets/wgan_",dname,"/df_", i, ".parquet"))
                                  plrd.out <- plrd(df$y, df$x, threshold = 0) # PLRD estimates and CI's.
                                  rdrob = ests.from.ci(rdrobust(df$y, df$x, c = 0, vce = "hc3")$ci) # rdrobust estimates and CI's
                                  rdh = RDHonest(df$y ~ df$x)$coefficients # RDHonest estimates and CI's

                                  w1 = 2*rdrob[,2] # widths of rdrobust CI's (conventional, bias corrected, robust)
                                  w2 = rdh$conf.high-rdh$conf.low # width of RDHonest CI
                                  w3 <- 2*plrd.out$tau.plusminus # width of PLRD CI

                                  c1 = apply(rdrob, 1, function(x) as.numeric(abs(x[1]-truth) <= x[2])) # coverages of rdrobust CI's (conventional, bias corrected, robust)
                                  c2 = as.numeric((rdh$conf.low<=truth)*(rdh$conf.high>=truth)) # coverage of RDHonest CI
                                  c3 <- as.numeric(abs(plrd.out$tau.hat - truth) <= w3/2) # coverage of PLRD CI
                                  return(c(c1, c2, c3, w1, w2, w3))
                                },
                                mc.cores = num_cores)
                     )
  colnames(results) = c("Conventional_coverage", "Bias.Corrected_coverage", "rdrobust_coverage",
                        "RDHonest_coverage", "plrd_coverage",
                        "Conventional_width", "Bias.Corrected_width", "rdrobust_width",
                        "RDHonest_width", "plrd_width")
  write.csv(results, paste0("./evaluations/expt_wgan_",dname,".csv"), row.names = FALSE)
}

wgan_table = c()
for(d in seq_along(data.names)){
  print(gt <- read.csv("ground_truths.csv")[d, ])
  dname = gt[1]
  res = read.csv(paste0("./evaluations/expt_wgan_",dname,".csv"))
  colnames(res) = c("Conventional_coverage", "Bias.Corrected_coverage", "rdrobust_coverage",
                    "RDHonest_coverage", "plrd_coverage",
                    "Conventional_width", "Bias.Corrected_width", "rdrobust_width",
                    "RDHonest_width", "plrd_width")
  wgan_table = rbind(wgan_table, colMeans(res))
}
rownames(wgan_table) = data.names
print(wgan_table[,c(1,6,2,7,3,8,4,9,5,10)], digits = 3)
