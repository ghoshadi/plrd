rm(list = ls())
if (!requireNamespace("RDHonest", quietly = TRUE)) install.packages("RDHonest")
if (!requireNamespace("rdrobust", quietly = TRUE)) install.packages("rdrobust")
if (!requireNamespace("glmnet", quietly = TRUE)) install.packages("glmnet")
if (!requireNamespace("splines", quietly = TRUE)) install.packages("splines")
if (!requireNamespace("parallel", quietly = TRUE)) install.packages("parallel")
if (!requireNamespace("R.matlab", quietly = TRUE)) install.packages("R.matlab")
if (!requireNamespace("readstata13", quietly = TRUE)) install.packages("readstata13")

library(RDHonest) # version 1.0.0
library(rdrobust) # version 2.2
library(glmnet) # version 4.1.8
library(splines) # version 4.4.2
library(parallel) # version 4.4.2
library(R.matlab) # version 3.7.0
library(plrd) # our package
library(readstata13)

mred = "#CC3311"; mgreen = "#009E73"

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

#=======================================================
# Lee_data
#=======================================================

set.seed(123)
Lee_data <- data.frame(readMat("davidlee.mat"))
y = Lee_data$y1; x = Lee_data$x; cutoff = 0
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))
plot(x, y, pch = 16, cex = 0.5,
     col = c(mred, mgreen)[1+as.numeric(x>cutoff)],
     xlab = "Difference in vote share (Democratic - Republican) in last election",
     ylab = "Democratic vote share",
     main = "Lee (2008) data (incumbency advantage in elections)")
abline(v = 0, col = "black", lty = 1, lwd = 1.5)

rdrob.out <- ests.from.ci(rdrobust(y, x, c = cutoff)$ci)
ak = RDHonest(y ~ x, cutoff = cutoff)
ak.hw = (ak$coefficients$conf.high - ak$coefficients$conf.low)/2
ak.out = c(ak$coefficients$estimate, ak.hw, ak$coefficients$conf.low, ak$coefficients$conf.high)
out <- plrd(y, x, cutoff)
plot(out)
plrd.out <- c(out$tau.hat, out$tau.plusminus, out$tau.hat-out$tau.plusminus, out$tau.hat+out$tau.plusminus)
print(rbind(rdrob.out, ak.out, plrd.out), digits = 3)

#=======================================================
# Auxiliary function for applying different RD methods
#=======================================================

applying_diff_RD_methods <- function(y, x, cutoff){
  par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))
  plot(x, y, pch = 16, cex = 0.5,
       col = c(mred, mgreen)[1+as.numeric(x>cutoff)])
  abline(v = cutoff, col = "black", lty = 1, lwd = 1.5)

  rdrob.out <- ests.from.ci(rdrobust(y, x, c = cutoff)$ci)
  ak = RDHonest(y ~ x, cutoff = cutoff)
  ak.hw = (ak$coefficients$conf.high - ak$coefficients$conf.low)/2
  ak.out = c(ak$coefficients$estimate, ak.hw, ak$coefficients$conf.low, ak$coefficients$conf.high)
  out <- plrd(y, x, cutoff); plot(out)
  plrd.out <- c(out$tau.hat, out$tau.plusminus, out$tau.hat-out$tau.plusminus, out$tau.hat+out$tau.plusminus)
  print(rbind(rdrob.out, ak.out, plrd.out), digits = 4)
}

#=======================================================
# Meyersson (2014) data, downloaded from https://github.com/Mixtape-Sessions/Regression-Discontinuity/tree/main
#=======================================================

polecon <- read.dta13("./polecon.dta")
# X: win margin for the Islamic party relative to the largest non-Islamic party
# the cutoff is therefore c = 0. The municipalities that fall below the cutoff,
# the control group, receive a secular mayor. Those above the cutoff, the
# treatment group, receive an Islamic mayor.
# The outcome of interest (variable Y in the dataset) is the school attainment for women
# who were (potentially) in high school during the period 1994-2000, measured with variables
# extracted from the 2000 census. The particular outcome we re-analyze is the share of the
# cohort of women ages 15 to 20 in 2000 who had completed high school by 2000.
applying_diff_RD_methods(y = polecon$Y, x = polecon$X, cutoff = 0)
# The results indicate that in municipalities where the Islamic party barely won,
# the educational attainment of women is roughly 3 percentage points higher than
# in municipalities where the party barely lost

#=======================================================
# Matsudaira data
#=======================================================

# Matsudaira Reading
mats_read <- read.csv("./mats_read.csv")
applying_diff_RD_methods(y = mats_read$y, x = - mats_read$x, cutoff = 0)

# Matsudaira Math
mats_math <- read.csv("./mats_math.csv")
applying_diff_RD_methods(y = mats_math$y, x = - mats_math$x, cutoff = 0)

#=======================================================
# Jacob-Lefgren data
#=======================================================

# Jacob-Lefgren Reading
jl_read <- read.csv("./jl_read.csv")
applying_diff_RD_methods(y = jl_read$y, x = - jl_read$x, cutoff = 0)

# Jacob-Lefgren Math
jl_math <- read.csv("./jl_math.csv")
applying_diff_RD_methods(y = jl_math$y, x = - jl_math$x, cutoff = 0)

#=======================================================
# Oreopoulos data
#=======================================================

Oreopoulos_data <- read.csv("./uk_analysis_sample.csv")
applying_diff_RD_methods(y = Oreopoulos_data$learn, x = Oreopoulos_data$x, cutoff = 46.99)

#=======================================================
# Ludwig-Miller headstart data
#=======================================================

LMdata = read.csv("./headstart.csv")
LMdata = LMdata[!is.na(LMdata$mort_age59_related_postHS),]
LMdata = LMdata[!is.na(LMdata$povrate60),]
applying_diff_RD_methods(y = LMdata$mort_age59_related_postHS, x = LMdata$povrate60, cutoff = 59.1984)



#=======================================================
# senate data
#=======================================================

senate_data <- read.csv("./senate.csv")
applying_diff_RD_methods(y = senate_data$y, x = senate_data$x, cutoff = 0)
