#' Bias-adjusted Gaussian confidence intervals.
#'
#' @param max.bias Worst-case bias of estimate.
#' @param sampling.se Sampling error of estimate.
#' @param alpha Coverage probability of confidence interval.
#'
#' @return Half-width of confidence interval.
#' @export
get.plusminus = function(max.bias, sampling.se, alpha = 0.95) {
    rel.bias = max.bias/sampling.se
    zz = stats::uniroot(function(z) stats::pnorm(rel.bias - z) +
                          stats::pnorm(-rel.bias - z) + alpha - 1,
                        c(0, rel.bias - stats::qnorm((1 - alpha)/3)))$root
    zz * sampling.se
}
#'
#' Function estimating the Lipschitz constant for the conditional response functions
#'
#' @param y The outcomes.
#' @param x The running variable.
#' @param threshold The threshold for treatment.
#' @param diff.curvatures Whether we consider different curvatures before and after the threshold.
#' @return Lipschitz constant for the conditional response functions.
#' @export
#'
get.Lipschitz.constant <- function(y, x, threshold, 
                                   diff.curvatures = F, alpha.B = 0.05){
  xc = x - threshold; w = xc >= 0
  df = data.frame(y=y, xc=xc, w=w)
  if(!diff.curvatures){
    cubic_reg <- stats::lm(y ~ w * xc + I(xc^2) + I(xc^3), data = df)
    s <- as.vector(summary(cubic_reg)$coeff["I(xc^3)", ])
    B <- 6 * max(abs(s[1] + qnorm(alpha.B/2) * s[2]), 
                 abs(s[1] - qnorm(alpha.B/2) * s[2]))
  } else {
    cubic_left <- stats::lm(y ~ (xc + I(xc^2) + I(xc^3)), 
                            data = df, subset = (xc<0))
    cubic_right <- stats::lm(y ~ (xc + I(xc^2) + I(xc^3)), 
                             data = df, subset = (xc>=0))
    s1 <- as.vector(summary(cubic_left)$coeff["I(xc^3)", ])
    s2 <- as.vector(summary(cubic_right)$coeff["I(xc^3)", ])
    B1 <- 6 * max(abs(s1[1] + qnorm(alpha.B/2) * s1[2]), 
                  abs(s1[1] - qnorm(alpha.B/2) * s1[2]))
    B2 <- 6 * max(abs(s2[1] + qnorm(alpha.B/2) * s2[2]), 
                  abs(s2[1] - qnorm(alpha.B/2) * s2[2]))
    B <- max(B1, B2)
  }
  eps = stats::sd(y)/100
  return(B_hat <- max(B, eps))
}


get.Lipschitz.constant.also.works <- function(y, x, threshold, diff.curvatures = F){
  w = x >= threshold; xc = x - threshold
  if(!diff.curvatures){
    cubic_reg <- stats::lm(y ~ w * xc + I(xc^2) + I(xc^3))
    B <- abs(6 * stats::coef(cubic_reg)["I(xc^3)"])
    quartic_reg <- stats::lm(y ~ w * xc + I(xc^2) + I(xc^3) + I(xc^4))
    B <- max(B, abs(6 * stats::coef(quartic_reg)["I(xc^3)"]))
  } else {
    cubic_reg <- stats::lm(y ~ w * (xc + I(xc^2) + I(xc^3)))
    B_left <- 6 * abs(stats::coef(cubic_reg)["I(xc^3)"])
    B_right <- 6 * abs(stats::coef(cubic_reg)["I(xc^3)"] +
                         stats::coef(cubic_reg)["wTRUE:I(xc^3)"])
    B <- max(B_left, B_right)
    quartic_reg <- stats::lm(y ~ w * (xc + I(xc^2) + I(xc^3) + I(xc^4)))
    if(anova(cubic_reg, quartic_reg)[["Pr(>F)"]][2] <= 0.5){
      B_left <- 6 * abs(stats::coef(quartic_reg)["I(xc^3)"])
      B_right <- 6 * abs(stats::coef(quartic_reg)["I(xc^3)"] +
                             stats::coef(quartic_reg)["wTRUE:I(xc^3)"])
      B_quartic <- max(B_left, B_right)
      B <- sqrt(B * B_quartic)
    }
  }
  eps = stats::sd(y)/100
  return(B_hat <- max(c(B, eps)))
}
#' Getting estimates from rdrobust output
#'
#' rdrobust(...)$ci gives CI's in (lower, upper) format; using which the following
#' function extracts (point estimate, half-width, ci lower, ci upper).
#' @param ci lower and upper bounds of CI's, supports an rdrobust(...)$ci object
#' @return (point estimate, half-width, ci lower, ci upper)
#' @export
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
#' @param x plrd object
#' @export
print.plrd = function(x, ...) {
    if (!is.null(x$tau.hat)) {
        print(paste0(100 * x$alpha, "% CI for tau: ",
                     signif(x$tau.hat, 2), " +/- ", signif(x$tau.plusminus, 2)))
    } else {
        print(paste0(100 * x$alpha, "% CI for tau: [point estimate] +/- ",
                     signif(x$tau.plusminus, 2)))
    }
}
#' @param object plrd object
#' @export
summary.plrd = function(object, ...) {
    unlist(object)[1:5]
}
#' @param x plrd object
#' @param percentage.cumulative.weights The percentage of the cumulative absolute weights user wants to keep (for visualization purposes only)
#' @export
plot.plrd = function(x, percentage.cumulative.weights = .99, ...) {
  full_df = data.frame(X = (x$X-x$threshold),
                       Y = (x$Y - x$tau.hat*x$W),
                       W = x$W,
                       gamma = x$gamma)
  # Finding a window [threshold - l, threshold + l] containing most of the weights
  l = with(full_df, uniroot(function(b) sum(abs(gamma)[abs(X) <= b]) / sum(abs(gamma)) - percentage.cumulative.weights,
                            lower = 0, upper = max(abs(X)))$root)
  df = subset(full_df, abs(full_df$X) <= l)
  xx = seq(from = -l, to = l, length = 200) + x$threshold
  coefs = stats::coef(stats::lm(Y ~ X + I(W*X) + I(X^2),  data = df))
  yy = x$tau.hat*as.numeric(xx>x$threshold) + cbind(1, xx-x$threshold, (xx-x$threshold)*as.numeric(xx>x$threshold), (xx-x$threshold)^2) %*% coefs
    args = list(...)
    if (is.null(dim(x$gamma))) {
        if (!"xlim" %in% names(args)) {
            args$xlim = range(x$X)
        }
        if (!"ylim" %in% names(args)) {
            args$ylim = range(x$Y)
        }
        if (!"xlab" %in% names(args)) {
            args$xlab = "X (running variable)"
        }
        if (!"ylab" %in% names(args)) {
            args$ylab = "Y (response)"
        }
        args$x = NA; args$y = NA
        graphics::layout(matrix(1:2, ncol = 1), heights = c(3, 2))
        graphics::par(mar = c(0, 4.5, 2, 2))
        do.call(graphics::plot, c(args, xaxt = "n"))
        graphics::points(x$X, x$Y,
                         col = c("#CC3311","#009E73")[as.numeric(x$X>x$threshold)+1],
                         cex = .5)
        graphics::lines(xx, yy, col = 'black', lwd = 3)
        graphics::abline(v = x$threshold, lwd = 1.5, lty = 2)
        graphics::abline(v = c(-l,l)+x$threshold, lwd = 1.5, lty = 3)
        graphics::par(mar = c(4.5, 4.5, 0, 2))
        xs0 <- x$gamma.fun.0[[1]]
        ys0 <- x$gamma.fun.0[[2]]
        xs1 <- x$gamma.fun.1[[1]]
        ys1 <- x$gamma.fun.1[[2]]
        plot(
          NA, type = "n",
          xlim = range(xs0, xs1),
          ylim = range(ys0, ys1),
          xlab = args$xlab,
          ylab = expression(hat(gamma)(X))
        )
        if (length (unique(c(xs0, xs1)) > 40)) {
          graphics::points(xs0, ys0, col = "#CC3311", pch = 20, cex = 0.5)
          graphics::points(xs1, ys1, col = "#009E73", pch = 20, cex = 0.5)
        } else {
          graphics::lines(xs0, ys0, col = "#CC3311", lwd = 1)
          graphics::lines(xs1, ys1, col = "#009E73", lwd = 1)
        }
        graphics::abline(v = (c(-l,l)+x$threshold),
                         lwd = 1.5, lty = 3)
        graphics::abline(v = x$threshold, lwd = 1.5, lty = 2)
        graphics::abline(h = 0, lwd = 1.5, lty = 2)
        print(paste("The black curves in the top panel are representative regression functions in our data-driven function class. The dashed line shows the threshold, and the dotted lines indicate the window containing 95% of the cumulative absolute weights. The bottom panel shows two sets of plrd weights because we use cross-fitting."))
    } else {
        stop("Corrupted object.")
    }
}
