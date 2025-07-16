#' Bias-adjusted Gaussian confidence intervals.
#'
#' @param max.bias Worst-case bias of estimate.
#' @param sampling.se Sampling error of estimate.
#' @param alpha Coverage probability of confidence interval.
#'
#' @return Half-width of confidence interval.
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
#' @param alpha.B Significance threshold we use in estimation of the Lipschitz constant.
#' @return Lipschitz constant for the conditional response functions.
#'
get.Lipschitz.constant <- function(y, x, threshold,
                                   diff.curvatures = F, alpha.B = 0.05){
  xc = x - threshold; w = xc >= 0
  df = data.frame(y=y, xc=xc, w=w)
  if(!diff.curvatures){
    cubic_reg <- stats::lm(y ~ w * xc + I(xc^2) + I(xc^3), data = df)
    s <- as.vector(summary(cubic_reg)$coeff["I(xc^3)", ])
    B <- 6 * max(abs(s[1] + stats::qnorm(alpha.B/2) * s[2]),
                 abs(s[1] - stats::qnorm(alpha.B/2) * s[2]))
  } else {
    cubic_left <- stats::lm(y ~ (xc + I(xc^2) + I(xc^3)),
                            data = df, subset = (xc<0))
    cubic_right <- stats::lm(y ~ (xc + I(xc^2) + I(xc^3)),
                             data = df, subset = (xc>=0))
    s1 <- as.vector(summary(cubic_left)$coeff["I(xc^3)", ])
    s2 <- as.vector(summary(cubic_right)$coeff["I(xc^3)", ])
    B1 <- 6 * max(abs(s1[1] + stats::qnorm(alpha.B/2) * s1[2]),
                  abs(s1[1] - stats::qnorm(alpha.B/2) * s1[2]))
    B2 <- 6 * max(abs(s2[1] + stats::qnorm(alpha.B/2) * s2[2]),
                  abs(s2[1] - stats::qnorm(alpha.B/2) * s2[2]))
    B <- max(B1, B2)
  }
  eps = stats::sd(y)/100
  return(B_hat <- max(B, eps))
}

#' Print a plrd object
#' @param x plrd object
#' @param digits number of digits to print
#' @param ... Additional arguments (currently ignored).
#' @export
print.plrd = function(x, digits = NULL, ...) {
  if(is.null(digits)) digits = 3
  cat(paste0(100 * x$alpha, "% CI for tau: [", signif(x$ci.lower, digits), ", ", signif(x$ci.upper, digits),"]\n"))
  # if (!is.null(x$tau.hat)) {
  #   print(paste0(100 * x$alpha, "% CI for tau: ",
  #                signif(x$tau.hat, 2), " +/- ", signif(x$tau.plusminus, 2)))
  # } else {
  #   print(paste0(100 * x$alpha, "% CI for tau: [point estimate] +/- ",
  #                signif(x$tau.plusminus, 2)))
  # }
}

#' plrd summary
#' @param object plrd object
#' @param ... Additional arguments (currently ignored).
#' @export
summary.plrd = function(object, ...) {
  unlist(object)[1:6]
}

#' Plot a plrd object
#' @param x plrd object
#' @param percentage.cumulative.weights The percentage of the cumulative absolute weights user wants to keep (for visualization purposes only)
#' @param ... Additional arguments (currently ignored).
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

#' Compute MSE-optimal Imbens-Kalyanaraman bandwidth for a sharp RD.
#'
#' This convenience function computes weights using the Imbens-Kalyanaraman bandwidth procedure.
#' The code does exactly what the MATLAB code available on the author's website does.
#'
#' @param Y The outcomes.
#' @param X The running variable.
#' @param threshold The threshold.
#' @param kernel The kernel type used to construct weights within the bandwidth.
#'
#' @return A list containing the sample weights along with optimal bandwidth.
#'
#' @references Imbens, G., and Kalyanaraman, K. (2012).
#'  Optimal Bandwidth Choice for the Regression Discontinuity Estimator.
#'  The Review of Economic Studies, 79(3).
#'
#' @examples
#' \donttest{
#' n = 1000; threshold = 0
#' X = runif(n, -1, 1)
#' W = as.numeric(X >= threshold)
#' Y = (1 + 2*W)*(1 + X^2) + 1 / (1 + exp(X)) + rnorm(n, sd = .5)
#' IK_bandwidth(Y, X, threshold)
#' }
#'
#' @export
IK_bandwidth <- function(Y, X, threshold, kernel = c("triangular", "uniform")) {
  if (threshold >= max(X) || threshold <= min(X)) {
    stop("RD threshold is outside the running variable range.")
  }
  kernel <- match.arg(kernel)

  x <- X - threshold
  n <- length(Y)

  h.silverman <- 1.84 * stats::sd(x) * n^(-1/5)

  i.plus.1 <- x >= 0 & x <= h.silverman
  i.min.1  <- x < 0 & x >= -h.silverman

  y.ave.plus.1 <- mean(Y[i.plus.1])
  y.ave.min.1  <- mean(Y[i.min.1])

  dy.plus.1 <- Y[i.plus.1] - y.ave.plus.1
  dy.min.1  <- Y[i.min.1]  - y.ave.min.1

  sigmas <- (sum(dy.plus.1^2) + sum(dy.min.1^2)) / (sum(i.plus.1) + sum(i.min.1))
  fc <- (sum(i.plus.1) + sum(i.min.1)) / (2 * n * h.silverman)

  # third derivative
  median.plus <- stats::median(x[x >= 0])
  median.min  <- stats::median(x[x < 0])

  middle <- x >= median.min & x <= median.plus
  x.mid <- x[middle]
  y.mid <- Y[middle]

  tt <- cbind(1, x.mid >= 0, x.mid, x.mid^2, x.mid^3)
  gamma <- solve(t(tt) %*% tt, t(tt) %*% y.mid)
  third.der <- 6 * gamma[5]

  # second derivatives
  h.plus.2 <- 3.56 * (sigmas / (fc * max(third.der^2, 0.01)))^(1/7) * sum(x >= 0)^(-1/7)
  h.min.2  <- 3.56 * (sigmas / (fc * max(third.der^2, 0.01)))^(1/7) * sum(x < 0)^(-1/7)

  i.min.3  <- x < 0 & x >= -h.min.2
  i.plus.3 <- x >= 0 & x <= h.plus.2
  if (sum(i.min.1) <= 2 || sum(i.plus.3) <= 2) {
    stop("Insufficient observations near discontinuity.")
  }

  # left second derivative
  x.min <- x[i.min.3]
  y.min <- Y[i.min.3]
  t.min <- cbind(1, x.min, x.min^2)
  beta.min <- solve(t(t.min) %*% t.min, t(t.min) %*% y.min)
  second.der.min <- 2 * beta.min[3]

  # right second derivative
  x.plus <- x[i.plus.3]
  y.plus <- Y[i.plus.3]
  t.plus <- cbind(1, x.plus, x.plus^2)
  beta.plus <- solve(t(t.plus) %*% t.plus, t(t.plus) %*% y.plus)
  second.der.plus <- 2 * beta.plus[3]

  # regularization terms
  r.plus <- 2160 * sigmas / (sum(i.plus.3) * h.plus.2^4)
  r.min  <- 2160 * sigmas / (sum(i.min.3)  * h.min.2^4)

  CK <- 3.4375
  denom <- (second.der.plus - second.der.min)^2 + r.plus + r.min
  h.opt <- CK * ((2 * sigmas) / (fc * denom))^(1/5) * n^(-1/5)

  # Compute weights
  dist <- abs((X - threshold) / h.opt)
  if (kernel == "triangular") {
    weights <- (1 - dist) * (dist <= 1) / h.opt
    weights <- weights / sum(weights)
  } else {
    weights <- (1 - dist) * (dist <= 1)
    weights[weights > 0 ] <- 1
  }

  list(
    weights = weights,
    bandwidth = h.opt
  )
}
