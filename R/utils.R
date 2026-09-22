#' Bias-adjusted Gaussian confidence intervals.
#'
#' @param max.bias Worst-case bias of estimate.
#' @param sampling.se Sampling error of estimate.
#' @param alpha Coverage probability of confidence interval.
#' @return Half-width of confidence interval.
#' @keywords internal
#' @noRd
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
#' @keywords internal
#' @noRd
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

#' Extract plrd coefficient
#' @param object plrd object
#' @param ... Additional arguments (currently ignored).
#' @export
coef.plrd = function(object, ...) {
  c(estimate = object$tau.hat)
}

#' Print a plrd object
#' @param x plrd object
#' @param digits number of digits to print
#' @param percentage.cumulative.weights The percentage of the cumulative absolute weights to determine effective sample size
#' @param ... Additional arguments passed to print methods.
#' @export
print.plrd = function(x, digits = max(3, getOption("digits") - 3), percentage.cumulative.weights = 0.99, ...) {
  cat(paste0("Partially linear regression discontinuity inference: \n"))
  cat(paste0("Threshold: ", signif(x$threshold, digits), "\n"))
  cat(paste0("Lipschitz constant: ", signif(x$Lipschitz.constant, digits), "\n"))
  cat(paste0("Max bias: ", signif(x$max.bias, digits), "\n"))
  cat(paste0("Sampling SE: ", signif(x$sampling.se, digits), "\n"))
  cat(paste0("Effective sample size: ", signif(length(x$gamma[abs(x$gamma)<=stats::quantile(abs(x$gamma),percentage.cumulative.weights)]), digits), "\n"))
  cat(paste0("Confidence level: ", x$alpha * 100, "%", "\n"))
  cat("\n")
  print(summary(x), digits = digits, ...)
  invisible(x)
}

#' plrd summary
#' @param object plrd object
#' @param ... Additional arguments (currently ignored).
#' @export
summary.plrd = function(object, ...) {
  out = data.frame(
    Estimate = object$tau.hat,
    `CI Lower` = object$ci.lower,
    `CI Upper` = object$ci.upper,
    `p-value` = object$pval,
    check.names = FALSE
  )
  rownames(out) = "RD Estimate"

  out
}

#' Fit a natural spline under a Lipschitz constraint
#'
#' @param y The outcomes.
#' @param x The centered running variable.
#' @param w Indicator for observations above the threshold.
#' @param spline.df Degrees of freedom of the natural spline.
#' @param B Lipschitz constant for the second derivative.
#' @param diff.curvatures Whether we consider different curvatures before and after the threshold.
#' @return A constrained spline fit and its prediction function.
#' @keywords internal
#' @noRd
fit_constrained_spline = function(y, x, w, spline.df, B, diff.curvatures) {
  if (spline.df < 2) stop("'spline.df' must be at least 2.")

  get_G = function(S) {
    knots = attr(S, "knots")
    bk = attr(S, "Boundary.knots")
    ak = sort(c(rep(bk, 4), knots))

    const = splines::splineDesign(
      ak, bk, ord = 4, derivs = c(2, 2)
    )[, -1, drop = FALSE]
    qrc = qr(t(const))

    breaks = c(bk[1], knots, bk[2])
    mid = (breaks[-length(breaks)] + breaks[-1]) / 2
    D3 = splines::splineDesign(
      ak, mid, ord = 4, derivs = rep(3, length(mid))
    )[, -1, drop = FALSE]

    list(
      G = t(qr.qty(qrc, t(D3)))[, -(1:2), drop = FALSE],
      breaks = breaks
    )
  }

  if (isTRUE(diff.curvatures)) {
    i0 = w == 0
    i1 = w == 1
    bk = range(x)

    # Separate spline bases
    S0 = splines::ns(x[i0], df = spline.df, Boundary.knots = bk)
    S1 = splines::ns(x[i1], df = spline.df, Boundary.knots = bk)
    p = spline.df

    s00 = drop(stats::predict(S0, 0))
    s10 = drop(stats::predict(S1, 0))

    Z0 = matrix(0, length(x), p)
    Z1 = matrix(0, length(x), p)
    Z0[i0, ] = sweep(S0, 2, s00)
    Z1[i1, ] = sweep(S1, 2, s10)
    Z = cbind(1, Z0, Z1)

    g0 = get_G(S0)
    g1 = get_G(S1)

    G0 = g0$G[g0$breaks[-length(g0$breaks)] < 0, , drop = FALSE]
    G1 = g1$G[g1$breaks[-1] > 0, , drop = FALSE]

    A = matrix(0, nrow(G0) + nrow(G1), ncol(Z))
    A[seq_len(nrow(G0)), 1 + seq_len(p)] = G0
    A[nrow(G0) + seq_len(nrow(G1)), 1 + p + seq_len(p)] = G1

  } else {
    S = splines::ns(x, df = spline.df)
    p = ncol(S)
    Z = cbind(1, S, w*x)

    g = get_G(S)
    A = matrix(0, nrow(g$G), ncol(Z))
    A[, 1 + seq_len(p)] = g$G
  }

  if (qr(Z)$rank < ncol(Z))
    stop("'spline.df' is too large for the available running-variable support.")

  beta = quadprog::solve.QP(
    Dmat = crossprod(Z),
    dvec = drop(crossprod(Z, y)),
    Amat = t(rbind(A, -A)),
    bvec = rep(-B, 2*nrow(A))
  )$solution

  predict.fun = function(x.new, w.new) {
    if (isTRUE(diff.curvatures)) {
      Z.new = cbind(
        1,
        (1-w.new) * sweep(stats::predict(S0, x.new), 2, s00),
        w.new * sweep(stats::predict(S1, x.new), 2, s10)
      )
    } else {
      Z.new = cbind(1, stats::predict(S, x.new), w.new*x.new)
    }
    drop(Z.new %*% beta)
  }

  list(coefficients = beta, predict = predict.fun)
}
#' Find the effective weight window
#'
#' Computes per-side distances from the threshold spanning the cumulative
#' absolute plrd weight mass up to (but not exceeding) the target share.
#' Treats the two sides independently which results in an asymmetric window.
#'
#' @param x plrd object
#' @param percentage.cumulative.weights Share of cumulative absolute weights to retain on each side.
#' @return A list with per-side distances (\code{l_below}, \code{l_above}) from the threshold.
#' @keywords internal
find_weight_window <- function(x, percentage.cumulative.weights = 0.99) {
  if (percentage.cumulative.weights <= 0 || percentage.cumulative.weights > 1) {
    stop("`percentage.cumulative.weights` must be in (0, 1].")
    }
  xs0 <- x$gamma.fun.0[[1]]; gs0 <- x$gamma.fun.0[[2]]  # below
  xs1 <- x$gamma.fun.1[[1]]; gs1 <- x$gamma.fun.1[[2]]  # above

  # Find farthest point still below target percentage (else nearest)
  weight_quantile_distance <- function(xs, gs) {
    ord <- order(abs(xs - x$threshold))
    cum <- cumsum(abs(gs[ord])) / sum(abs(gs))
    idx <- which(cum < percentage.cumulative.weights)
    i   <- if (length(idx)) idx[length(idx)] else 1
    abs(xs[ord][i] - x$threshold)
  }

  l_below <- weight_quantile_distance(xs0, gs0)
  l_above <- weight_quantile_distance(xs1, gs1)

  list(
    l_below = l_below,
    l_above = l_above
  )
}

#' Plot a plrd object
#'
#' @param x plrd object
#' @param type Type of plot the user wants to see. We offer three options: "default", "weights", and "combined". The "default" option shows the scatterplot of the original data with black curves that are representative regression functions in our data-driven function class. The dashed line shows the threshold, and the dotted lines indicate the window containing a percentage (99% by default) of the cumulative absolute plrd weights. The "weights" option plots plrd weights (note, two sets of plrd weights because we use cross-fitting). The "combined" option is a fancy plot combining both the scatterplot and the plot of plrd weights.
#' @param percentage.cumulative.weights The percentage of the cumulative absolute weights user wants to keep (for visualization purposes only)
#' @param spline.df Degrees of freedom of the natural spline used to plot a representative member of the data-driven function class, constrained to satisfy the smoothness condition.
#' @param ... Additional graphical arguments to customize the plot, such as \code{xlim}, \code{ylab}, \code{main}, etc.
#' @export
plot.plrd = function(x, type = "default", percentage.cumulative.weights = .99, spline.df = 3, ...) {
  op <- graphics::par(no.readonly = TRUE)
  threshold <- x$threshold
  ge.threshold <- as.numeric(x$X >= threshold) # Introduce to avoid any issues for implementation of fuzzy RD.
  full_df = data.frame(Xc = (x$X - threshold),
                       Y0 = (x$Y - x$tau.hat * ge.threshold),
                       ge.threshold  = ge.threshold)

  # Fit spline satisfying the smoothness restriction used by plrd
  fit = fit_constrained_spline(
    y = full_df$Y0,
    x = full_df$Xc,
    w = full_df$ge.threshold,
    spline.df = spline.df,
    B = x$Lipschitz.constant,
    diff.curvatures = x$diff.curvatures
  )

  # Plot model only in a window [threshold - l below, threshold + l above] containing the specified cumulative absolute weights (effective support)
  windows.effective.support  = find_weight_window(x, percentage.cumulative.weights)
  x_lo = threshold - windows.effective.support$l_below
  x_hi = threshold + windows.effective.support$l_above

  # Generate grid on running variable for x coordinates with corresponding y coordinates within effective support
  step = min(threshold - x_lo, x_hi - threshold) / 200
  xx_left  = seq(x_lo, threshold, by = step)
  xx_right = seq(threshold, x_hi, by = step)

  yy_left  = fit$predict(xx_left - threshold, 0)
  yy_right = fit$predict(xx_right - threshold, 1) + x$tau.hat

  # Extract weight coordinates
  xs0 <- x$gamma.fun.0[[1]]
  ys0 <- x$gamma.fun.0[[2]]
  xs1 <- x$gamma.fun.1[[1]]
  ys1 <- x$gamma.fun.1[[2]]

  args = list(...)
  if (is.null(dim(x$gamma))) {
    if (!"xlim" %in% names(args)) args$xlim = if (type == "weights") range(xs0, xs1) else range(x$X)
    if (!"ylim" %in% names(args)) args$ylim = if (type == "weights") range(ys0, ys1) else range(x$Y)
    if (!"xlab" %in% names(args)) args$xlab = "X (running variable)"
    if (!"ylab" %in% names(args)) args$ylab = if (type == "weights") expression("plrd weights " ~ hat(gamma)(X)) else "Y (response)"
    args$x = NA; args$y = NA
    if(type == "default"){
      graphics::layout(matrix(1))
      graphics::par(mar = c(4.5, 4.5, 2, 2))
      do.call(graphics::plot, args)
      graphics::points(x$X, x$Y,
                       col = c("#CC3311","#009E73")[as.numeric(x$X >= threshold)+1],
                       cex = .5)
      graphics::lines(xx_left,  yy_left,  col = 'black', lwd = 3)
      graphics::lines(xx_right, yy_right, col = 'black', lwd = 3)
      graphics::abline(v = threshold, lwd = 1.5, lty = 2)
      graphics::abline(v = c(x_lo, x_hi), lwd = 1.5, lty = 3)
    } else if (type == "weights"){
      graphics::layout(matrix(1))
      graphics::par(mar = c(4.5, 4.5, 2, 2))
      do.call(graphics::plot, utils::modifyList(args, list(type = "n")))
      if (length (unique(c(xs0, xs1))) > 40) {
        graphics::points(xs0, ys0, col = "#CC3311", pch = 20, cex = 0.5)
        graphics::points(xs1, ys1, col = "#009E73", pch = 20, cex = 0.5)
      } else {
        graphics::lines(xs0, ys0, col = "#CC3311", lwd = 1)
        graphics::lines(xs1, ys1, col = "#009E73", lwd = 1)
      }
      graphics::abline(v = c(x_lo, x_hi),
                       lwd = 1.5, lty = 3)
      graphics::abline(v = threshold, lwd = 1.5, lty = 2)
      graphics::abline(h = 0, lwd = 1.5, lty = 2)
    } else if (type == "combined"){
      graphics::layout(matrix(1:2, ncol = 1), heights = c(4, 3.5))
      graphics::par(mar = c(0, 4.5, 2, 2))
      do.call(graphics::plot, utils::modifyList(args, list(xaxt = "n", yaxt = "n")))
      graphics::axis(2, las = 1)
      graphics::points(x$X, x$Y,
                       col = c("#CC3311","#009E73")[as.numeric(x$X >= threshold)+1],
                       cex = .5)
      graphics::lines(xx_left,  yy_left,  col = 'black', lwd = 3)
      graphics::lines(xx_right, yy_right, col = 'black', lwd = 3)
      graphics::abline(v = threshold, lwd = 1.5, lty = 2)
      graphics::abline(v = c(x_lo, x_hi), lwd = 1.5, lty = 3)
      graphics::par(mar = c(4.5, 4.5, 0, 2))
      plot(
        NA, type = "n",
        xlim = args$xlim,
        ylim = range(ys0, ys1),
        xlab = args$xlab,
        ylab = expression(hat(gamma)(X)),
        yaxt="n"
      )
      graphics::axis(2, at = pretty(range(ys0, ys1), n = 4),
                     las = 1, cex.axis = 0.9)
      if (length (unique(c(xs0, xs1))) > 40) {
        graphics::points(xs0, ys0, col = "#CC3311", pch = 20, cex = 0.5)
        graphics::points(xs1, ys1, col = "#009E73", pch = 20, cex = 0.5)
      } else {
        graphics::lines(xs0, ys0, col = "#CC3311", lwd = 1)
        graphics::lines(xs1, ys1, col = "#009E73", lwd = 1)
      }
      graphics::abline(v = c(x_lo, x_hi),
                       lwd = 1.5, lty = 3)
      graphics::abline(v = threshold, lwd = 1.5, lty = 2)
      graphics::abline(h = 0, lwd = 1.5, lty = 2)
    } else {
      stop("Please select plot type among 'default', 'weights', or 'combined'.")
    }
    on.exit(graphics::par(op), add = TRUE) # restore par
    on.exit(graphics::layout(1), add = TRUE) # reset layout to 1 panel
  } else {
    stop("Corrupted object.")
  }
  invisible(list(
    main = list(
      x_coordinates = c(xx_left, xx_right),
      y_coordinates = c(yy_left, yy_right)),
    gamma = list(
      x_coordinates = c(xs0, xs1),
      y_coordinates = c(ys0, ys1))
  ))
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
#' set.seed(42)
#' n = 1000; threshold = 0
#' X = runif(n, -1, 1)
#' W = as.numeric(X >= threshold)
#' Y = (1 + 2*W)*(1 + X^2) + 1 / (1 + exp(X)) + rnorm(n, sd = .5)
#' out = IK_bandwidth(Y, X, threshold)
#'
#' @export
IK_bandwidth <- function(Y, X, threshold,
                         kernel = c("triangular", "uniform", "epanechnikov")) {
  if (length(Y) != length(X)) stop("'Y' and 'X' must have the same length.")
  if (threshold >= max(X) || threshold <= min(X))
    stop("RD threshold is outside the running variable range.")

  kernel <- match.arg(kernel)
  x <- X - threshold; n <- length(x)
  left <- x < 0; right <- !left

  # Density and conditional variances at the threshold
  h1 <- 1.84 * stats::sd(x) * n^(-1/5)
  i.min <- left & x >= -h1
  i.plus <- right & x <= h1
  n1 <- c(sum(i.min), sum(i.plus))
  if (any(n1 <= 1)) stop("Insufficient observations near discontinuity.")

  sigma2 <- c(stats::var(Y[i.min]), stats::var(Y[i.plus]))
  fc <- sum(n1) / (2 * n * h1)

  # Pilot third derivative and bandwidths for second derivatives
  m3 <- 6 * unname(stats::lm.fit(cbind(1, right, x, x^2, x^3), Y)$coefficients[5])

  h2 <- 7200^(1/7) *
    (sigma2 / (fc * m3^2))^(1/7) *
    c(sum(left), sum(right))^(-1/7)

  i.min <- left & x >= -h2[1]
  i.plus <- right & x <= h2[2]
  n2 <- c(sum(i.min), sum(i.plus))
  if (any(n2 <= 2)) stop("Insufficient observations near discontinuity.")

  m2 <- c(
    2 * unname(stats::lm.fit(cbind(1, x[i.min], x[i.min]^2), Y[i.min])$coefficients[3]),
    2 * unname(stats::lm.fit(cbind(1, x[i.plus], x[i.plus]^2), Y[i.plus])$coefficients[3])
  )

  # Regularization and optimal bandwidth
  r <- 2160 * sigma2 / (n2 * h2^4)
  CK <- switch(kernel,
               triangular   = 480^(1/5),
               uniform      = 144^(1/5),
               epanechnikov = (284160 / 847)^(1/5)
  )

  h.opt <- CK *
    (sum(sigma2) / (fc * ((m2[2] - m2[1])^2 + sum(r))))^(1/5) *
    n^(-1/5)

  # Kernel weights, normalized to sum to one
  u <- abs(x / h.opt)
  weights <- switch(kernel,
                    triangular   = pmax(1 - u, 0),
                    uniform      = as.numeric(u <= 1),
                    epanechnikov = pmax(1 - u^2, 0)
  )

  if (!any(left & u <= 1) || !any(right & u <= 1))
    stop("Insufficient observations in the calculated bandwidth.")

  list(
    bandwidth = unname(h.opt),
    weights = weights / sum(weights)
  )
}
