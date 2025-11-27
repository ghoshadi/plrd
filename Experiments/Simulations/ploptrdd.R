#' plrd: Partially Linear Regression Discontinuity
#'
#' Description: TBA
#'
#' @param Y The outcomes.
#' @param X The running variable.
#' @param cutoff Threshold "c" at which CATE is to be estimated. Default is c = 0.
#' @param max.window A reasonable window for estimating the Lipschitz bound accurately. Default is the whole range.
#' @param max.second.derivative A Lipschitz bound for the second derivative of $\mu_w(x) = E(Y(w) | X = x)$.
#' @param model.for.estim.deriv Model for estimating the Lipschitz bound. Default is a polynomial regression.
#' @param diff.curvatures Set true if user believes the curvatures are different before and after the threshold c.
#' @param alpha Coverage probability of confidence intervals.
#' @param bin.width Bin width for discrete approximation.
#' @param num.bucket Number of bins for discrete approximation. Can only be used if bin.width = NULL.
#' @param use.spline Whether non-parametric components should be modeled as quadratic splines
#'                   in order to reduce the number of optimization parameters, and potentially
#'                   improving computational performance.
#' @param spline.df Number of degrees of freedom (per running variable) used for spline computation.
#' @param use.homoskedastic.variance Whether confidence intervals should be built assuming homoskedasticity.
#' @param try.elnet.for.sigma.sq Whether an elastic net on a spline basis should be used for estimating sigma^2.
#' @param verbose whether the optimizer should print progress information
#'
#' @return A trained plrd object. (similar to optrdd (Imbens & Wager, 2017) object)
#'
#' @references Imbens, G., & Wager, S. (2017).
#' Optimized Regression Discontinuity Designs.
#' arXiv preprint arXiv:1705.01677.
#'
#' @examples
#' # Simple regression discontinuity with discrete X
#' n = 4000; threshold = 0
#' X = sample(seq(-1, 1, by = .05), n, replace = TRUE)
#' W = as.numeric(X >= threshold)
#' Y = 0.4 * W + 1 / (1 + exp(2 * X)) + 0.2 * rnorm(n)
#' print(out <- plrd(Y, X)); plot(out)
#'
#' @export
ploptrdd <- function (Y = NULL, X, cutoff = 0,
                  max.window = NULL,
                  max.second.derivative = NULL,
                  sigma.sq = NULL,
                  model.for.estim.deriv = NULL,
                  diff.curvatures = FALSE,
                  alpha = 0.95,
                  bin.width = NULL,
                  num.bucket = 400,
                  use.spline = TRUE,
                  spline.df = 40,
                  use.homoskedastic.variance = TRUE,
                  try.elnet.for.sigma.sq = FALSE,
                  verbose = TRUE) {
  if(is.null(max.window)){
    if(is.null(dim(X))) n = length(X)
    if(!is.null(dim(X))) n = dim(X)[1]
    max.window = max(abs(X-cutoff))
  }
  datsub.indices <- which(abs(X - cutoff) <= max.window)
  X = X[datsub.indices]; Y = Y[datsub.indices]; n = length(Y)
  # cross-fitting setup
  if(is.null(dim(X))) X = matrix(X, ncol = 1)
  W = as.numeric(X>cutoff)
  folds <- sample(rep(1:2, length.out = n))
  test_idx <- which(folds == 1)
  train_idx <- which(folds != 1)
  # Homoscedastic estimation of sigma.sq to be used for train_idx
  sigma.sq_train = summary(lm(Y[test_idx] ~ X[test_idx] * W[test_idx]))$sigma^2
  # Elastic net for sigma.sq (in case user wants it)
  if (try.elnet.for.sigma.sq) {
    if (ncol(X) > 1) {
      stop("Elastic net for sigma squared not implemented with more than 1 running variable.")
    }
    linear.params = 1 + 2 * ncol(X)
    elnet.df = 7
    ridge.mat = cbind(W[test_idx], X[test_idx,], W[test_idx] * X[test_idx,], matrix(0, length(W[test_idx]), 2 * elnet.df))
    ridge.mat[W[test_idx] == 0, linear.params + 1:elnet.df] = splines::ns(X[test_idx][W[test_idx] == 0], df = elnet.df)
    ridge.mat[W[test_idx] == 1, linear.params + elnet.df + 1:elnet.df] = splines::ns(X[test_idx][W[test_idx] == 1], df = elnet.df)
    elnet = glmnet::cv.glmnet(ridge.mat, Y[test_idx],
                              penalty.factor = c(rep(0, linear.params), rep(1, 2 * elnet.df)),
                              keep = TRUE, alpha = 0.5)
    elnet.hat = elnet$fit.preval[, !is.na(colSums(elnet$fit.preval)),
                                 drop = FALSE][, elnet$lambda == elnet$lambda.1se]
    sigma.sq.elnet = mean((elnet.hat - Y[test_idx])^2)
    sigma.sq_train = min(sigma.sq_train, sigma.sq.elnet)
  }
  # Homoscedastic estimation of sigma.sq to be used for test_idx
  sigma.sq_test = summary(lm(Y[train_idx] ~ X[train_idx] * W[train_idx]))$sigma^2
  # Elastic net for sigma.sq (in case user wants it)
  if (try.elnet.for.sigma.sq) {
    if (ncol(X) > 1) {
      stop("Elastic net for sigma squared not implemented with more than 1 running variable.")
    }
    linear.params = 1 + 2 * ncol(X)
    elnet.df = 7
    ridge.mat = cbind(W[train_idx], X[train_idx,], W[train_idx] * X[train_idx,], matrix(0, length(W[train_idx]), 2 * elnet.df))
    ridge.mat[W[train_idx] == 0, linear.params + 1:elnet.df] = splines::ns(X[train_idx][W[train_idx] == 0], df = elnet.df)
    ridge.mat[W[train_idx] == 1, linear.params + elnet.df + 1:elnet.df] = splines::ns(X[train_idx][W[train_idx] == 1], df = elnet.df)
    elnet = glmnet::cv.glmnet(ridge.mat, Y[train_idx],
                              penalty.factor = c(rep(0, linear.params), rep(1, 2 * elnet.df)),
                              keep = TRUE, alpha = 0.5)
    elnet.hat = elnet$fit.preval[, !is.na(colSums(elnet$fit.preval)),
                                 drop = FALSE][, elnet$lambda == elnet$lambda.1se]
    sigma.sq.elnet = mean((elnet.hat - Y[train_idx])^2)
    sigma.sq_test = min(sigma.sq_test, sigma.sq.elnet)
  }

  # Estimating max.second.derivative (the Lipschitz constant for second derivative)
  # B_train and B_test are to be used for train_idx and test_idx, respectively
  if(is.null(max.second.derivative)){
    B_test = get.max.deriv(Y[train_idx], X[train_idx], cutoff, d=2)
    B_train = get.max.deriv(Y[test_idx], X[test_idx], cutoff, d=2)
    if(verbose) print(c(B_train,B_test))
  }
  if(!is.null(max.second.derivative)) B_train = B_test = max.second.derivative

  out_train = ploptrdd.optim(Y[train_idx], X[train_idx], cutoff, max.window,
                   max.second.derivative = B_train, sigma.sq = sigma.sq_train,
                   model.for.estim.deriv, diff.curvatures,
                   alpha, bin.width, num.bucket, use.spline, spline.df,
                   use.homoskedastic.variance, try.elnet.for.sigma.sq, verbose)

  out_test = ploptrdd.optim(Y[test_idx], X[test_idx], cutoff, max.window,
                   max.second.derivative = B_test, sigma.sq = sigma.sq_test,
                   model.for.estim.deriv, diff.curvatures,
                   alpha, bin.width, num.bucket, use.spline, spline.df,
                   use.homoskedastic.variance, try.elnet.for.sigma.sq, verbose)

  tau.hat = (out_train$tau.hat + out_test$tau.hat)/2
  max.bias = (out_train$max.bias + out_test$max.bias)/2
  gamma = rep(0, n)
  gamma[test_idx] = out_test$gamma/2; gamma[train_idx] = out_train$gamma/2
  sigma.sq = (sigma.sq_test + sigma.sq_train)/2
  if (use.homoskedastic.variance) {
    #se.hat.tau = sqrt(sum(gamma^2 * sigma.sq))
    se.hat.tau = sqrt(sum(gamma[train_idx]^2 * sigma.sq_train) + sum(gamma[test_idx]^2 * sigma.sq_test))
  }
  else {
    regr.df = data.frame(X = X, W = W, Y = Y, gamma.sq = gamma^2)
    regr.df = regr.df[regr.df$gamma.sq != 0, ]
    Y.fit = stats::lm(Y ~ X * W, data = regr.df, weights = regr.df$gamma.sq)
    self.influence = stats::lm.influence(Y.fit)$hat
    Y.resid.sq = (regr.df$Y - stats::predict(Y.fit))^2
    se.hat.tau = sqrt(sum(Y.resid.sq * regr.df$gamma.sq/(1 - self.influence)))
    if (!try.elnet.for.sigma.sq & se.hat.tau^2 < 0.8 * sum(regr.df$gamma.sq) * sigma.sq) {
      warning(paste("Implicit estimate of sigma^2 may be too pessimistic,",
                    "resulting in test but needlessly long confidence intervals.",
                    "Try the option `try.elnet.for.sigma.sq = TRUE` for potentially improved performance."))
    }
  }

  tau.plusminus = get.plusminus(max.bias, se.hat.tau, alpha)

  rho.hat = rep(0, n)
  rho.hat[test_idx] = B_train*out_test$rho; rho.hat[train_idx] = B_test*out_train$rho
  Lipschitz.constant = (B_train + B_test)/2

  # Creating the list to return
  ret = list(tau.hat = tau.hat,
             tau.plusminus = tau.plusminus,
             alpha = alpha,
             max.bias = max.bias,
             sampling.se = se.hat.tau,
             gamma = gamma,
             Lipschitz.constant = max.second.derivative,
             rho = rho.hat
             )
  class(ret) = "plrd"
  return(ret)
}

ploptrdd.optim <- function (Y = NULL, X, cutoff = 0,
                  max.window = NULL,
                  max.second.derivative = NULL,
                  sigma.sq = NULL,
                  model.for.estim.deriv = NULL,
                  diff.curvatures = FALSE,
                  alpha = 0.95,
                  bin.width = NULL,
                  num.bucket = 400,
                  use.spline = TRUE,
                  spline.df = 40,
                  use.homoskedastic.variance = TRUE,
                  try.elnet.for.sigma.sq = FALSE,
                  verbose = TRUE) {
  if(is.null(max.window)){
    max.window = max(abs(X - cutoff))
  }
  datsub.indices <- which(abs(X - cutoff) <= max.window)
  X = X[datsub.indices]; Y = Y[datsub.indices]
  W = as.numeric(X > cutoff)
  # Homoscedastic estimation of sigma.sq
  if(is.null(dim(X))) X = matrix(X, ncol = 1)
  if(is.null(sigma.sq)){
    sigma.sq = summary(lm(Y ~ X * W))$sigma^2
    # Elastic net for sigma.sq (in case user wants it)
    if (try.elnet.for.sigma.sq) {
      if (ncol(X) > 1) {
        stop("Elastic net for sigma squared not implemented with more than 1 running variable.")
      }
      linear.params = 1 + 2 * ncol(X)
      elnet.df = 7
      ridge.mat = cbind(W, X, W * X, matrix(0, length(W), 2 * elnet.df))
      ridge.mat[W == 0, linear.params + 1:elnet.df] = splines::ns(X[W == 0, ], df = elnet.df)
      ridge.mat[W == 1, linear.params + elnet.df + 1:elnet.df] = splines::ns(X[W == 1, ], df = elnet.df)
      elnet = glmnet::cv.glmnet(ridge.mat, Y,
                                penalty.factor = c(rep(0, linear.params), rep(1, 2 * elnet.df)),
                                keep = TRUE, alpha = 0.5)
      elnet.hat = elnet$fit.preval[, !is.na(colSums(elnet$fit.preval)),
                                   drop = FALSE][, elnet$lambda == elnet$lambda.1se]
      sigma.sq.elnet = mean((elnet.hat - Y)^2)
      sigma.sq = min(sigma.sq, sigma.sq.elnet)
    }
  }

  # Estimating max.second.derivative (the Lipschitz constant for second derivative)
  # if(is.null(max.second.derivative)){
  #   pre.idx = (X < cutoff)
  #   max.second.derivative = max(get.max.derivative(Y[pre.idx], (X-cutoff)[pre.idx], d = 2, model = model.for.estim.deriv),
  #                              get.max.derivative(Y[-pre.idx], (X-cutoff)[-pre.idx], d = 2, model = model.for.estim.deriv))
  # }
  if(is.null(max.second.derivative)){
    pre.idx = (X < cutoff)
    max.second.derivative = 2*max(get.max.deriv(Y[pre.idx], (X-cutoff)[pre.idx], cutoff, d = 2),
                               get.max.deriv(Y[!pre.idx], (X-cutoff)[!pre.idx], cutoff, d = 2))
  }


  # Binning X to num.bucket many bins
  bin.width = (max(X[, 1]) - min(X[, 1]))/num.bucket
  breaks = seq(min(X[, 1]) - bin.width/2, max(X[, 1]) + bin.width, by = bin.width)
  xx.grid = breaks[-1] - bin.width/2
  num.bucket = length(xx.grid) # +1 than before
  xx.grid = matrix(xx.grid, ncol = 1) # x.0 thru x.n, the bins are (x.i, x.{i+1}]
  xx.centered = matrix(xx.grid - cutoff, ncol = 1)
  zeroed.idx = max(which(xx.centered < 0)) + c(0, 1) # where the grid crosses the cutoff

  idx.to.bucket = as.numeric(cut(X[, 1], breaks = breaks)) # putting X into bins (only 25 are non-empty)
  fact = factor(idx.to.bucket, levels = as.character(1:num.bucket)) # labeling as factors
  bucket.map = Matrix::sparse.model.matrix(~fact + 0, transpose = TRUE)

  X.counts = as.numeric(bucket.map %*% rep(1, length(W))) # number of X's in each bin
  W.counts = as.numeric(bucket.map %*% W) # number of X's in each bin past the cutoff
  realized.idx.0 = which(X.counts > W.counts) # ids of bins b4 c that got some X in them
  realized.idx.1 = which(W.counts > 0) # ids of bins after c that got some X in them
  num.realized.0 = length(realized.idx.0)
  num.realized.1 = length(realized.idx.1)
  selector.0 = Matrix::Diagonal(num.bucket, 1)[realized.idx.0, ] # sparse NxN diag with 1 at realized.idx.0
  selector.1 = Matrix::Diagonal(num.bucket, 1)[realized.idx.1, ] # sparse NxN diag with 1 at realized.idx.1
  centering.matrix = Matrix::sparseMatrix(dims = c(length(zeroed.idx), num.bucket),
                                          i = 1:length(zeroed.idx), j = zeroed.idx,
                                          x = rep(1, length(zeroed.idx))) # sparse 2xN with Id.2 in the zeroed place
  if(use.spline){
    # Setting up B-spline basis; spline.df = 40 as default
    basis.mat.raw = splines::bs(xx.grid, degree = 3, df = spline.df, intercept = TRUE)
    class(basis.mat.raw) = "matrix"
    basis.mat = Matrix::Matrix(basis.mat.raw)
  }
  else{
    basis.mat = diag(num.bucket)
  }

  # Difference matrix for third derivative
  # D2 = Matrix::bandSparse(n = num.bucket - 3, m = num.bucket,
  #                         k = c(0, 1, 2, 3),
  #                         diag = list(rep(1, num.bucket),
  #                                     rep(-3, num.bucket),
  #                                     rep(3, num.bucket - 1),
  #                                     rep(-1, num.bucket - 2)))
  D2 = Matrix::bandSparse(n = num.bucket - 2, m = num.bucket, k = c(0, 1, 2),
                          diag = list(rep(1, num.bucket), rep(-2, num.bucket))[c(1, 2, 1)])


  D2 = D2 %*% basis.mat # If f = basis.mat %*% coef, then D2 %*% coef gives the second derivatives
  selector.0 = selector.0 %*% basis.mat # Required for the G.i's before the threshold
  selector.1 = selector.1 %*% basis.mat # Required for the G.i's after the threshold
  centering.matrix = centering.matrix %*% basis.mat # Required for f'(c) = f''(c) = 0
  num.df = ncol(basis.mat)

  # Setting up the quadratic programming problem
  num.lambda = 5 + as.numeric(diff.curvatures)
  num.params = num.realized.0 + num.realized.1 + num.lambda + num.df

  Dmat.diagonal = c((X.counts[realized.idx.0] - W.counts[realized.idx.0])/2/sigma.sq,
                    (W.counts[realized.idx.1])/2/sigma.sq,
                    1/2,
                    rep(0, num.lambda - 1 + num.df))

  dvec = c(rep(0, num.realized.0 + num.realized.1 + 1), 1, -1, rep(0, num.lambda - 3 + num.df))

  # Constraint matrix
  Amat = rbind(cbind(Matrix::Diagonal(num.realized.0, -1),
                     Matrix::Matrix(0, num.realized.0, num.realized.1),
                     0, 0, 1,
                     xx.centered[realized.idx.0, ],
                     -xx.centered[realized.idx.0, ],
                     selector.0), # equality constraints for G.i's before the threshold
               cbind(Matrix::Matrix(0, num.realized.1, num.realized.0),
                     Matrix::Diagonal(num.realized.1, -1),
                     0, 1, 0,
                     xx.centered[realized.idx.1, ],
                     xx.centered[realized.idx.1, ],
                     selector.1), # equality constraints for G.i's after the threshold
               cbind(Matrix::Matrix(0, length(zeroed.idx),
                                    num.realized.0 + num.realized.1 + num.lambda),
                     centering.matrix), # 2 equality constraints to ensure f'(c) = f(c) = 0
               cbind(Matrix::Matrix(0, 2 * nrow(D2), num.realized.0 + num.realized.1),
                     bin.width^2*max.second.derivative,
                     matrix(0, 2 * nrow(D2), num.lambda - 1),
                     rbind(D2, -D2)
                     ) # constraints for third derivative: h^2 B (lambda.1/B) + D2 * f ≥ 0 and h^2 B (lambda.1/B) - D2 * f ≥ 0
  )

  meq = num.realized.0 + num.realized.1 + length(zeroed.idx)
  bvec = rep(0, nrow(Amat))
  gamma.0 = rep(0, num.bucket)
  gamma.1 = rep(0, num.bucket)

  if(verbose){
    print(paste0("Running quadrprog with problem of size: ", dim(Amat)[1], " x ", dim(Amat)[2], "..."))
  }
  # Solve quadratic programming problem
  soln = quadprog::solve.QP(Matrix::Diagonal(num.params, Dmat.diagonal + 1e-6),
                            -dvec, Matrix::t(Amat), bvec, meq = meq)

  gamma.0[realized.idx.0] = -soln$solution[1:num.realized.0]/sigma.sq/2
  gamma.1[realized.idx.1] = -soln$solution[num.realized.0 + 1:num.realized.1]/sigma.sq/2
  t.hat = soln$solution[num.realized.0 + num.realized.1 + 1]/(2 * max.second.derivative) # the solution is lambda.1/B, hence solution/2B gives t.hat

  # Final gamma using on gamma.0 and gamma.1
  gamma = rep(0, length(W))
  gamma[W == 0] = as.numeric(Matrix::t(bucket.map[, which(W == 0)]) %*% gamma.0)
  gamma[W == 1] = as.numeric(Matrix::t(bucket.map[, which(W == 1)]) %*% gamma.1)
  gamma[W == 0] = -gamma[W == 0]/sum(gamma[W == 0])
  gamma[W == 1] = gamma[W == 1]/sum(gamma[W == 1])

  # Calculating max.bias and other statistics
  max.bias = max.second.derivative * t.hat
  tau.hat = sum(gamma * Y)

  lambda.minus.one = soln$solution[num.realized.0 + num.realized.1 + 2:(6+as.numeric(diff.curvatures))]
  if(diff.curvatures){
    rho.hat = ((-2)*sigma.sq*gamma -cbind(W, 1-W, X, (2*W-1)*X, W*X^2, (1-W)*X^2) %*% lambda.minus.one)/(soln$solution[num.realized.0 + num.realized.1+1]*max.second.derivative)

  }
  else{
    rho.hat = ((-2)*sigma.sq*gamma -cbind(W, 1-W, X, (2*W-1)*X, X^2) %*% lambda.minus.one)/(soln$solution[num.realized.0 + num.realized.1+1]*max.second.derivative)
  }
  # need to change above line if diff.curvatures = TRUE
  #tau.plusminus = get.plusminus(max.bias, se.hat.tau, alpha)

  # Creating the list to return
  ret = list(tau.hat = tau.hat,
             max.bias = max.bias,
             gamma = gamma,
             Lipschitz.constant = max.second.derivative,
             rho = rho.hat
  )
#  class(ret) = "plrd"
  return(ret)
}
#'
#' Bounds on higher order derivatives of the conditional response function
#'
#' @param y The outcomes.
#' @param x The running variable.
#' @param cutoff cutoff for treatment.
#' @param d Order of the derivative (default is d = 3).
#' @return Bound on d-th derivative of the conditional response function.
#' @export
#'
get.max.deriv <- function(y, x, cutoff, d = 2){
  fit_left <- stats::lm(y[x < cutoff] ~ poly((x-cutoff)[x < cutoff], d, raw = TRUE))
  fit_right <- stats::lm(y[x > cutoff] ~ poly((x-cutoff)[x > cutoff], d, raw = TRUE))
  eps = stats::sd(y)/100
  B_left <- abs(factorial(d) * stats::coef(fit_left)[(d+1)])
  B_right <- abs(factorial(d) * stats::coef(fit_right)[(d+1)])
  return(B_hat <- max(c(B_left, B_right, eps)))
}
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
