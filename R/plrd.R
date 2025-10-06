#' plrd: Partially Linear Regression Discontinuity Inference
#'
#' Optimized estimation and bias-aware inference for treatment effects identified by regression discontinuities
#' under partially linear model as proposed by Ghosh, Imbens and Wager (2025).
#'
#' @param Y The outcomes.
#' @param X The running variable.
#' @param threshold Sharp Regression Discontinuity threshold "c" at which CATE is to be estimated.
#' @param max.window A reasonable window for estimating the Lipschitz bound on curvature accurately. Default is the whole range.
#' @param W The treatment indicator. Currently plrd only supports sharp RDD, hence W is set as NULL by default.
#' @param Lipschitz.constant A Lipschitz bound for the second derivative of mu_w(x) = E(Y(w) | X = x).
#' @param sigma.sq An estimate of the homoskedastic variance of Y conditional on X
#' @param alpha Coverage probability of confidence intervals.
#' @param diff.curvatures Set true if user believes the curvatures are different before and after the threshold c.
#' @param signif.curvature Significance level for testing whether curvature is different before and after the threshold.
#' @param bin.width Bin width for discrete approximation.
#' @param num.bucket Number of bins for discrete approximation. Can only be used if bin.width = NULL.
#' @param use.spline Whether non-parametric components should be modeled as cubic splines
#'                   in order to reduce the number of optimization parameters, and potentially
#'                   improving computational performance.
#' @param spline.df Number of degrees of freedom (per running variable) used for spline computation.
#' @param seed Random seed for reproducibility.
#' @param use.homoskedastic.variance Whether confidence intervals should be built assuming homoskedasticity. As default we use FALSE (i.e., use heteroscedastic standard errors).
#' @param verbose Whether the optimizer should print progress information.
#' @return A trained plrd object containing the PLRD estimator, half width of confidence interval and other auxiliary objects such as MSE-optimal weights, worst-case bias, Lipschitz constant, etc.
#'
#' @references Ghosh, A., Imbens, G., & Wager, S. (2025).
#' PLRD: Partially Linear Regression Discontinuity Inference.
#' arXiv preprint arXiv:2503.09907.
#'
#' @examples
#' \donttest{
#' # Simple example of regression discontinuity design
#' set.seed(42)
#' n = 1000; threshold = 0
#' X = runif(n, -1, 1)
#' W = as.numeric(X >= threshold)
#' Y = (1 + 2*W)*(1 + X^2) + 1 / (1 + exp(X)) + rnorm(n, sd = .5)
#' out = plrd(Y, X, threshold)
#' print(out)
#' plot(out)
#' }
#'
#' @export
plrd <- function (Y, X, threshold, W = NULL,
                  max.window = NULL,
                  Lipschitz.constant = NULL,
                  sigma.sq = NULL,
                  alpha = 0.95,
                  diff.curvatures = FALSE,
                  signif.curvature = 0.001,
                  bin.width = NULL,
                  num.bucket = 400,
                  use.spline = TRUE,
                  spline.df = NULL,
                  seed = 42,
                  use.homoskedastic.variance = FALSE,
                  verbose = FALSE)
  {
  # Input checks:
  X = as.vector(X); Y = as.vector(Y)
  if(length(X)!=length(Y)){
    stop("Error: Running variable (X) and response (Y) cannot have different lengths!")
  }
  if((threshold < range(X)[1])||(threshold > range(X)[2])){
    stop("Error: Threshold must be within the range of the running variable.")
  }
  if(anyNA(Y) || anyNA(X)){
    stop("Error: Please remove any NA values present in X or Y.")
  }

  # Setting global max window (if max.window is not provided by the user):
  if(is.null(max.window)){
    max.window = max(abs(X - threshold))
  }

  # Computing W for sharp RDD (unless user is providing the treatment indicators):
  if(is.null(W)){
    W = as.integer(X>=threshold)
  }

  # Determining whether different curvatures should be used before and after the threshold
  Xc = X - threshold
  df = data.frame(Y = Y, Xc = Xc, W = W)
  if(is.null(diff.curvatures)){
      fit1 = stats::lm(Y ~ W * I(Xc) + I(Xc^2) + I(Xc^3),
                   data = df, subset = abs(Xc) <= max.window) # W interacts with only 1 and X - c
      fit2 = stats::lm(Y ~ W * (I(Xc) + I(Xc^2) + I(Xc^3)),
                   data = df, subset = abs(Xc) <= max.window) # W interacts with 1, X-c, (X-c)^2, (X-c)^3
      diff.curvatures = (as.numeric(stats::anova(fit1, fit2)[["Pr(>F)"]][2]) <= signif.curvature)
  }

  scale.Y = diff(range(df$Y)); Y.scaled = df$Y/scale.Y
  X.scaled = threshold + (df$Xc)/max.window
  if(is.null(dim(X.scaled))) X.scaled = matrix(X.scaled, ncol = 1)
  n = length(Y)

  # cross-fitting setup
  withr::with_seed(seed = seed,
                   if(is.null(Lipschitz.constant)){
                     # The following loop guards against unfortunate sample splits where B.fold1 and B.fold2 are vastly different
                     B.ratio = Inf; trial = 0; max.trial = 300  # Setting arbitrary B.fold1 and B.fold2 to initialize the following loop
                     while((B.ratio > 1.05)&(trial < max.trial)){
                       folds.new <- sample(rep(1:2, length.out = n))
                       fold2.idx <- which(folds.new == 1)
                       fold1.idx <- which(folds.new != 1)
                       # Estimating Lipschitz.constant (the Lipschitz constant for second derivative)
                       # B.fold1 and B.fold2 are to be used for fold1.idx and fold2.idx, respectively
                       B.fold2 = get.Lipschitz.constant(Y.scaled[fold1.idx], X.scaled[fold1.idx], threshold, diff.curvatures)
                       B.fold1 = get.Lipschitz.constant(Y.scaled[fold2.idx], X.scaled[fold2.idx], threshold, diff.curvatures)
                       if(B.ratio > max(B.fold1 / B.fold2, B.fold2 / B.fold1)){
                         folds = folds.new
                         B.ratio = max(B.fold1 / B.fold2, B.fold2 / B.fold1)
                       }
                       trial = trial + 1
                     }
                     fold2.idx <- which(folds == 1)
                     fold1.idx <- which(folds != 1)
                     B.fold2 = get.Lipschitz.constant(Y.scaled[fold1.idx], X.scaled[fold1.idx], threshold, diff.curvatures)
                     B.fold1 = get.Lipschitz.constant(Y.scaled[fold2.idx], X.scaled[fold2.idx], threshold, diff.curvatures)
                     if(trial==max.trial){
                       print(paste("The estimated Lipschitz constants for curvatures differ vastly across the sample splits, even after trying a bunch of different sample splits. We are proceeding with a more conservative estimate: taking the maximum across sample splits."))
                       B.conservative = max(B.fold1, B.fold2)
                       B.fold1 = B.fold2 = B.conservative
                     }
                   } else {
                     B.fold1 = B.fold2 = Lipschitz.constant
                     folds <- sample(rep(1:2, length.out = n))
                     fold2.idx <- which(folds == 1)
                     fold1.idx <- which(folds != 1)
                   }
  )
  if(verbose) print(paste("Estimated Lipschitz constant for curvature for cross-fitting folds:",
                          round(B.fold1*(max.window^3)*scale.Y,3),
                          round(B.fold2*(max.window^3)*scale.Y,3))) # Adjusted for scaling

  # Estimation of sigma.sq to be used for fold1.idx
  sigma.sq.fold1 = summary(stats::lm(Y.scaled[fold2.idx] ~ X.scaled[fold2.idx] * W[fold2.idx]))$sigma^2

  # Estimation of sigma.sq to be used for fold2.idx
  sigma.sq.fold2 = summary(stats::lm(Y.scaled[fold1.idx] ~ X.scaled[fold1.idx] * W[fold1.idx]))$sigma^2

  if(max(sigma.sq.fold1/sigma.sq.fold2, sigma.sq.fold2/sigma.sq.fold1) > 1.05){
    sigma.sq = max(sigma.sq.fold1, sigma.sq.fold2)
    sigma.sq.fold1 = sigma.sq.fold2 <- sigma.sq
  }

  if(verbose) print(paste("Estimated variance proxy for cross-fitting folds:",
                          round(sigma.sq.fold1*(scale.Y^2),3),
                          round(sigma.sq.fold2*(scale.Y^2),3)))

  if(is.null(spline.df)) spline.df = 40

  out_train = plrd.optim(Y = Y.scaled[fold1.idx],
                         X = X.scaled[fold1.idx],
                         threshold = threshold,
                         Lipschitz.constant = B.fold1,
                         sigma.sq = sigma.sq.fold1,
                         diff.curvatures = diff.curvatures,
                         alpha = alpha,
                         bin.width = bin.width,
                         num.bucket = num.bucket,
                         use.spline = use.spline,
                         spline.df = spline.df,
                         use.homoskedastic.variance = use.homoskedastic.variance,
                         verbose = verbose)

  out_test = plrd.optim(Y = Y.scaled[fold2.idx],
                        X = X.scaled[fold2.idx],
                        threshold = threshold,
                        Lipschitz.constant = B.fold2,
                        sigma.sq = sigma.sq.fold2,
                        diff.curvatures = diff.curvatures,
                        alpha = alpha,
                        bin.width = bin.width,
                        num.bucket = num.bucket,
                        use.spline = use.spline,
                        spline.df = spline.df,
                        use.homoskedastic.variance = use.homoskedastic.variance,
                        verbose = verbose)

  tau.hat = (out_train$tau.hat + out_test$tau.hat)/2 * scale.Y
  max.bias = (out_train$max.bias + out_test$max.bias)/2 * scale.Y
  Lipschitz.constant = (B.fold1+B.fold2)/2 * (max.window^3) * scale.Y
  gamma = rep(0, n)
  gamma[fold2.idx] = out_test$gamma/2; gamma[fold1.idx] = out_train$gamma/2
  sigma.sq = (sigma.sq.fold2 + sigma.sq.fold1)/2 * (scale.Y^2)
  if (use.homoskedastic.variance) {
    se.hat.tau = sqrt(sum(gamma^2 * sigma.sq))
  } else {
    regr.df = data.frame(X = X.scaled - threshold, W = W, Y = Y.scaled, gamma.sq = gamma^2)
    regr.df = regr.df[regr.df$gamma.sq != 0, ]
    Y.fit = stats::lm(Y ~ X * W, data = regr.df, weights = regr.df$gamma.sq)
    self.influence = stats::lm.influence(Y.fit)$hat
    Y.resid.sq = (regr.df$Y - stats::predict(Y.fit))^2
    se.hat.tau = sqrt(sum(Y.resid.sq * regr.df$gamma.sq/(1 - self.influence))) * scale.Y
  }

  tau.plusminus = get.plusminus(max.bias, se.hat.tau, alpha)
  gamma.fun.0 = rbind(out_train$gamma.fun.0, out_test$gamma.fun.0)
  gamma.fun.1 = rbind(out_train$gamma.fun.1, out_test$gamma.fun.1)
  gamma.fun.0[[1]] = (gamma.fun.0[[1]] - threshold)*max.window+threshold
  gamma.fun.1[[1]] = (gamma.fun.1[[1]] - threshold)*max.window+threshold

  rel.bias <- max.bias / se.hat.tau
  z.obs    <- abs(tau.hat) / se.hat.tau
  pval = stats::pnorm(rel.bias - z.obs) + stats::pnorm(-rel.bias - z.obs)

  ci.lower = tau.hat - tau.plusminus
  ci.upper = tau.hat + tau.plusminus

  # Creating the list to return
  ret = list(tau.hat = tau.hat,
             max.bias = max.bias,
             sampling.se = se.hat.tau,
             ci.lower = ci.lower,
             ci.upper = ci.upper,
             pval = pval,
             alpha = alpha,
             tau.plusminus = tau.plusminus,
             gamma = gamma,
             gamma.fun.0 = gamma.fun.0,
             gamma.fun.1 = gamma.fun.1,
             Lipschitz.constant = Lipschitz.constant,
             max.window = max.window,
             threshold = threshold,
             Y = Y,
             X = X,
             W = W,
             diff.curvatures = diff.curvatures
  )
  class(ret) = "plrd"
  return(ret)
}

#' plrd.optim: The optimizer we use to implement plrd with cross-fitting
#'
#' @param Y The outcomes.
#' @param X The running variable.
#' @param threshold Sharp Regression Discontinuity threshold "c" at which CATE is to be estimated.
#' @param W The treatment indicators, default is set as NULL for sharp RDD.
#' @param Lipschitz.constant A Lipschitz bound for the second derivative of mu_w(x) = E(Y(w) | X = x).
#' @param sigma.sq An estimate of the homoskedastic variance of Y conditional on X
#' @param diff.curvatures Set true if user believes the curvatures are different before and after the threshold c.
#' @param alpha Coverage probability of confidence intervals.
#' @param bin.width Bin width for discrete approximation.
#' @param num.bucket Number of bins for discrete approximation. Can only be used if bin.width = NULL.
#' @param use.spline Whether non-parametric components should be modeled as quadratic splines
#'                   in order to reduce the number of optimization parameters, and potentially
#'                   improving computational performance.
#' @param spline.df Number of degrees of freedom (per running variable) used for spline computation.
#' @param use.homoskedastic.variance Whether confidence intervals should be built assuming homoskedasticity.
#' @param verbose whether the optimizer should print progress information
#'
#' @return Parts of the desired plrd object for the specific sample split
#'
#' @keywords internal
#' @noRd
plrd.optim <- function (Y = NULL, X = NULL, threshold = NULL, W = NULL,
                        Lipschitz.constant = NULL,
                        sigma.sq = NULL,
                        diff.curvatures = FALSE,
                        alpha = NULL,
                        bin.width = NULL,
                        num.bucket = NULL,
                        use.spline = TRUE,
                        spline.df = NULL,
                        use.homoskedastic.variance = FALSE,
                        verbose = FALSE) {
  if(is.null(W)){
    W = as.numeric(X >= threshold)
  }

  # Homoscedastic estimation of sigma.sq
  if(is.null(dim(X))) X = matrix(X, ncol = 1)
  if(is.null(sigma.sq)){
    sigma.sq = summary(stats::lm(Y ~ X * W))$sigma^2
  }

  if(is.null(Lipschitz.constant)){
    Lipschitz.constant = get.Lipschitz.constant(Y, X, threshold)
  }

  # Binning X to num.bucket many bins
  bin.width = (max(X[, 1]) - min(X[, 1]))/num.bucket
  breaks = seq(min(X[, 1]) - bin.width/2, max(X[, 1]) + bin.width, by = bin.width)
  xx.grid = breaks[-1] - bin.width/2
  num.bucket = length(xx.grid) # +1 than before
  xx.grid = matrix(xx.grid, ncol = 1) # x.0 thru x.n, the bins are (x.i, x.{i+1}]
  xx.centered = matrix(xx.grid - threshold, ncol = 1)
  zeroed.idx = max(which(xx.centered < 0)) + c(0, 1, 2) # where the grid crosses the threshold

  idx.to.bucket = as.numeric(cut(X[, 1], breaks = breaks)) # putting X into bins (only 25 are non-empty)
  fact = factor(idx.to.bucket, levels = as.character(1:num.bucket)) # labeling as factors
  bucket.map = Matrix::sparse.model.matrix(~fact + 0, transpose = TRUE)

  X.counts = as.numeric(bucket.map %*% rep(1, length(W))) # number of X's in each bin
  W.counts = as.numeric(bucket.map %*% W) # number of X's in each bin past the threshold
  realized.idx.0 = which(X.counts > W.counts) # ids of bins b4 c that got some X in them
  realized.idx.1 = which(W.counts > 0) # ids of bins after c that got some X in them
  num.realized.0 = length(realized.idx.0)
  num.realized.1 = length(realized.idx.1)
  selector.0 = Matrix::Diagonal(num.bucket, 1)[realized.idx.0, ] # sparse NxN diag with 1 at realized.idx.0
  selector.1 = Matrix::Diagonal(num.bucket, 1)[realized.idx.1, ] # sparse NxN diag with 1 at realized.idx.1
  centering.matrix = Matrix::sparseMatrix(dims = c(length(zeroed.idx), num.bucket),
                                          i = 1:length(zeroed.idx), j = zeroed.idx,
                                          x = rep(1, length(zeroed.idx))) # sparse 3xN with Id.3 in the zeroed place
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
  D3 = Matrix::bandSparse(n = num.bucket - 3, m = num.bucket,
                          k = c(0, 1, 2, 3),
                          diag = list(rep(1, num.bucket),
                                      rep(-3, num.bucket),
                                      rep(3, num.bucket - 1),
                                      rep(-1, num.bucket - 2)))

  D3 = D3 %*% basis.mat # If f = basis.mat %*% coef, then D3 %*% coef gives the third derivatives
  selector.0 = selector.0 %*% basis.mat # Required for the G.i's before the threshold
  selector.1 = selector.1 %*% basis.mat # Required for the G.i's after the threshold
  centering.matrix = centering.matrix %*% basis.mat # Required for f(c) = f'(c) = f''(c) = 0
  num.df = ncol(basis.mat)

  # Setting up the quadratic programming problem
  num.lambda = 4 + 2 + as.numeric(diff.curvatures)
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
                     (xx.centered[realized.idx.0, ])^2,
                     if(diff.curvatures){
                       -(xx.centered[realized.idx.0, ])^2
                     },
                     selector.0
  ), # equality constraints for G.i's before the threshold
  cbind(Matrix::Matrix(0, num.realized.1, num.realized.0),
        Matrix::Diagonal(num.realized.1, -1),
        0, 1, 0,
        xx.centered[realized.idx.1, ],
        xx.centered[realized.idx.1, ],
        (xx.centered[realized.idx.1, ])^2,
        if(diff.curvatures){
          (xx.centered[realized.idx.1, ])^2
        },
        selector.1
  ), # equality constraints for G.i's after the threshold
  cbind(Matrix::Matrix(0, length(zeroed.idx), num.realized.0 + num.realized.1 + num.lambda),
        centering.matrix
  ), # 3 equality constraints to ensure f''(c) = f'(c) = f(c) = 0
  cbind(Matrix::Matrix(0, 2 * nrow(D3), num.realized.0 + num.realized.1),
        bin.width^3*Lipschitz.constant,
        matrix(0, 2 * nrow(D3), num.lambda - 1),
        rbind(D3, -D3)
  ) # constraints for third derivative: h^3 B (lambda.1/B) + D3 * f ≥ 0 and h^3 B (lambda.1/B) - D3 * f ≥ 0
  )

  meq = num.realized.0 + num.realized.1 + length(zeroed.idx)
  bvec = rep(0, nrow(Amat))
  gamma.0 = rep(0, num.bucket)
  gamma.1 = rep(0, num.bucket)

  if(verbose){
    print(paste0("Running quadprog with problem of size: ", dim(Amat)[1], " x ", dim(Amat)[2], "..."))
  }
  # Solve quadratic programming problem
  soln = quadprog::solve.QP(Matrix::Diagonal(num.params, Dmat.diagonal + 1e-6),
                            -dvec, Matrix::t(Amat), bvec, meq = meq)

  gamma.0[realized.idx.0] = -soln$solution[1:num.realized.0]/sigma.sq/2
  gamma.1[realized.idx.1] = -soln$solution[num.realized.0 + 1:num.realized.1]/sigma.sq/2
  t.hat = soln$solution[num.realized.0 + num.realized.1 + 1]/(2 * Lipschitz.constant) # the solution is lambda.1/B, hence solution/2B gives t.hat

  # Final gamma using on gamma.0 and gamma.1
  gamma = rep(0, length(W))
  gamma[W == 0] = as.numeric(Matrix::t(bucket.map[, which(W == 0)]) %*% gamma.0)
  gamma[W == 1] = as.numeric(Matrix::t(bucket.map[, which(W == 1)]) %*% gamma.1)
  gamma[W == 0] = -gamma[W == 0]/sum(gamma[W == 0])
  gamma[W == 1] = gamma[W == 1]/sum(gamma[W == 1])

  # Calculating max.bias and other statistics
  max.bias = Lipschitz.constant * t.hat
  tau.hat = sum(gamma * Y)

  # Creating the list to return
  ret = list(tau.hat = tau.hat,
             max.bias = max.bias,
             gamma = gamma,
             Lipschitz.constant = Lipschitz.constant,
             gamma.fun.0 = data.frame(xx=xx.grid[realized.idx.0,],
                                      gamma=gamma.0[realized.idx.0]),
             gamma.fun.1 = data.frame(xx=xx.grid[realized.idx.1,],
                                      gamma=gamma.1[realized.idx.1])
  )
  #  class(ret) = "plrd"
  return(ret)
}
