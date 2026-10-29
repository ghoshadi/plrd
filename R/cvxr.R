#' Solve the plrd minimax problem in its primal form with CVXR
#'
#' Solves problem (31) of Ghosh, Imbens and Wager (2025) directly on the grid set
#' up by plrd.optim. The supremum over the function class has a closed form: on the
#' grid, f is determined by its third differences v = D3 f (with |v| <= h^3) and
#' the three centering conditions C f = 0. Since M = rbind(C, D3) is square and
#' invertible, sup_f <g, f> = h^3 * ||(M^{-T} g)[-(1:3)]||_1, an L1 norm that is
#' affine in the weights. No dual, spline basis or ridge is needed.
#'
#' @param xx.centered Grid points minus the threshold (num.bucket x 1).
#' @param bin.width Grid spacing h.
#' @param zeroed.idx The three grid indices where the grid crosses the threshold.
#' @param counts.0,counts.1 Number of control / treated observations in each bin.
#' @param realized.idx.0,realized.idx.1 Bins holding at least one control / treated observation.
#' @param Lipschitz.constant Lipschitz bound B for the second derivative.
#' @param sigma.sq Homoskedastic variance.
#' @param diff.curvatures Whether the curvatures differ on the two sides of the threshold.
#' @param verbose Whether the solver should print progress information.
#' @return A list with per-bin weights gamma.0 and gamma.1 (length num.bucket) and
#'   t.hat, the worst-case bias divided by Lipschitz.constant.
#' @keywords internal
#' @noRd
plrd.primal.cvxr <- function(xx.centered, bin.width, zeroed.idx,
                             counts.0, counts.1,
                             realized.idx.0, realized.idx.1,
                             Lipschitz.constant, sigma.sq,
                             diff.curvatures = FALSE, verbose = FALSE) {
  num.bucket = length(xx.centered)

  # Bias operator: g |-> h^3 * (M^{-T} g)[-(1:3)], M = rbind(centering, D3)
  D3 = Matrix::bandSparse(n = num.bucket - 3, m = num.bucket,
                          k = c(0, 1, 2, 3),
                          diag = list(rep(1, num.bucket),
                                      rep(-3, num.bucket),
                                      rep(3, num.bucket - 1),
                                      rep(-1, num.bucket - 2)))
  centering.matrix = Matrix::sparseMatrix(dims = c(length(zeroed.idx), num.bucket),
                                          i = 1:length(zeroed.idx), j = zeroed.idx,
                                          x = rep(1, length(zeroed.idx)))
  M = as.matrix(rbind(centering.matrix, D3))
  bias.op = bin.width^3 * solve(t(M))[-(1:3), , drop = FALSE]

  # One weight per realized bin; the bin mass multiplies it everywhere
  w0 = counts.0[realized.idx.0]; w1 = counts.1[realized.idx.1]
  x0 = xx.centered[realized.idx.0]; x1 = xx.centered[realized.idx.1]
  A0 = bias.op[, realized.idx.0, drop = FALSE] %*% diag(w0, length(w0))
  A1 = bias.op[, realized.idx.1, drop = FALSE] %*% diag(w1, length(w1))

  # Scaling, for conditioning only (the problem is unchanged): the weights enter as
  # u = n * gamma, the bias bound as s = B * n * t, and the operator is normalized to
  # unit max entry. Without it CLARABEL can stop short ("optimal_inaccurate") when the
  # threshold is far from the center of the scaled running variable.
  n.obs = sum(w0) + sum(w1)
  op.scale = 1 / max(abs(A0), abs(A1))
  A0 = op.scale * A0; A1 = op.scale * A1
  u0 = CVXR::Variable(length(w0)); u1 = CVXR::Variable(length(w1)); s = CVXR::Variable(1)
  bias = if (diff.curvatures) {
    CVXR::p_norm(A0 %*% u0, 1) + CVXR::p_norm(A1 %*% u1, 1)
  } else {
    CVXR::p_norm(A0 %*% u0 + A1 %*% u1, 1)
  }
  constraints = list(
    sum(w1 * u1) == n.obs,
    sum(w0 * u0) == -n.obs,
    sum(w1 * x1 * u1) == 0,
    sum(w0 * x0 * u0) == 0,
    (Lipschitz.constant / op.scale) * bias <= s
  )
  constraints = c(constraints,
                  if (diff.curvatures) {
                    list(sum(w1 * x1^2 * u1) == 0, sum(w0 * x0^2 * u0) == 0)
                  } else {
                    list(sum(w1 * x1^2 * u1) + sum(w0 * x0^2 * u0) == 0)
                  })
  # n^2 times sigma.sq * sum(gamma^2) + B^2 * t^2
  objective = sigma.sq * (CVXR::sum_squares(sqrt(w0) * u0) + CVXR::sum_squares(sqrt(w1) * u1)) +
    CVXR::square(s)

  prob = CVXR::Problem(CVXR::Minimize(objective), constraints)
  # CLARABEL is named explicitly: it ships with CVXR, and automatic selection could
  # otherwise pick a different (possibly commercial) solver depending on the machine.
  # tol_gap_rel = 1e-7 (CLARABEL's default is 1e-8) is set only to avoid a spurious
  # "Solution may be inaccurate" warning. With the threshold in the tail of X (e.g.
  # settings 6 and 7 of tests/testthat/test_plot.R), CLARABEL can stall at a relative
  # gap of 2-3e-8 with feasibility already met (residuals <= 4e-9) and report
  # "optimal_inaccurate"; the solutions in those cases match PIQP's to 5e-9 relative
  # in objective. 1e-7 is roughly single precision.
  CVXR::psolve(prob, solver = "CLARABEL", tol_gap_rel = 1e-7, verbose = verbose)
  status = CVXR::status(prob)
  if (!(status %in% c("optimal", "optimal_inaccurate"))) {
    stop("Error: the optimization problem was not solved (solver status: ", status, ").")
  }

  gamma.0 = rep(0, num.bucket); gamma.1 = rep(0, num.bucket)
  gamma.0[realized.idx.0] = as.numeric(CVXR::value(u0)) / n.obs
  gamma.1[realized.idx.1] = as.numeric(CVXR::value(u1)) / n.obs
  t.hat = as.numeric(CVXR::value(s)) / (Lipschitz.constant * n.obs)
  list(gamma.0 = gamma.0, gamma.1 = gamma.1, t.hat = t.hat)
}
