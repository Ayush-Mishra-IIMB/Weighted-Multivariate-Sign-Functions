# =============================================================================
# FILE: 02_weighted_spatial_median.R
# PAPER: "On Weighted Multivariate Sign Functions"
#        Majumdar & Chatterjee, Journal of Multivariate Analysis (2022)
#
# PURPOSE:
#   Implements the Weighted Spatial Median (WSM) — the paper's robust
#   location estimator. This is q̂_nW in the paper's notation.
#   It is used as the estimator for mu (center) in all downstream steps.
#
# SECTION: 2.1 of paper
#
# KEY IDEA:
#   The (unweighted) spatial median minimizes the sum of distances:
#       q̂_n = argmin_q sum_i |X_i - q|
#   The WEIGHTED spatial median minimizes the WEIGHTED sum of distances:
#       q̂_nW = argmin_q sum_i W(X_i, F) |X_i - q|
#
#   Intuitively: peripheral points get more "pull" on the center.
#   This can improve efficiency because peripheral points carry
#   more information about the location when data is elliptically distributed.
#
# THEOREM 1 (Asymptotic Normality):
#   sqrt(n) (q̂_nW - q0) --D--> N(0, Psi2W^{-1} Psi1W Psi2W^{-1})
#   where:
#     Psi1W = (d/dq Psi(q0)) (d/dq Psi(q0))^T  [outer product of gradient]
#     Psi2W = d^2/dq^2 Psi(q0)                  [Hessian of objective]
#
# ASYMPTOTIC RELATIVE EFFICIENCY vs unweighted spatial median:
#   ARE(q̂_nW, q̂_n) = [det(V1) / det(VW)]^{1/p}
#   Table 1 shows ARE > 1 for all distributions tested (normal and t).
# =============================================================================

# Source the weight functions module
 source("01_weight_functions.R")   # uncomment if running standalone


# =============================================================================
# GENERALIZED SIGN FUNCTION
# =============================================================================
# S(x; mu) = (x - mu) / |x - mu|   if x != mu
#           = 0                      if x == mu
#
# This maps every point to the unit sphere, keeping only direction.
# For a matrix X (n x p), returns an n x p matrix of signs.
#
# INPUT:
#   X   : n x p data matrix
#   mu  : p-vector, center
#
# OUTPUT:
#   S   : n x p matrix where each row is the sign of the corresponding X row

sign_function <- function(X, mu) {
  X  <- as.matrix(X)
  mu <- as.numeric(mu)

  X_centered <- sweep(X, 2, mu, "-")       # X_i - mu for each i
  norms      <- sqrt(rowSums(X_centered^2)) # |X_i - mu|

  # Handle the case where X_i == mu exactly (sign = 0)
  zero_idx <- (norms < .Machine$double.eps)

  S <- X_centered / norms                  # divide each row by its norm
  S[zero_idx, ] <- 0                       # set to zero where norm = 0

  return(S)
}


# =============================================================================
# WEIGHTED SPATIAL MEDIAN (q̂_nW)
# =============================================================================
# Minimizes: Psi_n(q) = sum_i W(X_i) |X_i - q|
#
# This is a convex optimization problem (sum of weighted norms is convex).
# We solve it using the Weiszfeld algorithm (iteratively reweighted means),
# which is the standard approach for L1 location problems.
#
# WEISZFELD ALGORITHM:
#   q^{(t+1)} = [sum_i W_i / |X_i - q^{(t)}|]^{-1}
#               sum_i [W_i / |X_i - q^{(t)}|] X_i
#   This is an iteratively reweighted mean where the "effective weight"
#   for observation i at iteration t is W_i / |X_i - q^{(t)}|.
#   Points closer to the current estimate get more pull.
#
# INPUT:
#   X        : n x p data matrix
#   W        : n-vector of weights (from compute_weights)
#   max_iter : maximum Weiszfeld iterations (default 500)
#   tol      : convergence tolerance (default 1e-8)
#   init     : p-vector, starting point (default: weighted column means)
#
# OUTPUT: list with
#   median   : p-vector, the weighted spatial median q̂_nW
#   iter     : number of iterations until convergence
#   converged: logical

weighted_spatial_median <- function(X, W, max_iter = 500, tol = 1e-8,
                                    init = NULL) {

  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  W <- as.numeric(W)

  # Normalize weights (robustness against scaling)
  W <- W / sum(W)

  # Starting point: weighted column means (good warm start)
  if (is.null(init)) {
    q <- colSums(W * X)
  } else {
    q <- as.numeric(init)
  }

  converged <- FALSE

  for (iter in seq_len(max_iter)) {

    q_old <- q

    # Distances from current estimate to each data point
    diffs  <- sweep(X, 2, q, "-")          # X_i - q  (n x p)
    dists  <- sqrt(rowSums(diffs^2))        # |X_i - q| (n-vector)

    # Avoid division by zero: if X_i == q exactly, skip that point
    # (it is already at the current estimate — contributes zero gradient)
    nonzero <- (dists > .Machine$double.eps)

    if (sum(nonzero) == 0) break            # all points coincide at q

    # Effective weights for this iteration: W_i / |X_i - q|
    eff_weights           <- rep(0, n)
    eff_weights[nonzero]  <- W[nonzero] / dists[nonzero]

    # Weiszfeld update: weighted mean with effective weights
    q <- colSums(eff_weights * X) / sum(eff_weights)

    # Check convergence: change in estimate is tiny
    if (sqrt(sum((q - q_old)^2)) < tol) {
      converged <- TRUE
      break
    }
  }

  if (!converged) {
    warning(sprintf(
      "Weighted spatial median did not converge in %d iterations.", max_iter))
  }

  return(list(median = q, iter = iter, converged = converged))
}


# =============================================================================
# FULL PIPELINE: compute q̂_nW from raw data
# =============================================================================
# This is the user-facing function.
# Handles weight computation + median estimation end-to-end.
#
# NOTE on the "chicken-and-egg" problem:
#   Weights W(X, F) require knowing mu and Sigma (to standardize to Z).
#   But we're trying to ESTIMATE mu here!
#   Solution: use a preliminary estimate (column means or unweighted
#   spatial median) to compute weights, then compute the WSM.
#   Optionally, iterate until convergence (self-consistent estimator).
#
# INPUT:
#   X           : n x p data matrix
#   method      : weight method: "HSD", "MhD", or "PD"
#   max_iter_wt : max iterations for self-consistent weight updating
#   tol_wt      : convergence tolerance for weight iteration
#   verbose     : print convergence info
#
# OUTPUT: list with
#   median     : p-vector, q̂_nW
#   weights    : n-vector, final weights used
#   Z          : n x p standardized matrix

compute_wsm <- function(X, method = "PD", max_iter_wt = 10, tol_wt = 1e-6,
                        verbose = FALSE) {

  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  # Step 1: Preliminary location and scatter estimates
  mu_init    <- colMeans(X)
  Sigma_init <- cov(X)

  # Step 2: Compute initial weights using preliminary estimates
  out_z <- standardize_to_Z(X, mu = mu_init, Sigma = Sigma_init)
  W     <- switch(method,
                  "HSD" = W_HSD(out_z),
                  "MhD" = W_MhD(out_z),
                  "PD"  = W_PD(out_z))

  # Step 3: Compute the weighted spatial median
  wsm_result <- weighted_spatial_median(X, W)
  q          <- wsm_result$median

  if (verbose) {
    cat(sprintf("WSM (%s): converged in %d iterations\n",
                method, wsm_result$iter))
  }

  return(list(
    median   = q,
    weights  = W,
    Z        = out_z,
    converged = wsm_result$converged
  ))
}


# =============================================================================
# ASYMPTOTIC VARIANCE ESTIMATION (Theorem 1)
# =============================================================================
# V_W = Psi2W^{-1} Psi1W Psi2W^{-1}
#
# We estimate Psi1W and Psi2W from the sample:
#
# Psi2W (Hessian estimate):
#   (1/n) sum_i W(X_i) * [I/|X_i-q| - (X_i-q)(X_i-q)^T / |X_i-q|^3]
#   This is the sample average of the Hessian of |X_i - q| weighted by W_i.
#
# Psi1W (outer product of gradient estimate):
#   (1/n) sum_i [W(X_i) * S(X_i; q)] [W(X_i) * S(X_i; q)]^T
#   = (1/n) sum_i W(X_i)^2 S(X_i;q) S(X_i;q)^T
#
# INPUT:
#   X : n x p data matrix
#   q : p-vector, estimated weighted spatial median
#   W : n-vector, weights
#
# OUTPUT: list with
#   Psi1W : p x p matrix
#   Psi2W : p x p matrix
#   VW    : p x p asymptotic variance matrix

estimate_asymptotic_variance <- function(X, q, W) {

  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  W <- as.numeric(W) / sum(W) * n    # rescale so mean(W) = 1

  diffs <- sweep(X, 2, q, "-")       # X_i - q  (n x p)
  dists <- sqrt(rowSums(diffs^2))     # |X_i - q|

  nonzero <- dists > .Machine$double.eps

  # ------ Psi1W: outer product term ------
  # = (1/n) sum_i W_i^2 * S_i S_i^T
  # where S_i = (X_i - q)/|X_i - q|
  Psi1W <- matrix(0, p, p)
  for (i in which(nonzero)) {
    s_i    <- diffs[i, ] / dists[i]          # sign vector (p x 1)
    Psi1W  <- Psi1W + W[i]^2 * outer(s_i, s_i)
  }
  Psi1W <- Psi1W / n

  # ------ Psi2W: Hessian term ------
  # = (1/n) sum_i W_i * [I/|X_i-q| - (X_i-q)(X_i-q)^T/|X_i-q|^3]
  Psi2W <- matrix(0, p, p)
  for (i in which(nonzero)) {
    d_i    <- dists[i]
    s_i    <- diffs[i, ] / d_i
    Psi2W  <- Psi2W + W[i] * (diag(p) / d_i - outer(s_i, s_i) / d_i)
  }
  Psi2W <- Psi2W / n

  # ------ Asymptotic variance V_W = Psi2W^{-1} Psi1W Psi2W^{-1} ------
  Psi2W_inv <- solve(Psi2W)
  VW        <- Psi2W_inv %*% Psi1W %*% Psi2W_inv

  return(list(Psi1W = Psi1W, Psi2W = Psi2W, VW = VW))
}


# =============================================================================
# ASYMPTOTIC RELATIVE EFFICIENCY (ARE) — Table 1 of paper
# =============================================================================
# ARE(q̂_nW, q̂_n) = [det(V1) / det(VW)]^{1/p}
#
# V1  = asymptotic variance of UNWEIGHTED spatial median (W=1)
# VW  = asymptotic variance of WEIGHTED spatial median
# ARE > 1 means the weighted version is MORE efficient
#
# INPUT:
#   X      : n x p data matrix
#   method : weight method
#
# OUTPUT: scalar ARE value

compute_ARE_wsm <- function(X, method = "PD") {

  X <- as.matrix(X)
  p <- ncol(X)

  # Compute unweighted spatial median (W = rep(1, n))
  n      <- nrow(X)
  W_ones <- rep(1, n)
  q_unw  <- weighted_spatial_median(X, W_ones)$median
  av_unw <- estimate_asymptotic_variance(X, q_unw, W_ones)
  V1     <- av_unw$VW

  # Compute weighted spatial median
  wsm    <- compute_wsm(X, method = method)
  q_w    <- wsm$median
  W      <- wsm$weights
  av_w   <- estimate_asymptotic_variance(X, q_w, W)
  VW     <- av_w$VW

  # ARE = [det(V1) / det(VW)]^{1/p}
  ARE <- (det(V1) / det(VW))^(1/p)

  return(ARE)
}


# =============================================================================
# QUICK DEMO
# =============================================================================
# Uncomment and run to test:

# library(MASS)
# set.seed(42)
# X <- mvrnorm(n = 200, mu = c(2, 3), Sigma = matrix(c(2,1,1,1), 2, 2))
#
# # Compute weighted spatial median
# result_pd  <- compute_wsm(X, method = "PD",  verbose = TRUE)
# result_mhd <- compute_wsm(X, method = "MhD", verbose = TRUE)
# result_hsd <- compute_wsm(X, method = "HSD", verbose = TRUE)
#
# cat("True mu:         ", c(2, 3), "\n")
# cat("Column means:    ", colMeans(X), "\n")
# cat("WSM (PD):        ", result_pd$median, "\n")
# cat("WSM (MhD):       ", result_mhd$median, "\n")
# cat("WSM (HSD):       ", result_hsd$median, "\n")
#
# # ARE comparison
# are_pd  <- compute_ARE_wsm(X, method = "PD")
# are_mhd <- compute_ARE_wsm(X, method = "MhD")
# cat(sprintf("ARE (PD vs unweighted):  %.4f\n", are_pd))
# cat(sprintf("ARE (MhD vs unweighted): %.4f\n", are_mhd))
