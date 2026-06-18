# =============================================================================
# FILE: 04_affine_equivariant.R
# PAPER: "On Weighted Multivariate Sign Functions"
#        Majumdar & Chatterjee, Journal of Multivariate Analysis (2022)
#
# PURPOSE:
#   Implements Σ* (Sigma-star) — the affine equivariant version of Σ̃.
#   Σ̃ is NOT affine equivariant; Σ* fixes this.
#
# SECTION: 3 of paper
#
# THE PROBLEM WITH Σ̃:
#   Affine equivariance means: if Y = AX + b, then T_Y = A T_X A^T.
#   Σ̃ fails this: for Y = cX (scalar scaling), Σ̃_Y = Σ̃_X unchanged.
#   This is a problem because dispersions should scale with data.
#
# THE FIX — Σ*:
#   Σ* is defined implicitly as the solution to:
#       Σ* = (p / VW(X)) * E[W^2(X) (X-mu)(X-mu)^T / (X-mu)^T Σ*^{-1} (X-mu)]
#   This is an M-estimator equation (like Huber's, Tyler's scatter).
#   The Mahalanobis distance (X-mu)^T Σ*^{-1} (X-mu) appears in denominator,
#   re-normalizing for scale — this restores affine equivariance.
#
# COMPUTATION:
#   No closed-form solution. We iterate:
#       Σ*^{(k+1)} = (p / VW) * (1/n) sum_i W^2(X_i) (X_i-mu)(X_i-mu)^T
#                    / (X_i-mu)^T [Σ*^{(k)}]^{-1} (X_i-mu)
#   until convergence (fixed-point iteration).
#
# THEOREM 5 (Asymptotic Variance):
#   ARE(γ̂*_i, γ̂_i; F) = V12^{-1}
#   where V12 = E[(alpha_{Σ*}(|Z|) * S12(Z;0))^2]
#   Table 2 shows ARE >> 1 for heavy-tailed distributions.
#
# TRADEOFF (Section 6):
#   Σ* is MORE efficient than Σ̃ but:
#   - Computationally more expensive (iterative)
#   - Less robust with high contamination
#   - Σ̃ with projection depth is recommended for general use
# =============================================================================

# Source dependencies
# source("01_weight_functions.R")
# source("02_weighted_spatial_median.R")
# source("03_weighted_sign_cov.R")


# =============================================================================
# AFFINE EQUIVARIANT SCATTER MATRIX Σ*
# =============================================================================
# Solves the fixed-point equation iteratively.
#
# The iteration is:
#   Σ*^{(k+1)} = (p / vw) * mean_i [ W_i^2 * outer(d_i, d_i) / md_i^{(k)} ]
# where:
#   d_i    = X_i - mu_hat  (centered observation)
#   md_i   = d_i^T [Σ*^{(k)}]^{-1} d_i  (squared Mahalanobis distance)
#   vw     = mean(W^2)  (normalizing constant = VW(X) in paper notation)
#   p      = dimension
#
# WHY (p / vw)?
#   The factor ensures E[Σ*] has trace p (standardization),
#   making the solution identifiable (otherwise Σ* and c*Σ* are both
#   solutions, so we pin scale by this normalization).
#
# INPUT:
#   X        : n x p data matrix
#   W        : n-vector of weights (from compute_weights)
#   mu_hat   : p-vector, location estimate
#   max_iter : maximum iterations
#   tol      : convergence tolerance (Frobenius norm of change)
#   init     : initial Sigma* (if NULL, uses WSCM from module 03)
#   verbose  : print iteration info
#
# OUTPUT: list with
#   Sigma_star : p x p affine equivariant scatter matrix
#   iter       : iterations until convergence
#   converged  : logical

affine_equivariant_scatter <- function(X, W, mu_hat, max_iter = 200,
                                       tol = 1e-7, init = NULL,
                                       verbose = FALSE) {

  X      <- as.matrix(X)
  n      <- nrow(X)
  p      <- ncol(X)
  W      <- as.numeric(W)

  # Centered observations: d_i = X_i - mu_hat
  D      <- sweep(X, 2, mu_hat, "-")   # n x p matrix of d_i

  # Normalizing constant: vw = mean(W^2)
  vw     <- mean(W^2)
  scale  <- p / vw                     # the (p / VW) factor from paper eq. (3)

  # Initial Sigma*: use WSCM as warm start (same eigenvectors, scale differs)
  if (is.null(init)) {
    Sigma_star <- weighted_sign_cov(X, W, mu_hat)
    # Make sure it's scaled to have trace p
    Sigma_star <- Sigma_star * p / sum(diag(Sigma_star))
  } else {
    Sigma_star <- init
  }

  converged <- FALSE

  for (iter in seq_len(max_iter)) {

    Sigma_old <- Sigma_star

    # Compute inverse of current Sigma*
    Sigma_inv <- tryCatch(
      solve(Sigma_star),
      error = function(e) {
        # If singular, use pseudoinverse via eigendecomposition
        eig <- eigen(Sigma_star, symmetric = TRUE)
        idx <- eig$values > .Machine$double.eps * max(eig$values)
        eig$vectors[, idx] %*%
          diag(1 / eig$values[idx]) %*%
          t(eig$vectors[, idx])
      }
    )

    # Squared Mahalanobis distances: md_i = d_i^T Sigma_inv d_i
    # Efficiently as rowSums((D %*% Sigma_inv) * D)
    md <- rowSums((D %*% Sigma_inv) * D)    # n-vector

    # Avoid division by zero
    md_safe <- pmax(md, .Machine$double.eps)

    # Update: Sigma*^{new} = scale * (1/n) sum_i W_i^2 outer(d_i,d_i) / md_i
    # Rewrite as: scale * t(D * (W^2 / md_safe)) %*% D / n
    effective_w <- W^2 / md_safe            # n-vector of effective weights
    Sigma_star  <- scale * t(D * effective_w) %*% D / n

    # Symmetrize
    Sigma_star <- (Sigma_star + t(Sigma_star)) / 2

    # Check convergence: Frobenius norm of change
    delta <- norm(Sigma_star - Sigma_old, type = "F")

    if (verbose) {
      cat(sprintf("  Iter %3d: Frobenius change = %.2e\n", iter, delta))
    }

    if (delta < tol) {
      converged <- TRUE
      break
    }
  }

  if (!converged) {
    warning(sprintf(
      "affine_equivariant_scatter: did not converge in %d iterations (delta=%.2e)",
      max_iter, delta))
  }

  return(list(Sigma_star = Sigma_star, iter = iter, converged = converged))
}


# =============================================================================
# FULL PIPELINE: affine equivariant robust PCA
# =============================================================================
# INPUT:
#   X      : n x p data matrix
#   method : weight method "HSD", "MhD", or "PD"
#   ...    : passed to affine_equivariant_scatter
#
# OUTPUT: list with
#   Sigma_star : p x p scatter matrix
#   Gamma_hat  : p x p eigenvector matrix
#   Lambda_hat : p-vector of eigenvalues
#   mu_hat     : location estimate
#   weights    : weights used

robust_pca_affine <- function(X, method = "PD", ...) {

  X <- as.matrix(X)

  # Step 1: Location estimate
  wsm_out <- compute_wsm(X, method = method)
  mu_hat  <- wsm_out$median
  W       <- wsm_out$weights

  # Step 2: Sigma*
  ae_out     <- affine_equivariant_scatter(X, W, mu_hat, ...)
  Sigma_star <- ae_out$Sigma_star

  # Step 3: Spectral decomposition
  eig        <- eigen(Sigma_star, symmetric = TRUE)
  Gamma_hat  <- eig$vectors
  Lambda_hat <- eig$values

  return(list(
    Sigma_star = Sigma_star,
    Gamma_hat  = Gamma_hat,
    Lambda_hat = Lambda_hat,
    mu_hat     = mu_hat,
    weights    = W,
    iter       = ae_out$iter,
    converged  = ae_out$converged
  ))
}


# =============================================================================
# TYLER'S M-ESTIMATOR (benchmark from Section 6)
# =============================================================================
# Tyler's scatter matrix is the special case where W^2(x) = p (constant),
# corresponding to the pure sign covariance approach normalized for scale.
# Its fixed-point equation is:
#   V = (p/n) sum_i (X_i-mu)(X_i-mu)^T / (X_i-mu)^T V^{-1} (X_i-mu)
#
# In our framework: vw = p (since W=sqrt(p) is constant => W^2=p, mean(W^2)=p)
# So scale = p/vw = p/p = 1.
#
# We implement it separately for clean comparison in simulations.
#
# INPUT:
#   X       : n x p data matrix
#   mu_hat  : p-vector, center estimate
#   max_iter, tol, verbose : as before
#
# OUTPUT: p x p Tyler's scatter matrix

tyler_scatter <- function(X, mu_hat, max_iter = 200, tol = 1e-7,
                          verbose = FALSE) {

  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  W <- rep(sqrt(p), n)    # constant weight W = sqrt(p) so W^2 = p

  out <- affine_equivariant_scatter(X, W, mu_hat, max_iter = max_iter,
                                    tol = tol, verbose = verbose)
  return(out$Sigma_star)
}


# =============================================================================
# QUICK DEMO
# =============================================================================
# Uncomment and run to test:

# library(MASS)
# set.seed(123)
# p <- 4
# Sigma_true <- diag(c(4, 3, 2, 1))
# X <- mvrnorm(n = 200, mu = rep(0, p), Sigma = Sigma_true)
#
# # Add 10% outliers
# X[1:20, ] <- matrix(rnorm(20 * p, mean = 10, sd = 0.5), 20, p)
#
# # Compare: affine equivariant Sigma* vs Sigma-tilde vs Tyler
# res_ae    <- robust_pca_affine(X, method = "PD", verbose = FALSE)
# res_wscm  <- robust_pca_wscm(X, method = "PD")
# mu_hat    <- colMeans(X)
# V_tyler   <- tyler_scatter(X, mu_hat)
# eig_tyler <- eigen(V_tyler, symmetric = TRUE)
#
# cat("True 1st eigenvector:   ", c(1, 0, 0, 0), "\n")
# cat("Sigma* (PD) 1st PC:     ", res_ae$Gamma_hat[, 1], "\n")
# cat("Sigma-tilde (PD) 1st PC:", res_wscm$Gamma_hat[, 1], "\n")
# cat("Tyler 1st PC:           ", eig_tyler$vectors[, 1], "\n")
# cat("Sample cov 1st PC:      ", eigen(cov(X))$vectors[, 1], "\n")
# cat(sprintf("Sigma* converged in %d iterations\n", res_ae$iter))
