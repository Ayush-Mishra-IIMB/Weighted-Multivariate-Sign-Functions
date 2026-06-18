# =============================================================================
# FILE: 03_weighted_sign_cov.R
# PAPER: "On Weighted Multivariate Sign Functions"
#        Majumdar & Chatterjee, Journal of Multivariate Analysis (2022)
#
# PURPOSE:
#   Implements Σ̃ (Sigma-tilde) — the paper's core contribution:
#   the Weighted Sign Covariance Matrix (WSCM).
#   Used for robust PCA (eigenvectors of Σ).
#
# SECTION: 2.2 and 2.3 of paper
#
# KEY MATH:
#   Define the weighted sign surrogate:
#       X̃ = W(X, F) * S(X; mu)
#   where S(X; mu) = (X-mu)/|X-mu| is the generalized sign function.
#
#   The population WSCM is:
#       Σ̃ = E[X̃ X̃^T] = E[W^2(X) S(X;mu) S(X;mu)^T]
#
#   THEOREM 2: The eigenvectors of Σ̃ are IDENTICAL to those of Σ!
#   So we can use Σ̃ to robustly estimate the DIRECTIONS (eigenvectors)
#   of the true covariance matrix, even though the eigenvalues differ.
#
# WHY DOES THIS WORK?
#   X = mu + R * Gamma * Lambda^{1/2} * U
#   where U ~ uniform on unit sphere, R ~ radial, independent of U.
#   Then S(X; mu) = Gamma * Lambda^{1/2} * Z / |Lambda^{1/2} Z|
#   Since W(X) = W(|Z|) is a function of R = |Z| only, and R is
#   independent of U, the weight W factors out of the expectation,
#   leaving the same eigenvector structure as the population covariance.
#
# SAMPLE VERSION (Theorem 3 — Asymptotic Normality):
#   Σ̂̃ = (1/n) sum_i W^2(X_i, Fn) S(X_i; mu_hat) S(X_i; mu_hat)^T
#   sqrt(n) [vec(Σ̂̃) - E W^2(X) vec(S(X;mu))] --D--> N(0, VW)
#
# CONNECTION TO EXISTING METHODS:
#   W = 1:                Σ̃ = Sign Covariance Matrix (SCM) [Visuri et al. 2000]
#   W = |X-mu|:           Σ̃ = sample covariance matrix (up to scale)
#   Our W (depth-based):  improves efficiency over SCM while keeping
#                         robustness (bounded influence function)
# =============================================================================

# Source dependencies
# source("01_weight_functions.R")
# source("02_weighted_spatial_median.R")


# =============================================================================
# OUTER PRODUCT SIGN MATRIX S(x; mu)
# =============================================================================
# S(x; mu) = S(x;mu) S(x;mu)^T = (x-mu)(x-mu)^T / |x-mu|^2
#
# This is a p x p rank-1 matrix for each observation.
# The sample WSCM averages these over all observations.
#
# INPUT:
#   X   : n x p data matrix
#   mu  : p-vector, center estimate
#
# OUTPUT:
#   List of n matrices, each p x p. (Or directly sum them.)

sign_outer_products <- function(X, mu) {
  X  <- as.matrix(X)
  n  <- nrow(X)
  p  <- ncol(X)

  diffs <- sweep(X, 2, mu, "-")      # X_i - mu
  dists <- sqrt(rowSums(diffs^2))    # |X_i - mu|
  nonzero <- dists > .Machine$double.eps

  # Initialize p x p accumulator
  S_list <- vector("list", n)
  for (i in seq_len(n)) {
    if (nonzero[i]) {
      s_i      <- diffs[i, ] / dists[i]       # unit vector (sign)
      S_list[[i]] <- outer(s_i, s_i)           # p x p outer product
    } else {
      S_list[[i]] <- matrix(0, p, p)           # zero matrix at center
    }
  }
  return(S_list)
}


# =============================================================================
# WEIGHTED SIGN COVARIANCE MATRIX: Σ̂̃
# =============================================================================
# Σ̂̃ = (1/n) sum_i W^2(X_i) S(X_i; mu_hat)
#    = (1/n) sum_i W^2_i * (X_i - mu_hat)(X_i - mu_hat)^T / |X_i - mu_hat|^2
#
# This is the SAMPLE version of Theorem 2's population parameter Σ̃.
# The eigenvectors of Σ̂̃ consistently estimate those of Σ.
#
# INPUT:
#   X       : n x p data matrix
#   W       : n-vector of weights W(X_i) (NOT squared — we square inside)
#   mu_hat  : p-vector, location estimate (use WSM from module 02)
#
# OUTPUT:
#   Sigma_tilde : p x p symmetric matrix

weighted_sign_cov <- function(X, W, mu_hat) {

  X      <- as.matrix(X)
  n      <- nrow(X)
  p      <- ncol(X)
  W      <- as.numeric(W)

  diffs  <- sweep(X, 2, mu_hat, "-")   # X_i - mu_hat  (n x p)
  dists  <- sqrt(rowSums(diffs^2))     # |X_i - mu_hat| (n-vector)

  nonzero <- dists > .Machine$double.eps

  # Accumulate: sum_i W_i^2 * outer(s_i, s_i)
  Sigma_tilde <- matrix(0, p, p)
  for (i in which(nonzero)) {
    s_i         <- diffs[i, ] / dists[i]          # unit sign vector
    Sigma_tilde <- Sigma_tilde + W[i]^2 * outer(s_i, s_i)
  }
  Sigma_tilde <- Sigma_tilde / n

  # Symmetrize (numerical safety)
  Sigma_tilde <- (Sigma_tilde + t(Sigma_tilde)) / 2

  return(Sigma_tilde)
}


# =============================================================================
# SIGN COVARIANCE MATRIX (SCM): unweighted benchmark
# =============================================================================
# SCM = (1/n) sum_i S(X_i; mu) S(X_i; mu)^T
#     = (1/n) sum_i (X_i-mu)(X_i-mu)^T / |X_i-mu|^2
#
# Special case of WSCM with W = 1.
# Used as benchmark in paper's simulations (Fig 2).
#
# INPUT:
#   X      : n x p data matrix
#   mu_hat : p-vector, center estimate
#
# OUTPUT:
#   SCM : p x p matrix

sign_covariance_matrix <- function(X, mu_hat) {
  n <- nrow(X)
  W_ones <- rep(1, n)
  return(weighted_sign_cov(X, W_ones, mu_hat))
}


# =============================================================================
# FULL PIPELINE: robust scatter + PCA
# =============================================================================
# Given raw data X:
#   1. Compute robust location mu_hat (WSM from module 02)
#   2. Compute weights W
#   3. Compute Σ̂̃
#   4. Compute spectral decomposition: Σ̂̃ = Γ̂ Λ̂̃ Γ̂^T
#   5. Return eigenvectors (robust PCA directions) and eigenvalues
#
# INPUT:
#   X      : n x p data matrix
#   method : weight method "HSD", "MhD", or "PD"
#   mu_hat : optional pre-computed location (if NULL, uses WSM)
#
# OUTPUT: list with
#   Sigma_tilde : p x p WSCM
#   Gamma_hat   : p x p matrix of eigenvectors (columns = PCs)
#   Lambda_hat  : p-vector of eigenvalues (in decreasing order)
#   mu_hat      : p-vector, location estimate used
#   weights     : n-vector, weights used

robust_pca_wscm <- function(X, method = "PD", mu_hat = NULL) {

  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  # Step 1: Location estimation
  if (is.null(mu_hat)) {
    wsm_out <- compute_wsm(X, method = method)
    mu_hat  <- wsm_out$median
    W       <- wsm_out$weights
  } else {
    # Use provided mu_hat to compute weights
    Z <- standardize_to_Z(X, mu = mu_hat)
    W <- switch(method,
                "HSD" = W_HSD(Z),
                "MhD" = W_MhD(Z),
                "PD"  = W_PD(Z))
  }

  # Step 2: Weighted Sign Covariance Matrix
  Sigma_tilde <- weighted_sign_cov(X, W, mu_hat)

  # Step 3: Spectral decomposition
  eig         <- eigen(Sigma_tilde, symmetric = TRUE)
  Gamma_hat   <- eig$vectors      # columns = eigenvectors
  Lambda_hat  <- eig$values       # eigenvalues (decreasing order by default)

  return(list(
    Sigma_tilde = Sigma_tilde,
    Gamma_hat   = Gamma_hat,
    Lambda_hat  = Lambda_hat,
    mu_hat      = mu_hat,
    weights     = W
  ))
}


# =============================================================================
# ASYMPTOTIC VARIANCE VW (Appendix B of paper)
# =============================================================================
# VW = V[W^2(X) vec(S(X;mu))]
#    = E[W^4(X) vec(S(X;mu)) vec(S(X;mu))^T] - vec(Σ̃) vec(Σ̃)^T
#
# In practice, we estimate this from the sample.
# This is needed for Theorem 4 (eigenvector asymptotics) and
# confidence intervals for principal components.
#
# INPUT:
#   X       : n x p data matrix
#   W       : n-vector of weights
#   mu_hat  : p-vector, location estimate
#
# OUTPUT: p^2 x p^2 matrix VW

estimate_VW <- function(X, W, mu_hat) {

  X  <- as.matrix(X)
  n  <- nrow(X)
  p  <- ncol(X)
  p2 <- p^2

  diffs  <- sweep(X, 2, mu_hat, "-")
  dists  <- sqrt(rowSums(diffs^2))
  nonzero <- dists > .Machine$double.eps

  # Sample average of W^4 * vec(S_i) vec(S_i)^T
  VW_raw <- matrix(0, p2, p2)
  for (i in which(nonzero)) {
    s_i     <- diffs[i, ] / dists[i]
    S_i     <- outer(s_i, s_i)           # p x p
    vs_i    <- as.vector(S_i)            # p^2 vector
    VW_raw  <- VW_raw + W[i]^4 * outer(vs_i, vs_i)
  }
  VW_raw <- VW_raw / n

  # Subtract vec(Σ̃) vec(Σ̃)^T
  Sigma_tilde <- weighted_sign_cov(X, W, mu_hat)
  vs_tilde    <- as.vector(Sigma_tilde)
  VW          <- VW_raw - outer(vs_tilde, vs_tilde)

  return(VW)
}


# =============================================================================
# PRINT SUMMARY OF ROBUST PCA
# =============================================================================

summarize_robust_pca <- function(pca_result, X) {

  p   <- ncol(X)
  cat("=== Robust PCA via Weighted Sign Covariance Matrix ===\n\n")
  cat(sprintf("Dimension p = %d,  n = %d\n\n", p, nrow(X)))

  cat("Location estimate (mu_hat):\n")
  cat(sprintf("  %s\n\n", paste(round(pca_result$mu_hat, 4), collapse = ", ")))

  cat("Eigenvalues of Σ̃ (proportional to PC variances):\n")
  Lambda <- pca_result$Lambda_hat
  pct    <- Lambda / sum(Lambda) * 100
  for (i in seq_len(p)) {
    cat(sprintf("  PC%d: lambda = %.4f  (%.1f%% variance)\n", i, Lambda[i], pct[i]))
  }

  cat("\nFirst eigenvector (1st PC direction):\n")
  cat(sprintf("  %s\n", paste(round(pca_result$Gamma_hat[, 1], 4), collapse = ", ")))
}


# =============================================================================
# QUICK DEMO
# =============================================================================
# Uncomment and run to test:

# library(MASS)
# set.seed(42)
# # Generate data with known covariance structure
# Sigma_true <- matrix(c(4,2,2,1), 2, 2)   # eigenvalues: 5 and 0; PC1 = (2,1)/sqrt(5)
# X <- mvrnorm(n = 300, mu = c(0, 0), Sigma = Sigma_true)
#
# # Add 10% outliers
# n_out <- 30
# X[1:n_out, ] <- mvrnorm(n_out, mu = c(10, 10), Sigma = diag(2))
#
# # Robust PCA
# res_pd  <- robust_pca_wscm(X, method = "PD")
# res_mhd <- robust_pca_wscm(X, method = "MhD")
# res_scm <- list(
#   Sigma_tilde = sign_covariance_matrix(X, colMeans(X)),
#   Gamma_hat   = eigen(sign_covariance_matrix(X, colMeans(X)))$vectors
# )
#
# cat("True 1st PC:    ", c(2,1)/sqrt(5), "\n")
# cat("PD 1st PC:      ", res_pd$Gamma_hat[,1], "\n")
# cat("MhD 1st PC:     ", res_mhd$Gamma_hat[,1], "\n")
# cat("SCM 1st PC:     ", eigen(res_scm$Sigma_tilde)$vectors[,1], "\n")
# cat("Sample cov PC:  ", eigen(cov(X))$vectors[,1], "\n")
