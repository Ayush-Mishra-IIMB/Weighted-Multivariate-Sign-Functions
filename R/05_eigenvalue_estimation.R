# =============================================================================
# FILE: 05_eigenvalue_estimation.R
# PAPER: "On Weighted Multivariate Sign Functions"
#        Majumdar & Chatterjee, Journal of Multivariate Analysis (2022)
#
# PURPOSE:
#   Estimates population eigenvalues lambda_i of Sigma robustly.
#   Σ̃ gives correct eigenvectors but WRONG eigenvalues (Theorem 2).
#   This module corrects that.
#
# SECTION: 4 of paper
#
# THE PROBLEM:
#   Theorem 2 shows: Σ̃ = Gamma Lambda-tilde Gamma^T
#   where lambda-tilde_i = lambda_i * E[W^2(X)] * E[U_i^2 / U^T Lambda U]
#   So the eigenvalues of Σ̃ are NOT the population eigenvalues lambda_i.
#   They are scaled/distorted versions of them.
#
# THE SOLUTION — Median-of-Small-Variances Estimator:
#   1. Compute Σ̂̃ and its eigenvectors Gamma-hat (correct directions).
#   2. Rotate data: S = Gamma-hat^T X  (project onto estimated PCs)
#   3. Divide S into k disjoint groups of size ~n/k.
#   4. For each group j and coordinate i, compute variance lambda_{i,j}^dagger.
#   5. Take coordinate-wise MEDIANS across groups: lambda_i^dagger = median_j(...)
#   6. Plug in: Sigma^dagger = Gamma-hat * diag(lambda^dagger) * Gamma-hat^T
#
# WHY MEDIAN OF VARIANCES?
#   - Within each small group, coordinate-wise variance is a consistent
#     estimator of lambda_i (the true eigenvalue).
#   - Taking the median across groups makes it robust:
#     even if some groups are contaminated by outliers, the majority
#     are clean, and the median is not affected.
#   - THEOREM 6: Sigma^dagger is a CONSISTENT estimator of Sigma
#     as n -> inf and k -> inf with n/k -> inf.
#
# CHOICE OF k:
#   Paper says: follow [34] (Minsker 2015 geometric median paper).
#   Practical rule: k ~ log(n) or k ~ sqrt(n) works well.
#   More groups = more robust; fewer groups = more efficient per group.
# =============================================================================

# Source dependencies
# source("01_weight_functions.R")
# source("02_weighted_spatial_median.R")
# source("03_weighted_sign_cov.R")


# =============================================================================
# MEDIAN-OF-SMALL-VARIANCES EIGENVALUE ESTIMATOR
# =============================================================================
# INPUT:
#   X         : n x p data matrix
#   Gamma_hat : p x p matrix of eigenvectors from Σ̂̃ (columns = eigenvectors)
#   k         : number of groups (default: floor(log(n)))
#   center    : if TRUE, subtract group mean within each group (recommended)
#
# OUTPUT: list with
#   lambda_dag : p-vector of robust eigenvalue estimates
#   Sigma_dag  : p x p consistent scatter estimate
#   k_used     : actual k used

estimate_eigenvalues_robust <- function(X, Gamma_hat, k = NULL, center = TRUE) {

  X         <- as.matrix(X)
  n         <- nrow(X)
  p         <- ncol(X)

  # Default k: log(n) is a practical rule from Minsker (2015)
  if (is.null(k)) {
    k <- max(2, floor(log(n)))
  }

  # Ensure k <= n/2 (each group needs at least 2 observations)
  k <- min(k, floor(n / 2))

  # Step 1: Rotate data to estimated PC coordinates
  # S = Gamma_hat^T X  =>  S_il = gamma_l^T X_i (projection onto lth PC)
  S <- X %*% Gamma_hat    # n x p  (each row is one observation in PC space)

  # Step 2: Randomly divide indices into k groups of size floor(n/k)
  group_size <- floor(n / k)
  idx        <- sample(seq_len(n))      # random permutation
  groups     <- vector("list", k)
  for (j in seq_len(k)) {
    groups[[j]] <- idx[((j-1)*group_size + 1) : (j * group_size)]
  }
  # Remaining indices (if n not divisible by k) go to last group
  remainder <- idx[(k * group_size + 1) : n]
  if (length(remainder) > 0) {
    groups[[k]] <- c(groups[[k]], remainder)
  }

  # Step 3: Coordinate-wise variance for each group
  # lambda_{i,j}^dagger = (1/|G_j|) sum_{l in G_j} (S_{li} - S̄_{G_j,i})^2
  lambda_mat <- matrix(0, nrow = k, ncol = p)   # k x p matrix of variances

  for (j in seq_len(k)) {
    S_j <- S[groups[[j]], , drop = FALSE]        # group j data in PC space

    if (center) {
      S_j_mean   <- colMeans(S_j)
      S_j_center <- sweep(S_j, 2, S_j_mean, "-")
    } else {
      S_j_center <- S_j
    }

    # Coordinate-wise variance (each column = one PC direction)
    lambda_mat[j, ] <- colMeans(S_j_center^2)
  }

  # Step 4: Coordinate-wise MEDIANS across groups
  lambda_dag <- apply(lambda_mat, 2, median)     # p-vector

  # Step 5: Reconstruct Sigma^dagger = Gamma_hat * diag(lambda_dag) * Gamma_hat^T
  Sigma_dag <- Gamma_hat %*% diag(lambda_dag) %*% t(Gamma_hat)

  # Symmetrize for numerical stability
  Sigma_dag <- (Sigma_dag + t(Sigma_dag)) / 2

  return(list(
    lambda_dag = lambda_dag,
    Sigma_dag  = Sigma_dag,
    k_used     = k,
    lambda_mat = lambda_mat   # k x p table of group variances (for diagnostics)
  ))
}


# =============================================================================
# FULL PIPELINE: robust covariance matrix estimation
# =============================================================================
# Combines WSCM eigenvectors with robust eigenvalue estimation.
# This gives a consistent estimate of the full covariance matrix Sigma.
#
# INPUT:
#   X      : n x p data matrix
#   method : weight method "HSD", "MhD", or "PD"
#   k      : number of groups (NULL = automatic)
#
# OUTPUT: list with
#   Sigma_dag   : p x p robust covariance estimate (consistent for Sigma)
#   Gamma_hat   : p x p eigenvectors (from Σ̃)
#   lambda_dag  : p-vector of robust eigenvalue estimates
#   Lambda_tilde: p-vector of eigenvalues of Σ̃ (for comparison)
#   mu_hat      : location estimate

robust_covariance <- function(X, method = "PD", k = NULL) {

  X <- as.matrix(X)

  # Step 1: Robust PCA via WSCM (gets correct eigenvectors)
  pca_out    <- robust_pca_wscm(X, method = method)
  Gamma_hat  <- pca_out$Gamma_hat

  # Step 2: Robust eigenvalue estimation
  eig_out    <- estimate_eigenvalues_robust(X, Gamma_hat, k = k)

  return(list(
    Sigma_dag    = eig_out$Sigma_dag,
    Gamma_hat    = Gamma_hat,
    lambda_dag   = eig_out$lambda_dag,
    Lambda_tilde = pca_out$Lambda_hat,   # eigenvalues of Σ̃ (not of Sigma)
    mu_hat       = pca_out$mu_hat,
    k_used       = eig_out$k_used
  ))
}


# =============================================================================
# COMPARISON: eigenvalue estimation accuracy
# =============================================================================
# Compare estimated vs true eigenvalues (when true values are known,
# e.g., in simulations).
#
# INPUT:
#   lambda_true : p-vector of true eigenvalues
#   lambda_dag  : p-vector of estimated eigenvalues
#   lambda_tilde: p-vector of raw Σ̃ eigenvalues (for comparison)

compare_eigenvalues <- function(lambda_true, lambda_dag, lambda_tilde = NULL) {

  p <- length(lambda_true)

  cat("=== Eigenvalue Estimation Comparison ===\n\n")
  cat(sprintf("%-5s  %-12s  %-12s", "PC", "True", "lambda_dag"))
  if (!is.null(lambda_tilde)) cat(sprintf("  %-12s", "lambda_tilde"))
  cat("\n")
  cat(strrep("-", 50), "\n")

  for (i in seq_len(p)) {
    cat(sprintf("%-5d  %-12.4f  %-12.4f", i, lambda_true[i], lambda_dag[i]))
    if (!is.null(lambda_tilde)) cat(sprintf("  %-12.4f", lambda_tilde[i]))
    cat("\n")
  }

  # Relative errors
  cat("\nRelative error of lambda_dag: ",
      round(abs(lambda_dag - lambda_true) / lambda_true, 4), "\n")
  if (!is.null(lambda_tilde)) {
    cat("Relative error of lambda_tilde:",
        round(abs(lambda_tilde - lambda_true) / lambda_true, 4), "\n")
  }
}


# =============================================================================
# QUICK DEMO
# =============================================================================
# Uncomment and run to test:

# library(MASS)
# set.seed(99)
# p <- 4
# lambda_true <- c(4, 3, 2, 1)
# Sigma_true  <- diag(lambda_true)
# X <- mvrnorm(n = 500, mu = rep(0, p), Sigma = Sigma_true)
#
# # Robust covariance estimation
# cov_result <- robust_covariance(X, method = "PD")
#
# compare_eigenvalues(
#   lambda_true  = lambda_true,
#   lambda_dag   = cov_result$lambda_dag,
#   lambda_tilde = cov_result$Lambda_tilde
# )
#
# # Frobenius norm error
# cat(sprintf("\n||Sigma_dag - Sigma_true||_F = %.4f\n",
#             norm(cov_result$Sigma_dag - Sigma_true, type = "F")))
# cat(sprintf("||Sigma_hat - Sigma_true||_F = %.4f  (sample cov)\n",
#             norm(cov(X) - Sigma_true, type = "F")))
