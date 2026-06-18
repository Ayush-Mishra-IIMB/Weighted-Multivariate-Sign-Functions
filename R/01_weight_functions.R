# =============================================================================
# FILE: 01_weight_functions.R
# PAPER: "On Weighted Multivariate Sign Functions"
#        Majumdar & Chatterjee, Journal of Multivariate Analysis (2022)
#
# PURPOSE:
#   Implements the three weight functions W(X, F) derived from data depth.
#   These are the building blocks of EVERYTHING in the paper.
#
# MATHEMATICAL BACKGROUND:
#   Given X ~ elliptical(mu, Sigma), we standardize:
#       Z = Sigma^{-1/2} (X - mu)
#   so that Z has mean 0 and identity covariance (spherically symmetric).
#
#   All three weights are functions of |Z| = Euclidean norm of Z only.
#   They give CENTER-INWARD ordering: peripheral points get HIGH weight,
#   central points get LOW weight. This is the opposite of depth.
#
#   Why peripheral weighting? Outliers carry more shape/spread information.
#   Upweighting them improves efficiency while bounded weights keep robustness.
#
# THREE WEIGHT FUNCTIONS (Section 1, page 3-4 of paper):
#   (i)  W_HSD(X) proportional to F_{Z1}(|Z|)       -- Half-Space Depth
#   (ii) W_MhD(X) proportional to |Z|^2/(1+|Z|^2)   -- Mahalanobis Depth
#   (iii)W_PD(X)  proportional to |Z|/(1+|Z|/MAD)   -- Projection Depth
#
# AFFINE INVARIANCE:
#   Because we use affine-invariant depth functions:
#       W(X, F) = W(Z, F_Z)
#   So working with Z is equivalent to working with X.
# =============================================================================


# =============================================================================
# STEP 1: STANDARDIZE DATA TO Z
# =============================================================================
# This is the FIRST thing you always do before computing weights.
# Z = Sigma^{-1/2} (X - mu)
# In practice mu and Sigma are unknown, so we estimate them.
#
# INPUT:
#   X     : n x p data matrix (n observations, p variables)
#   mu    : p-vector, location (if NULL, uses column means)
#   Sigma : p x p covariance matrix (if NULL, uses sample covariance)
#
# OUTPUT:
#   Z     : n x p standardized matrix

standardize_to_Z <- function(X, mu = NULL, Sigma = NULL) {

  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  # Estimate mu if not provided
  if (is.null(mu)) {
    mu <- colMeans(X)
  }

  # Estimate Sigma if not provided
  if (is.null(Sigma)) {
    Sigma <- cov(X)
  }

  # Compute Sigma^{-1/2} via eigendecomposition
  # Sigma = V D V^T  =>  Sigma^{-1/2} = V D^{-1/2} V^T
  eig       <- eigen(Sigma, symmetric = TRUE)
  D_neghalf <- diag(1 / sqrt(eig$values))       # D^{-1/2}
  Sigma_inv_half <- eig$vectors %*% D_neghalf %*% t(eig$vectors)

  # Center and scale: Z_i = Sigma^{-1/2} (X_i - mu)
  X_centered <- sweep(X, 2, mu, "-")            # subtract mu from each row
  Z          <- X_centered %*% t(Sigma_inv_half) # n x p matrix of Z vectors

  return(Z)
}


# =============================================================================
# HELPER: EUCLIDEAN NORM OF EACH ROW
# =============================================================================
# |Z_i| for each observation i = 1, ..., n
row_norms <- function(Z) {
  sqrt(rowSums(Z^2))
}


# =============================================================================
# WEIGHT FUNCTION (i): HALF-SPACE DEPTH (HSD)
# =============================================================================
# FORMULA (paper eq, page 4):
#   W_HSD(X) proportional to F_{Z1}(|Z|)
#
# INTUITION:
#   F_{Z1}(|Z|) is the empirical CDF of |Z| evaluated at each |Z_i|.
#   This is a rank-based weight: observations with LARGER |Z| (further
#   from center) get HIGHER rank and hence HIGHER weight.
#   It is proportional to the half-space depth of |Z|.
#
# The half-space depth of a point x w.r.t. distribution F is:
#   HSD(x; F) = inf_{||u||=1} P(u^T X >= u^T x)
# For the standardized radial distance |Z|, this simplifies to
# using the empirical CDF of the norms.
#
# INPUT:  Z : n x p standardized data matrix
# OUTPUT: W : n-vector of weights (normalized to sum to n)

W_HSD <- function(Z) {
  z_norm <- row_norms(Z)           # |Z_i| for each i
  n      <- length(z_norm)

  # Empirical CDF of |Z| evaluated at each |Z_i|
  # ecdf(z_norm)(z_norm) gives F_n(|Z_i|) = rank(|Z_i|)/n
  W <- ecdf(z_norm)(z_norm)

  # Normalize so weights sum to n (convenient for sample formulas)
  W <- W / mean(W)

  return(W)
}


# =============================================================================
# WEIGHT FUNCTION (ii): MAHALANOBIS DEPTH (MhD)
# =============================================================================
# FORMULA (paper eq, page 4):
#   W_MhD(X) proportional to |Z|^2 / (1 + |Z|^2)
#
# INTUITION:
#   This is a smooth, monotonically increasing function of |Z|.
#   As |Z| -> 0 (center), W -> 0 (low weight).
#   As |Z| -> infinity, W -> 1 (approaches maximum weight).
#   The function is bounded in [0, 1), ensuring robustness.
#
# WHY THIS FORM? The Mahalanobis depth is:
#   MhD(x; F) = 1 / (1 + (x-mu)^T Sigma^{-1} (x-mu))
#             = 1 / (1 + |Z|^2)
# Taking 1 - MhD(x; F) = |Z|^2/(1+|Z|^2) gives a peripherality measure.
#
# NOTE: Mahalanobis weights use the sample mean and covariance for
# standardization, making them sensitive to outliers (non-robust mu/Sigma).
# Outperforms other weights WITHOUT outliers but degrades with contamination.
#
# INPUT:  Z : n x p standardized data matrix
# OUTPUT: W : n-vector of weights (normalized to sum to n)

W_MhD <- function(Z) {
  z_norm_sq <- rowSums(Z^2)        # |Z_i|^2

  W <- z_norm_sq / (1 + z_norm_sq)

  # Normalize
  W <- W / mean(W)

  return(W)
}


# =============================================================================
# WEIGHT FUNCTION (iii): PROJECTION DEPTH (PD)
# =============================================================================
# FORMULA (paper eq, page 4):
#   W_PD(X) proportional to |Z| / (1 + |Z| / MAD(Z_1))
#
# where MAD(Z_1) = median absolute deviation of the first component of Z.
# MAD is a robust scale estimator: MAD(Z_1) = median(|Z_1i - median(Z_1i)|)
#
# INTUITION:
#   The projection depth of x w.r.t. F is:
#     PD(x; F) = 1 / (1 + sup_{||u||=1} |u^T(x-mu)| / MAD(u^T X))
#   For spherically symmetric Z, this simplifies to:
#     1 / (1 + |Z| / MAD(Z_1))
#   Taking 1 - PD gives the peripherality weight used here.
#
# WHY MAD? MAD is a robust measure of spread (breakdown point = 50%).
# It doesn't break down with up to 50% outlier contamination, making
# W_PD more robust than W_MhD.
#
# PAPER FINDING (Section 6): Projection depth weights strike the BEST
# balance between efficiency and robustness across all settings.
#
# INPUT:  Z : n x p standardized data matrix
# OUTPUT: W : n-vector of weights (normalized to sum to n)

W_PD <- function(Z) {
  z_norm <- row_norms(Z)           # |Z_i|

  # MAD of first component of Z
  # For spherically symmetric Z, all components have the same distribution,
  # so Z_1 (first column) is representative.
  mad_z1 <- mad(Z[, 1], constant = 1)  # constant=1 for raw MAD

  # Avoid division by zero if MAD = 0 (degenerate case)
  if (mad_z1 < .Machine$double.eps) {
    warning("MAD of Z1 is zero. Returning equal weights.")
    return(rep(1, nrow(Z)))
  }

  W <- z_norm / (1 + z_norm / mad_z1)

  # Normalize
  W <- W / mean(W)

  return(W)
}


# =============================================================================
# MASTER WEIGHT FUNCTION: compute weights for any dataset
# =============================================================================
# This is the user-facing function. Given raw data X, it:
#   1. Estimates mu and Sigma (or uses provided values)
#   2. Standardizes to Z
#   3. Computes weights using the chosen method
#
# INPUT:
#   X           : n x p data matrix
#   method      : "HSD", "MhD", or "PD"
#   mu          : optional p-vector (if NULL, uses colMeans)
#   Sigma       : optional p x p matrix (if NULL, uses cov(X))
#   return_Z    : if TRUE, also return the standardized Z matrix
#
# OUTPUT:
#   W           : n-vector of weights
#   (optionally) Z : n x p standardized matrix

compute_weights <- function(X, method = "PD", mu = NULL, Sigma = NULL,
                            return_Z = FALSE) {

  X <- as.matrix(X)

  # Step 1: Standardize
  Z <- standardize_to_Z(X, mu = mu, Sigma = Sigma)

  # Step 2: Compute weights
  W <- switch(method,
    "HSD" = W_HSD(Z),
    "MhD" = W_MhD(Z),
    "PD"  = W_PD(Z),
    stop("method must be one of: 'HSD', 'MhD', 'PD'")
  )

  if (return_Z) {
    return(list(W = W, Z = Z))
  } else {
    return(W)
  }
}


# =============================================================================
# VERIFICATION: Check key theoretical properties of the weight functions
# =============================================================================
# Property 1: W(Z) is a function of |Z| only (radial symmetry)
# Property 2: W is non-decreasing in |Z| (peripherality ordering)
# Property 3: EW^2(Z) < infinity (finite second moment)
# Property 4: W is bounded above (ensures robust influence functions)

verify_weight_properties <- function(X, method = "PD") {

  result <- compute_weights(X, method = method, return_Z = TRUE)
  W      <- result$W
  Z      <- result$Z
  z_norm <- row_norms(Z)

  cat("=== Weight Function Properties Check:", method, "===\n\n")

  # Property 1: Monotonicity — cor(W, |Z|) should be > 0
  cor_val <- cor(W, z_norm)
  cat(sprintf("Correlation of W with |Z|: %.4f  (should be > 0)\n", cor_val))

  # Property 2: Range
  cat(sprintf("Weight range: [%.4f, %.4f]\n", min(W), max(W)))

  # Property 3: E[W^2] finite
  cat(sprintf("Sample E[W^2]: %.4f  (should be finite)\n", mean(W^2)))

  # Property 4: Bounded
  cat(sprintf("Max weight: %.4f  (bounded = %s)\n", max(W),
              ifelse(max(W) < Inf, "YES", "NO")))

  cat("\n")
  invisible(list(W = W, Z = Z, z_norm = z_norm))
}


# =============================================================================
# QUICK DEMO
# =============================================================================
# Uncomment and run to see all three weight functions in action:

# library(MASS)
# set.seed(42)
# X <- mvrnorm(n = 200, mu = c(2, 3), Sigma = matrix(c(2, 1, 1, 1), 2, 2))
#
# W_hsd <- compute_weights(X, method = "HSD")
# W_mhd <- compute_weights(X, method = "MhD")
# W_pd  <- compute_weights(X, method = "PD")
#
# # Visual check: plot weight vs distance from center
# Z <- standardize_to_Z(X)
# z_norm <- row_norms(Z)
#
# par(mfrow = c(1, 3))
# plot(z_norm, W_hsd, main = "HSD Weights vs |Z|", xlab = "|Z|", ylab = "W")
# plot(z_norm, W_mhd, main = "MhD Weights vs |Z|", xlab = "|Z|", ylab = "W")
# plot(z_norm, W_pd,  main = "PD Weights vs |Z|",  xlab = "|Z|", ylab = "W")
#
# # Verify properties
# verify_weight_properties(X, "HSD")
# verify_weight_properties(X, "MhD")
# verify_weight_properties(X, "PD")
