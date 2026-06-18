# =============================================================================
# FILE: 06_influence_functions.R
# PAPER: "On Weighted Multivariate Sign Functions"
#        Majumdar & Chatterjee, Journal of Multivariate Analysis (2022)
#
# PURPOSE:
#   Computes and visualizes influence functions (IF) for eigenvectors of:
#     - Sample covariance matrix (unbounded IF = non-robust)
#     - Sign Covariance Matrix / SCM (bounded but inlier-effect problem)
#     - Tyler's scatter matrix (bounded but inlier-effect problem)
#     - Σ̃ with HSD, MhD, PD weights (bounded + no inlier effect)
#   Reproduces Figure 3 of the paper.
#
# SECTION: 5 of paper
#
# WHAT IS AN INFLUENCE FUNCTION?
#   IF(x0; T, F) = lim_{eps->0} [T((1-eps)F + eps*delta_x0) - T(F)] / eps
#   Measures how much a single contaminating point at x0 shifts the
#   functional T. If ||IF|| is bounded, T is robust. If unbounded, not robust.
#
# KEY RESULTS:
#
# Proposition 1 — IF of γ̃_i (ith eigenvector of Σ̃):
#   IF(x0; γ̃_i, F) = W^2(x0) * sum_{k≠i} S_{ik}(x0;mu) / (λ̃_i - λ̃_k) * gamma_k
#   where S_{ik}(x0;mu) = S_i(x0;mu) * S_k(x0;mu) is the (i,k) element
#   of the outer product sign matrix S(x0;mu)S(x0;mu)^T.
#
#   INTERPRETATION:
#   - The influence depends on W^2(x0) * off-diagonal elements of sign matrix.
#   - Since |S_{ik}| <= 1 always (bounded by unit sphere geometry) AND
#     W is bounded for HSD/MhD/PD, the IF is BOUNDED => ROBUST.
#   - Points near center (small |x0|) have small W^2(x0) => low influence
#     (no "inlier effect" unlike SCM and Tyler's).
#
# Proposition 3 — IF of γ*_i (ith eigenvector of Σ*):
#   IF(x0; γ*_i, F) = alpha_{Σ*}(|x0|) * sum_{k≠i} S_{ik}(x0;mu)/(λ_i-λ_k) * gamma_k
#   where alpha_{Σ*}(|x0|) = p(p+2)u(|z0|) / E[pu(|Z|) + u'(|Z|)]
#   Also bounded when W (and hence u) is bounded.
# =============================================================================

# Source dependencies
# source("01_weight_functions.R")
# source("02_weighted_spatial_median.R")
# source("03_weighted_sign_cov.R")


# =============================================================================
# INFLUENCE FUNCTION OF ith EIGENVECTOR OF Σ̃ (Proposition 1)
# =============================================================================
# IF(x0; γ̃_i, F) = W^2(x0) * sum_{k≠i} S_{ik}(x0;mu)/(λ̃_i-λ̃_k) * gamma_k
#
# INPUT:
#   x0         : p-vector, the contamination point
#   i          : which eigenvector (1 = first PC)
#   Gamma      : p x p eigenvector matrix (from Σ̃)
#   Lambda_tilde: p-vector of eigenvalues of Σ̃
#   mu         : p-vector, center
#   W_x0       : scalar weight W(x0) at the point x0
#
# OUTPUT: p-vector (the IF at x0)

IF_eigvec_Sigma_tilde <- function(x0, i, Gamma, Lambda_tilde, mu, W_x0) {

  x0 <- as.numeric(x0)
  mu <- as.numeric(mu)
  p  <- length(x0)

  # Sign vector at x0: s = (x0 - mu) / |x0 - mu|
  d    <- x0 - mu
  d_norm <- sqrt(sum(d^2))
  if (d_norm < .Machine$double.eps) return(rep(0, p))  # at center, IF = 0
  s    <- d / d_norm                                    # unit sign vector

  # Outer product sign matrix S(x0;mu) = s s^T
  S_mat <- outer(s, s)                                  # p x p

  # Sum over k != i:
  # sum_{k≠i} S_{ik}(x0;mu) / (λ̃_i - λ̃_k) * gamma_k
  IF_val <- rep(0, p)
  for (k in seq_len(p)) {
    if (k == i) next
    denom <- Lambda_tilde[i] - Lambda_tilde[k]
    if (abs(denom) < .Machine$double.eps) next    # skip if eigenvalues equal

    S_ik  <- S_mat[i, k]                           # (i,k) element of sign matrix
    IF_val <- IF_val + (S_ik / denom) * Gamma[, k]
  }

  return(W_x0^2 * IF_val)
}


# =============================================================================
# INFLUENCE FUNCTION OF SAMPLE COVARIANCE EIGENVECTOR (benchmark, eq. 2)
# =============================================================================
# AV(sqrt(n) * gamma_hat_i) = sum_{k≠i} lambda_i*lambda_k/(lambda_i-lambda_k)^2
#                              * gamma_k gamma_k^T
#
# The influence function itself:
# IF(x0; gamma_i, F) = sum_{k≠i} [(x0-mu)_i * (x0-mu)_k - Sigma_{ik}]
#                      / (lambda_i - lambda_k) * gamma_k
# where subscripts _i, _k refer to the ith and kth coordinates in PC space.
#
# INPUT:
#   x0      : p-vector
#   i       : eigenvector index
#   Gamma   : p x p eigenvectors of Sigma
#   Lambda  : p-vector of true eigenvalues
#   mu      : p-vector center
#   Sigma   : p x p covariance matrix
#
# OUTPUT: p-vector IF at x0

IF_eigvec_SampleCov <- function(x0, i, Gamma, Lambda, mu, Sigma) {

  x0 <- as.numeric(x0)
  mu <- as.numeric(mu)
  p  <- length(x0)

  d    <- x0 - mu                       # centered point
  z    <- t(Gamma) %*% d               # in PC coordinates (p-vector)

  IF_val <- rep(0, p)
  for (k in seq_len(p)) {
    if (k == i) next
    denom <- Lambda[i] - Lambda[k]
    if (abs(denom) < .Machine$double.eps) next

    # z[i]*z[k] - Sigma_{ik} (in PC space, Sigma_{ik} = 0 for i≠k)
    contrib <- (z[i] * z[k]) / denom
    IF_val  <- IF_val + contrib * Gamma[, k]
  }

  return(IF_val)
}


# =============================================================================
# INFLUENCE FUNCTION OF SCM EIGENVECTOR
# =============================================================================
# IF(x0; gamma_{S,i}, F) = sum_{k≠i} S_{ik}(x0;mu) / (lambda_{S,i} - lambda_{S,k})
#                          * gamma_k
# where lambda_{S,i} = E_Z[lambda_i * z_i^2 / sum_j lambda_j * z_j^2]
# (the eigenvalues of the SCM are different from those of Sigma)
#
# For a spherical distribution with Sigma = diag(lambda_1,...,lambda_p):
#   lambda_{S,i} ~ lambda_i / sum_j lambda_j  (approximately proportional)
# In general, this requires numerical computation.
#
# For simplicity in plotting (Fig 3), we use the approximate form with
# lambda_{S,i} replaced by the empirical eigenvalues of the SCM.

IF_eigvec_SCM <- function(x0, i, Gamma_scm, Lambda_scm, mu) {

  x0 <- as.numeric(x0)
  mu <- as.numeric(mu)
  p  <- length(x0)

  d      <- x0 - mu
  d_norm <- sqrt(sum(d^2))
  if (d_norm < .Machine$double.eps) return(rep(0, p))
  s      <- d / d_norm
  S_mat  <- outer(s, s)

  IF_val <- rep(0, p)
  for (k in seq_len(p)) {
    if (k == i) next
    denom <- Lambda_scm[i] - Lambda_scm[k]
    if (abs(denom) < .Machine$double.eps) next

    S_ik  <- S_mat[i, k]
    IF_val <- IF_val + (S_ik / denom) * Gamma_scm[, k]
  }

  return(IF_val)
}


# =============================================================================
# INFLUENCE FUNCTION OF TYLER'S SCATTER EIGENVECTOR
# =============================================================================
# IF(x0; gamma_{T,i}, F) = (p+2) * sum_{k≠i} S_{ik}(x0;mu)/(lambda_i-lambda_k)
#                          * gamma_k

IF_eigvec_Tyler <- function(x0, i, Gamma, Lambda, mu) {

  x0 <- as.numeric(x0)
  mu <- as.numeric(mu)
  p  <- length(x0)

  d      <- x0 - mu
  d_norm <- sqrt(sum(d^2))
  if (d_norm < .Machine$double.eps) return(rep(0, p))
  s      <- d / d_norm
  S_mat  <- outer(s, s)

  IF_val <- rep(0, p)
  for (k in seq_len(p)) {
    if (k == i) next
    denom <- Lambda[i] - Lambda[k]
    if (abs(denom) < .Machine$double.eps) next

    S_ik  <- S_mat[i, k]
    IF_val <- IF_val + (S_ik / denom) * Gamma[, k]
  }

  return((p + 2) * IF_val)
}


# =============================================================================
# REPRODUCE FIGURE 3: IF norm plots over a grid of x0 points
# =============================================================================
# For bivariate data (p=2): plot ||IF(x0)|| as a surface over a grid of x0.
# The paper plots this for the first eigenvector (i=1) of each method.
#
# INPUT:
#   Sigma   : 2x2 true covariance matrix (paper uses diag(2,1))
#   mu      : 2-vector center (paper uses c(0,0))
#   n_grid  : grid resolution
#   method  : weight method for Sigma-tilde
#   n_sim   : number of simulated samples for estimating quantities

plot_IF_norms <- function(Sigma = diag(c(2, 1)), mu = c(0, 0),
                          n_grid = 30, method = "MhD", n_sim = 5000) {

  library(MASS)

  # Generate data for estimating distributions
  X <- mvrnorm(n = n_sim, mu = mu, Sigma = Sigma)

  # True spectral decomposition
  eig_true   <- eigen(Sigma, symmetric = TRUE)
  Gamma_true <- eig_true$vectors
  Lambda_true <- eig_true$values

  # Estimated Sigma-tilde and SCM
  wsm       <- compute_wsm(X, method = method)
  pca_wscm  <- robust_pca_wscm(X, method = method, mu_hat = wsm$median)
  scm       <- sign_covariance_matrix(X, wsm$median)
  eig_scm   <- eigen(scm, symmetric = TRUE)

  # Grid of contamination points
  grid_seq  <- seq(-3, 3, length.out = n_grid)
  grid      <- expand.grid(x1 = grid_seq, x2 = grid_seq)

  # Compute ||IF(x0)|| for each method at each grid point
  n_pts <- nrow(grid)
  IF_norms <- data.frame(
    x1         = grid$x1,
    x2         = grid$x2,
    SampleCov  = numeric(n_pts),
    SCM        = numeric(n_pts),
    Tyler      = numeric(n_pts),
    WSCM       = numeric(n_pts)
  )

  for (j in seq_len(n_pts)) {
    x0 <- c(grid$x1[j], grid$x2[j])

    # Weight at x0 for WSCM
    Z_x0   <- standardize_to_Z(t(x0), mu = wsm$median, Sigma = cov(X))
    W_x0   <- switch(method,
                     "HSD" = W_HSD(Z_x0),
                     "MhD" = W_MhD(Z_x0),
                     "PD"  = W_PD(Z_x0))

    # Sample cov IF
    if_sc <- IF_eigvec_SampleCov(x0, 1, Gamma_true, Lambda_true, mu, Sigma)
    IF_norms$SampleCov[j] <- sqrt(sum(if_sc^2))

    # SCM IF
    if_scm <- IF_eigvec_SCM(x0, 1, eig_scm$vectors, eig_scm$values, wsm$median)
    IF_norms$SCM[j] <- sqrt(sum(if_scm^2))

    # Tyler IF
    if_ty <- IF_eigvec_Tyler(x0, 1, Gamma_true, Lambda_true, wsm$median)
    IF_norms$Tyler[j] <- sqrt(sum(if_ty^2))

    # WSCM IF
    if_wscm <- IF_eigvec_Sigma_tilde(x0, 1,
                                     pca_wscm$Gamma_hat,
                                     pca_wscm$Lambda_hat,
                                     wsm$median, W_x0)
    IF_norms$WSCM[j] <- sqrt(sum(if_wscm^2))
  }

  # Plotting as perspective plots (like Fig 3)
  par(mfrow = c(2, 2), mar = c(2, 2, 3, 1))

  plot_persp <- function(vals, title, zlim = NULL) {
    M <- matrix(vals, n_grid, n_grid)
    if (is.null(zlim)) zlim <- range(M, finite = TRUE)
    persp(grid_seq, grid_seq, M,
          main = title, xlab = "x1", ylab = "x2", zlab = "||IF||",
          theta = 30, phi = 20, col = "lightblue",
          zlim = zlim, ticktype = "detailed")
  }

  plot_persp(IF_norms$SampleCov, "(a) Sample Covariance")
  plot_persp(IF_norms$SCM,       "(b) SCM")
  plot_persp(IF_norms$Tyler,     "(c) Tyler")
  plot_persp(IF_norms$WSCM,      paste0("(d) Sigma-tilde (", method, ")"))

  invisible(IF_norms)
}


# =============================================================================
# QUICK DEMO
# =============================================================================
# Uncomment to reproduce a simplified version of Figure 3:

# # Requires sourcing all prior modules
# par(mfrow=c(2,2))
# IF_data <- plot_IF_norms(
#   Sigma  = diag(c(2, 1)),
#   mu     = c(0, 0),
#   n_grid = 25,
#   method = "MhD",
#   n_sim  = 2000
# )
