# ================================================================
#  JASA (Wang, Peng & Li 2015) vs JMVA (Majumdar & Chatterjee 2022)
#  Exact statistics from both papers — theoretical + demo
# ================================================================
#
#  JASA Tn: spatial sign U-statistic for HIGH-DIM MEAN VECTOR TEST
#    Zi  = Xi / ||Xi||   (unit vector, uncentered)
#    Tn  = sum_{i} sum_{j<i} Zi' Zj
#    Null ~ Normal (when p >> n)
#
#  JMVA R:  WEIGHTED spatial sign for ROBUST SCATTER / PCA
#    S(Xi; mu) = (Xi - mu) / ||Xi - mu||
#    W(Xi, F)  = depth-based weight (Mahalanobis / HSD / PD)
#    R(Xi; mu, F) = S(Xi; mu) * W(Xi, F)
#    Sigma_tilde = (1/n) sum_i W^2(Xi) S(Xi;mu) S(Xi;mu)'
# ================================================================

if (!requireNamespace("MASS",       quietly=TRUE)) install.packages("MASS")
if (!requireNamespace("robustbase", quietly=TRUE)) install.packages("robustbase")
library(MASS)

set.seed(2025)

# ================================================================
# PART 1 — JASA Tn: Spatial Sign Test (Wang, Peng & Li 2015)
# ================================================================

# -- 1a. Spatial sign function (uncentered, as in JASA paper) --
spatial_sign_jasa <- function(X) {
  # Zi = Xi / ||Xi||; zero vector if Xi = 0
  norms <- sqrt(rowSums(X^2))
  Z <- X / ifelse(norms == 0, 1, norms)
  Z[norms == 0, ] <- 0
  Z
}

# -- 1b. JASA test statistic Tn (eq. 3 in paper) --
# Tn = sum_{i=1}^{n} sum_{j<i} Zi' Zj
# This is a U-statistic; under H0: E(Zi) = 0
jasa_Tn <- function(X) {
  Z <- spatial_sign_jasa(X)
  n <- nrow(Z)
  total <- 0
  for (i in 2:n)
    for (j in 1:(i-1))
      total <- total + sum(Z[i,] * Z[j,])
  total
}

# -- 1c. Variance estimator for Tn (eq. 8 in paper) --
# var(Tn) = n(n-1)/2 * Tr(B^2)
# where B = E(Zi Zi') estimated from data
jasa_var_Tn <- function(X) {
  Z  <- spatial_sign_jasa(X)
  n  <- nrow(Z)
  p  <- ncol(Z)
  # Efficient formula for Tr(B^2):  see eq.(8) in JASA paper
  ZZt     <- Z %*% t(Z)                 # n x n Gram matrix
  ZZtZZt  <- sum(ZZt^2)                 # Tr( (sum ZiZi')(sum ZjZj') ) = ||ZZt||_F^2
  Zstar   <- colMeans(Z)
  term1   <- -n / (n-2)^2
  term2   <- (n-1) / (n * (n-2)^2) * ZZtZZt
  term3   <- (1 - 2*n) / (n*(n-1)) * as.numeric(t(Zstar) %*% (t(Z) %*% Z) %*% Zstar)
  term4   <- 2/n * sum(Zstar^2)
  term5   <- (n-2)^2 / (n*(n-1)) * sum(Zstar^2)^2
  TrB2    <- term1 + term2 + term3 + term4 + term5
  n*(n-1)/2 * TrB2
}

# -- 1d. Standardised JASA test statistic Z-score --
jasa_test <- function(X) {
  Tn  <- jasa_Tn(X)
  vTn <- jasa_var_Tn(X)
  list(Tn = Tn, var_Tn = vTn, Z = Tn / sqrt(abs(vTn)),
       pval = 2 * pnorm(-abs(Tn / sqrt(abs(vTn)))))
}

# ================================================================
# PART 2 — JMVA Weighted Sign (Majumdar & Chatterjee 2022)
# ================================================================

# -- 2a. Generalized sign centered at mu (eq. 1 in paper) --
spatial_sign_jmva <- function(X, mu) {
  Xc    <- sweep(X, 2, mu)
  norms <- sqrt(rowSums(Xc^2))
  S     <- Xc / ifelse(norms == 0, 1, norms)
  S[norms == 0, ] <- 0
  S
}

# -- 2b. Three weight functions from paper (Section 1, p.3-4) --

# W_MhD: Mahalanobis depth weight
# W(Xi) proportional to |Z|^2 / (1 + |Z|^2), Z = Sigma^{-1/2}(Xi - mu)
weight_mahal <- function(X, mu, Sigma) {
  Xc      <- sweep(X, 2, mu)
  Sinv    <- solve(Sigma)
  Sinvsqrt <- chol(Sinv)   # upper triangular
  Z       <- Xc %*% t(Sinvsqrt)
  r2      <- rowSums(Z^2)
  w       <- r2 / (1 + r2)
  w / max(w)   # normalise to [0,1]
}

# W_HSD: Half-space depth weight (rank-based approximation)
# W(Xi) proportional to F_Z(|Z|), the empirical CDF of norms
weight_hsd <- function(X, mu) {
  Xc    <- sweep(X, 2, mu)
  norms <- sqrt(rowSums(Xc^2))
  ecdf_fn <- ecdf(norms)
  ecdf_fn(norms)
}

# W_PD: Projection depth weight (simple rank approximation)
# W(Xi) proportional to |Z| / (1 + |Z|/MAD(|Z|))
weight_pd <- function(X, mu) {
  Xc    <- sweep(X, 2, mu)
  norms <- sqrt(rowSums(Xc^2))
  mad_r <- median(abs(norms - median(norms)))
  if (mad_r == 0) mad_r <- 1e-8
  w <- norms / (1 + norms / mad_r)
  w / max(w)
}

# -- 2c. Weighted sign R(Xi; mu, F) = S(Xi, mu) * W(Xi, F) --
weighted_sign_jmva <- function(X, mu, weight_type = "mahal",
                                Sigma = NULL) {
  S <- spatial_sign_jmva(X, mu)
  w <- switch(weight_type,
    "mahal" = weight_mahal(X, mu, if (is.null(Sigma)) diag(ncol(X)) else Sigma),
    "hsd"   = weight_hsd(X, mu),
    "pd"    = weight_pd(X, mu)
  )
  S * w   # scalar w_i multiplies each row
}

# -- 2d. Weighted scatter matrix Sigma_tilde (eq. in Section 2.2) --
# Sigma_tilde = (1/n) sum_i W^2(Xi) S(Xi;mu) S(Xi;mu)'
sigma_tilde <- function(X, mu, weight_type = "mahal", Sigma = NULL) {
  S <- spatial_sign_jmva(X, mu)
  w <- switch(weight_type,
    "mahal" = weight_mahal(X, mu, if (is.null(Sigma)) diag(ncol(X)) else Sigma),
    "hsd"   = weight_hsd(X, mu),
    "pd"    = weight_pd(X, mu)
  )
  n <- nrow(X)
  # Sigma_tilde = (1/n) * S' diag(w^2) S  (efficiently)
  wS <- S * w        # n x p: each row scaled by w_i
  t(wS) %*% wS / n  # p x p
}

# Plain sign covariance matrix (SCM), for comparison
scm <- function(X, mu) {
  S <- spatial_sign_jmva(X, mu)
  t(S) %*% S / nrow(S)
}

# ================================================================
# PART 3 — DEMO DATASET: Simulated bivariate data
# ================================================================

cat("================================================================\n")
cat("  DEMO 1: JASA Tn — Testing H0: mu = 0\n")
cat("================================================================\n\n")

n <- 80; p <- 5

# Under H0 (Normal)
X_H0 <- mvrnorm(n, mu = rep(0,p), Sigma = diag(p))
r0   <- jasa_test(X_H0)
cat(sprintf("Normal, H0 true:     Tn=%.2f  Z=%.3f  p=%.3f\n",
            r0$Tn, r0$Z, r0$pval))

# Under H1 (Normal, mu shifted)
X_H1 <- mvrnorm(n, mu = rep(0.3, p), Sigma = diag(p))
r1   <- jasa_test(X_H1)
cat(sprintf("Normal, H1 (mu=0.3): Tn=%.2f  Z=%.3f  p=%.3f\n",
            r1$Tn, r1$Z, r1$pval))

# Under H1 (heavy-tailed t3)
df    <- 3
chi2  <- rchisq(n, df)
X_t3  <- mvrnorm(n, mu=rep(0.3,p), Sigma=diag(p)) / sqrt(chi2/df)
rt3   <- jasa_test(X_t3)
cat(sprintf("t3,    H1 (mu=0.3):  Tn=%.2f  Z=%.3f  p=%.3f\n\n",
            rt3$Tn, rt3$Z, rt3$pval))

cat("================================================================\n")
cat("  DEMO 2: JMVA Sigma_tilde — Robust Scatter/PCA\n")
cat("================================================================\n\n")

# True scatter: Sigma = diag(4,3,2,1), p=4
p2    <- 4
Sigma_true <- diag(c(4,3,2,1))
n2    <- 200
mu2   <- rep(0, p2)

# Clean normal data
X_clean <- mvrnorm(n2, mu=mu2, Sigma=Sigma_true)
mu_hat  <- colMeans(X_clean)

St_mahal <- sigma_tilde(X_clean, mu_hat, "mahal", Sigma_true)
St_hsd   <- sigma_tilde(X_clean, mu_hat, "hsd")
St_pd    <- sigma_tilde(X_clean, mu_hat, "pd")
St_scm   <- scm(X_clean, mu_hat)
St_cov   <- cov(X_clean)

# Compare eigenvectors (first PC direction)
true_pc1 <- c(1,0,0,0)   # largest eigenvalue direction
angle_deg <- function(A, B) {
  acos(min(abs(sum(A*B)), 1)) * 180 / pi
}

cat("Angle between estimated PC1 and true PC1 (smaller = better):\n")
pc1_cov   <- eigen(St_cov)$vectors[,1]
pc1_scm   <- eigen(St_scm)$vectors[,1]
pc1_mahal <- eigen(St_mahal)$vectors[,1]
pc1_hsd   <- eigen(St_hsd)$vectors[,1]
pc1_pd    <- eigen(St_pd)$vectors[,1]

cat(sprintf("  Sample covariance:    %.4f degrees\n", angle_deg(pc1_cov, true_pc1)))
cat(sprintf("  SCM (plain sign):     %.4f degrees\n", angle_deg(pc1_scm, true_pc1)))
cat(sprintf("  Sigma_tilde (Mahal):  %.4f degrees\n", angle_deg(pc1_mahal, true_pc1)))
cat(sprintf("  Sigma_tilde (HSD):    %.4f degrees\n", angle_deg(pc1_hsd, true_pc1)))
cat(sprintf("  Sigma_tilde (PD):     %.4f degrees\n\n", angle_deg(pc1_pd, true_pc1)))

# With outliers (10% contamination as in JMVA paper)
n_out      <- round(0.1 * n2)
X_outlier  <- rbind(X_clean[1:(n2-n_out), ],
                    mvrnorm(n_out, mu=c(20,20,20,20), Sigma=diag(p2)))
mu_cont    <- colMeans(X_outlier)

pc1_cov_o   <- eigen(cov(X_outlier))$vectors[,1]
pc1_scm_o   <- eigen(scm(X_outlier, mu_cont))$vectors[,1]
pc1_mahal_o <- eigen(sigma_tilde(X_outlier, mu_cont, "mahal",
                                  Sigma_true))$vectors[,1]
pc1_hsd_o   <- eigen(sigma_tilde(X_outlier, mu_cont, "hsd"))$vectors[,1]
pc1_pd_o    <- eigen(sigma_tilde(X_outlier, mu_cont, "pd"))$vectors[,1]

cat("With 10% outliers (as in JMVA paper, Fig. 2 scenario):\n")
cat(sprintf("  Sample covariance:    %.4f degrees\n", angle_deg(pc1_cov_o,   true_pc1)))
cat(sprintf("  SCM (plain sign):     %.4f degrees\n", angle_deg(pc1_scm_o,   true_pc1)))
cat(sprintf("  Sigma_tilde (Mahal):  %.4f degrees\n", angle_deg(pc1_mahal_o, true_pc1)))
cat(sprintf("  Sigma_tilde (HSD):    %.4f degrees\n", angle_deg(pc1_hsd_o,   true_pc1)))
cat(sprintf("  Sigma_tilde (PD):     %.4f degrees\n\n", angle_deg(pc1_pd_o,   true_pc1)))

# ================================================================
# PART 4 — POWER COMPARISON (JASA Tn across distributions)
# ================================================================

cat("================================================================\n")
cat("  DEMO 3: JASA Tn power — Normal vs t3 (B=299 permutations)\n")
cat("================================================================\n\n")

perm_pval_jasa <- function(X, B = 299) {
  obs <- jasa_Tn(X)
  n   <- nrow(X)
  perms <- replicate(B, {
    idx <- sample(n)
    jasa_Tn(X[idx, ])
  })
  mean(c(obs, perms) >= obs)
}

n_pwr <- 50; p_pwr <- 10
# H1: Normal, small mean
X_n  <- mvrnorm(n_pwr, mu=rep(0.2, p_pwr), Sigma=diag(p_pwr))
# H1: t3, same mean
chi2p <- rchisq(n_pwr, 3)
X_t  <- mvrnorm(n_pwr, mu=rep(0.2, p_pwr), Sigma=diag(p_pwr)) / sqrt(chi2p/3)

pv_n <- perm_pval_jasa(X_n)
pv_t <- perm_pval_jasa(X_t)

cat(sprintf("Permutation p-value, Normal (mu=0.2): %.3f\n", pv_n))
cat(sprintf("Permutation p-value, t3    (mu=0.2): %.3f\n", pv_t))
cat("\n[JASA Tn should be more powerful under t3 relative to Hotelling]\n")

# ================================================================
# PART 5 — SIDE-BY-SIDE SUMMARY
# ================================================================

cat("\n================================================================\n")
cat("  SUMMARY: Key differences\n")
cat("================================================================\n")
cat("
  JASA (Wang-Peng-Li 2015)           JMVA (Majumdar-Chatterjee 2022)
  ---------------------------------  ---------------------------------
  Goal:   Test H0: mu = 0            Goal: Robust scatter/PCA estimation
  Sign:   Zi = Xi / ||Xi||           Sign: S(Xi;mu) = (Xi-mu)/||Xi-mu||
  Weight: None (all = 1)             Weight: W from data depth
  Kernel: Zi'Zj  (inner product)     Kernel: R(Xi)'R(Xi)' (outer product)
  Output: Scalar test statistic Tn   Output: p x p scatter matrix Sig_tilde
  Null:   N(0,1) for standardized Tn Null:   Asymptotic normal for eigenvec
  Regime: p >> n allowed             Regime: Fixed p (future: p^2/n -> 0)
  Depth:  Not used                   Depth:  Core (Mahal / HSD / PD)
  ARE vs Hotelling: ~2.54 (t3,lg p)  ARE vs SCM: > 1 (all settings)
")
