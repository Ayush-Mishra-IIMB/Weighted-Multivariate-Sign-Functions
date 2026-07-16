# ================================================================
#  SINGLE FILE: JASA + JMVA + TRUE MEAN COMPARISON
#  Wang, Peng & Li (2015) JASA  vs  Majumdar & Chatterjee (2022) JMVA
#
#  HOW TO RUN:  source("full_comparison_single_file.R")
#  RUNTIME:     ~25-35 min (B=500 for JASA, B=1000 for JMVA)
#  OUTPUT:      Console tables matching paper's Tables 1,2,3 and Fig 2
#
#  STRUCTURE:
#   SECTION 0 — packages + seed
#   SECTION 1 — all shared functions (signs, weights, statistics)
#   SECTION 2 — JASA simulation  (Tables 1, 2, 3 replication)
#   SECTION 3 — JMVA simulation  (Fig 2 FSE replication)
#   SECTION 4 — TRUE MEAN comparison (both papers, same data)
# ================================================================

# ================================================================
# SECTION 0: Setup
# ================================================================

if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")
library(MASS)

set.seed(2025)   # fix seed for reproducibility

# ================================================================
# SECTION 1: Shared functions — used by ALL three sections below
# ================================================================

# ------------------------------------------------------------------
# 1A. JASA functions  (Wang, Peng & Li 2015)
# ------------------------------------------------------------------

# Spatial sign: Zi = Xi / ||Xi||  (uncentered, as in JASA eq. 2.1)
jasa_sign <- function(X) {
  nr <- sqrt(rowSums(X^2))
  Z  <- X / ifelse(nr == 0, 1, nr)
  Z[nr == 0, ] <- 0
  Z
}

# JASA test statistic Tn (fast formula, eq. 8 in paper)
# Identity: sum_{i} sum_{j<i} Zi'Zj = (||sum Zi||^2 - n) / 2
# Valid because ||Zi|| = 1 for all i
jasa_Tn <- function(X) {
  Z  <- jasa_sign(X)
  cs <- colSums(Z)
  (sum(cs^2) - nrow(Z)) / 2
}

# JASA variance estimator: var(Tn) = n(n-1)/2 * Tr(B^2)
# Using the efficient cross-validation form (eq. 8)
jasa_var <- function(X) {
  Z     <- jasa_sign(X)
  n     <- nrow(Z)
  ZZt   <- Z %*% t(Z)             # n x n Gram matrix
  ZtZ   <- t(Z) %*% Z             # p x p
  Zs    <- colMeans(Z)            # p-vector: Z-bar
  t1 <- -n / (n - 2)^2
  t2 <- (n - 1) / (n * (n - 2)^2) * sum(ZZt^2)
  t3 <- (1 - 2*n) / (n*(n-1))    * as.numeric(t(Zs) %*% ZtZ %*% Zs)
  t4 <- 2/n                       * sum(Zs^2)
  t5 <- (n-2)^2 / (n*(n-1))      * sum(Zs^2)^2
  n*(n-1)/2 * (t1 + t2 + t3 + t4 + t5)
}

# JASA standardised statistic and two-sided p-value
jasa_pval <- function(X) {
  Tn <- jasa_Tn(X)
  vn <- jasa_var(X)
  Zs <- Tn / sqrt(abs(vn))
  2 * pnorm(-abs(Zs))
}

# Chen-Qin (CQ 2010) test  — the paper's main competitor
# CQ = sum_{i != j} Xi' Xj  (off-diagonal sum, raw observations)
cq_pval <- function(X) {
  n    <- nrow(X)
  cq   <- sum(colSums(X)^2) - sum(rowSums(X^2))
  S2   <- var(X)
  trS2 <- sum(S2^2)
  var_cq <- 2 * n * (n - 1) * trS2
  2 * pnorm(-abs(cq / sqrt(abs(var_cq))))
}

# ------------------------------------------------------------------
# 1B. JMVA functions  (Majumdar & Chatterjee 2022)
# ------------------------------------------------------------------

# Generalised sign centred at mu: S(Xi; mu) = (Xi-mu)/|Xi-mu|
jmva_sign <- function(X, mu) {
  Xc <- sweep(X, 2, mu)
  nr <- sqrt(rowSums(Xc^2))
  S  <- Xc / ifelse(nr == 0, 1, nr)
  S[nr == 0, ] <- 0
  S
}

# Weight function W_MhD (Mahalanobis depth, paper p.3)
# W(Xi) proportional to |Z|^2 / (1 + |Z|^2),  Z = Sigma^{-1/2}(Xi - mu)
w_mahal <- function(X, mu, Sigma) {
  Xc <- sweep(X, 2, mu)
  Sq <- chol(solve(Sigma))        # upper Cholesky of Sigma^{-1}
  Z  <- Xc %*% t(Sq)
  r2 <- rowSums(Z^2)
  w  <- r2 / (1 + r2)
  w / max(w)
}

# Weight function W_HSD (half-space depth, paper p.3)
# W(Xi) proportional to F_Z(|Z|), the empirical CDF of centred norms
w_hsd <- function(X, mu) {
  nr <- sqrt(rowSums(sweep(X, 2, mu)^2))
  ecdf(nr)(nr)
}

# Weight function W_PD (projection depth, paper p.4)
# W(Xi) proportional to |Z| / (1 + |Z|/MAD(|Z|))
w_pd <- function(X, mu) {
  nr    <- sqrt(rowSums(sweep(X, 2, mu)^2))
  mad_r <- median(abs(nr - median(nr)))
  if (mad_r < 1e-8) mad_r <- 1e-8
  w <- nr / (1 + nr / mad_r)
  w / max(w)
}

# Weighted sign scatter matrix:
# Sigma_tilde = (1/n) sum_i W^2(Xi) S(Xi;mu) S(Xi;mu)'   (paper eq. 2.2)
sigma_tilde <- function(X, mu, wtype, Sigma_known = NULL) {
  S <- jmva_sign(X, mu)
  w <- switch(wtype,
    "mahal" = w_mahal(X, mu, if (is.null(Sigma_known)) diag(ncol(X)) else Sigma_known),
    "hsd"   = w_hsd(X, mu),
    "pd"    = w_pd(X, mu)
  )
  wS <- S * w                     # scale each row by its weight
  t(wS) %*% wS / nrow(X)
}

# Plain sign covariance matrix (SCM): W ≡ 1
scm <- function(X, mu) {
  S <- jmva_sign(X, mu)
  t(S) %*% S / nrow(S)
}

# ------------------------------------------------------------------
# 1C. Data generators used across all sections
# ------------------------------------------------------------------

# JASA distributions  (centred at mu, covariance Sigma)
gen_jasa <- function(n, mu, Sigma, dist) {
  p <- length(mu)
  if (dist == "normal") {
    return(mvrnorm(n, mu = mu, Sigma = Sigma))
  }
  if (dist == "t3") {
    chi2 <- rchisq(n, 3)
    Z    <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma) / sqrt(chi2 / 3)
    return(sweep(Z, 2, mu, "+"))
  }
  if (dist == "mix") {           # 0.9*N(mu,Sigma) + 0.1*N(mu,9*Sigma)
    idx <- sample(1:2, n, replace = TRUE, prob = c(0.9, 0.1))
    X   <- matrix(0, n, p)
    n1  <- sum(idx == 1); n2 <- sum(idx == 2)
    if (n1 > 0) X[idx == 1, ] <- mvrnorm(n1, mu = mu, Sigma = Sigma)
    if (n2 > 0) X[idx == 2, ] <- mvrnorm(n2, mu = mu, Sigma = 9 * Sigma)
    return(X)
  }
}

# JASA covariance structures (paper Section 3.1)
make_Sigma <- function(p, type) {
  if (type == 1) {               # compound symmetry, off-diag = 0.2
    S <- matrix(0.2, p, p); diag(S) <- 1; return(S)
  }
  if (type == 2) {               # AR(1), rho = 0.8
    return(outer(1:p, 1:p, function(i, j) 0.8^abs(i - j)))
  }
  if (type == 3) {               # Srivastava et al. (2013) setting
    d <- 2 + (p - 1:p + 1) / p
    R <- outer(1:p, 1:p, function(i, j)
          ifelse(i == j, 1, (-1)^(i+j) * 0.2^(abs(i-j) / 0.1)))
    D <- diag(d)
    return(D %*% R %*% D)
  }
}

# JASA mean vectors (paper Section 3.1)
make_mu <- function(p, type) {
  if (type == 0) return(rep(0, p))        # null
  if (type == 1) return(rep(0.25, p))     # dense alternative
  if (type == 2) {                         # mixed alternative
    mu <- rep(0, p)
    mu[(floor(p/3)+1):floor(2*p/3)] <-  0.25
    mu[(floor(2*p/3)+1):p]          <- -0.25
    return(mu)
  }
}

# JMVA: multivariate t generator
gen_t <- function(n, p, df, Sigma) {
  chi2 <- rchisq(n, df)
  mvrnorm(n, mu = rep(0, p), Sigma = Sigma) / sqrt(chi2 / df)
}

# JMVA metric: prediction angle between estimated and true PC1
pred_angle <- function(est_vec, true_vec) {
  acos(min(abs(sum(est_vec * true_vec)), 1))
}

# ================================================================
# SECTION 2: JASA Replication — Tables 1, 2, 3
# ================================================================
cat("\n")
cat("================================================================\n")
cat("SECTION 2  |  JASA Wang-Peng-Li (2015)\n")
cat("           |  Empirical size (mu0) and power (mu1, mu2)\n")
cat("           |  TRUE MEAN used throughout — no estimation error\n")
cat("================================================================\n")
cat("Paper settings: n=20,50  p=1000,2000  B=1000 (we use B=500)\n")
cat("Sigma1: compound sym  Sigma2: AR(1)  Sigma3: Srivastava\n")
cat("mu0=(0,...,0)  mu1=(0.25,...,0.25)  mu2=(0,+0.25,-0.25 thirds)\n\n")

B_jasa <- 500    # paper uses 1000; increase for closer match

run_jasa_cell <- function(n, p, mu_type, sigma_type, dist, B = B_jasa) {
  mu    <- make_mu(p, mu_type)      # TRUE mean (known analytically)
  Sigma <- make_Sigma(p, sigma_type)
  rej_new <- rej_cq <- 0
  for (b in 1:B) {
    X <- gen_jasa(n, mu, Sigma, dist)
    if (jasa_pval(X) < 0.05) rej_new <- rej_new + 1
    if (cq_pval(X)   < 0.05) rej_cq  <- rej_cq  + 1
  }
  c(New = round(rej_new/B, 3), CQ = round(rej_cq/B, 3))
}

print_jasa_table <- function(dist_code, dist_label) {
  cat(sprintf("\n--- %s ---\n", dist_label))
  cat(sprintf("%-4s %-4s %5s %5s | %6s %6s | Note\n",
              "Sig","mu","n","p","New","CQ"))
  cat(rep("-", 58), "\n", sep = "")
  for (sig in 1:2) {
    for (mu_t in 0:2) {
      for (n in c(20, 50)) {
        p   <- 1000
        res <- run_jasa_cell(n, p, mu_t, sig, dist_code)
        note <- if (mu_t == 0) {
          ok <- abs(res["New"] - 0.05) < 0.025
          sprintf("SIZE  -> target ~0.05  [%s]", ifelse(ok, "OK", "CHECK"))
        } else {
          sprintf("POWER -> should exceed 0.05")
        }
        cat(sprintf("Sig%d mu%d %5d %5d | %6.3f %6.3f | %s\n",
                    sig, mu_t, n, p, res["New"], res["CQ"], note))
      }
    }
    cat("\n")
  }
}

print_jasa_table("normal", "Example 1: Multivariate Normal        [Table 1]")
print_jasa_table("t3",     "Example 2: Multivariate t3 (heavy)    [Table 2]")
print_jasa_table("mix",    "Example 3: Scale mixture 0.9N+0.1N(9S)[Table 3]")

cat("\nJASA key claim to verify:\n")
cat("  Table 2 (t3): New >> CQ for all settings\n")
cat("  Table 1 (Normal): New ≈ CQ (no loss under normality)\n")

# ================================================================
# SECTION 3: JMVA Replication — Figure 2 (FSE)
# ================================================================
cat("\n")
cat("================================================================\n")
cat("SECTION 3  |  JMVA Majumdar-Chatterjee (2022)\n")
cat("           |  Finite sample efficiency (FSE) of first eigenvector\n")
cat("           |  TRUE Sigma=diag(4,3,2,1), TRUE PC1=(1,0,0,0)\n")
cat("================================================================\n")
cat("Paper settings: p=4  n in {50,100,...,500}  B=10000 (we use 1000)\n")
cat("Methods: CovMatrix(base), SCM, Sigma~-PD, Sigma~-HSD, Sigma~-Mah\n")
cat("FSE = MSPA(CovMatrix) / MSPA(method)  ->  FSE > 1 beats CovMatrix\n\n")

# Paper's exact JMVA settings
p_jmva     <- 4
Sigma_jmva <- diag(c(4, 3, 2, 1))  # TRUE scatter (known)
true_PC1   <- c(1, 0, 0, 0)         # TRUE first eigenvector (known)
B_jmva     <- 1000                   # paper: 10000

run_mspa <- function(gen_fn, true_vec, B, Sigma_true = NULL) {
  a_cov <- a_scm <- a_pd <- a_hsd <- a_mah <- numeric(B)
  for (b in 1:B) {
    X      <- gen_fn()
    mu_hat <- colMeans(X)
    a_cov[b] <- pred_angle(eigen(cov(X))$vectors[, 1],                       true_vec)
    a_scm[b] <- pred_angle(eigen(scm(X, mu_hat))$vectors[, 1],               true_vec)
    a_pd[b]  <- pred_angle(eigen(sigma_tilde(X, mu_hat, "pd"))$vectors[,1],  true_vec)
    a_hsd[b] <- pred_angle(eigen(sigma_tilde(X, mu_hat, "hsd"))$vectors[,1], true_vec)
    a_mah[b] <- pred_angle(
                  eigen(sigma_tilde(X, mu_hat, "mahal", Sigma_true))$vectors[,1],
                  true_vec)
  }
  list(cov = mean(a_cov^2), scm = mean(a_scm^2),
       pd  = mean(a_pd^2),  hsd = mean(a_hsd^2), mah = mean(a_mah^2))
}

# --- Clean data (no outliers) ---
cat("--- Clean data (no outliers) ---\n")
cat(sprintf("%-7s %-5s | %9s | %6s | %7s | %8s | %8s\n",
            "Dist","n","Cov(base)","SCM","Sig~-PD","Sig~-HSD","Sig~-Mah"))
cat(rep("-", 65), "\n", sep = "")

for (dist_name in c("Normal", "t3", "t5", "t10")) {
  for (n in c(50, 100, 200, 500)) {
    gen <- local({
      dn <- dist_name; nn <- n
      function() {
        if (dn == "Normal") return(mvrnorm(nn, rep(0, p_jmva), Sigma_jmva))
        df <- as.integer(sub("t", "", dn))
        gen_t(nn, p_jmva, df, Sigma_jmva)
      }
    })
    m   <- run_mspa(gen, true_PC1, B_jmva, Sigma_jmva)
    fse <- function(x) sprintf("%.3f", m$cov / x)
    cat(sprintf("%-7s n=%-4d | %9s | %6s | %7s | %8s | %8s\n",
                dist_name, n, "1.000",
                fse(m$scm), fse(m$pd), fse(m$hsd), fse(m$mah)))
  }
  cat("\n")
}

# --- Contaminated data: 10% outliers, paper's exact shift ---
cat("--- With 10% contamination  [1000 x max_ij(x_ij), first coord] ---\n")
cat("Paper Section 6.1 exact: n_out = round(0.10*n), shift first column\n\n")
cat(sprintf("%-7s %-5s | %9s | %6s | %7s | %8s\n",
            "Dist","n","Cov(base)","SCM","Sig~-PD","Sig~-HSD"))
cat(rep("-", 52), "\n", sep = "")

for (dist_name in c("Normal", "t3")) {
  for (n in c(100, 200, 500)) {
    gen_c <- local({
      dn <- dist_name; nn <- n
      function() {
        if (dn == "Normal")
          X <- mvrnorm(nn, rep(0, p_jmva), Sigma_jmva)
        else
          X <- gen_t(nn, p_jmva, 3, Sigma_jmva)
        # Paper's exact outlier construction (Section 6.1)
        n_out <- round(0.10 * nn)
        idx   <- sample(nn, n_out)
        shift <- 1000 * max(abs(X))       # 1000 * max_ij(x_ij)
        X[idx, 1] <- X[idx, 1] + shift    # shift first coordinate only
        X
      }
    })
    m   <- run_mspa(gen_c, true_PC1, B_jmva, Sigma_jmva)
    fse <- function(x) sprintf("%.3f", m$cov / x)
    cat(sprintf("%-7s n=%-4d | %9s | %6s | %7s | %8s\n",
                dist_name, n, "1.000",
                fse(m$scm), fse(m$pd), fse(m$hsd)))
  }
  cat("\n")
}

cat("JMVA key claim to verify:\n")
cat("  Clean t3/t5: FSE(Sig~-PD) > FSE(SCM) > 1\n")
cat("  Normal clean: FSE ≈ 1 (small gain)\n")
cat("  Contaminated 10%: Sig~-PD still competitive or better than SCM\n")

# ================================================================
# SECTION 4: TRUE MEAN COMPARISON — Both papers on same data
# ================================================================
cat("\n")
cat("================================================================\n")
cat("SECTION 4  |  TRUE MEAN COMPARISON\n")
cat("           |  JASA Tn  +  JMVA Sigma_tilde  on the same dataset\n")
cat("           |  mu and PC1 are analytically known (no estimation)\n")
cat("================================================================\n")
cat("Setting: p=4, n=100, Sigma=diag(4,3,2,1), B=500\n")
cat("TRUE mu  = (0,0,0,0)   -> JASA H0 is true -> size check\n")
cat("TRUE PC1 = (1,0,0,0)   -> JMVA eigenvector target\n\n")

n_cmp    <- 100
B_cmp    <- 500
Sig_cmp  <- diag(c(4, 3, 2, 1))
true_mu  <- rep(0, 4)               # analytically known: mu = 0
true_pc  <- c(1, 0, 0, 0)          # analytically known: first eigenvector

cat(sprintf("%-8s | %12s | %12s | %12s | %12s\n",
            "Dist",
            "JASA size",
            "FSE(PD/Cov)",
            "FSE(PD/SCM)",
            "Interpretation"))
cat(rep("-", 75), "\n", sep = "")

for (dist_name in c("Normal", "t3", "t5")) {

  rej_jasa <- 0
  a_cov <- a_scm <- a_pd <- numeric(B_cmp)

  for (b in 1:B_cmp) {

    # Generate data at TRUE mu = (0,0,0,0)
    if (dist_name == "Normal") {
      X <- mvrnorm(n_cmp, mu = true_mu, Sigma = Sig_cmp)
    } else {
      df   <- as.integer(sub("t", "", dist_name))
      chi2 <- rchisq(n_cmp, df)
      X    <- mvrnorm(n_cmp, mu = true_mu, Sigma = Sig_cmp) / sqrt(chi2/df)
    }

    # JASA: test H0: mu=0  (TRUE -> reject rate should be ~0.05)
    if (jasa_pval(X) < 0.05) rej_jasa <- rej_jasa + 1

    # JMVA: recover PC1 = (1,0,0,0)  (TRUE -> angle should be small)
    mu_hat   <- colMeans(X)
    a_cov[b] <- pred_angle(eigen(cov(X))$vectors[, 1],                    true_pc)
    a_scm[b] <- pred_angle(eigen(scm(X, mu_hat))$vectors[, 1],            true_pc)
    a_pd[b]  <- pred_angle(eigen(sigma_tilde(X, mu_hat, "pd"))$vectors[,1], true_pc)
  }

  size     <- rej_jasa / B_cmp
  fse_pd_c <- mean(a_cov^2) / mean(a_pd^2)  # PD vs sample cov
  fse_pd_s <- mean(a_scm^2) / mean(a_pd^2)  # PD vs SCM

  size_ok  <- abs(size - 0.05) < 0.025
  fse_ok   <- fse_pd_c > 1

  interp <- sprintf("size[%s] FSE[%s]",
                    ifelse(size_ok, "OK", "!!"),
                    ifelse(fse_ok, "PD wins", "check"))

  cat(sprintf("%-8s | %12.3f | %12.3f | %12.3f | %s\n",
              dist_name, size, fse_pd_c, fse_pd_s, interp))
}

# ================================================================
# FINAL SUMMARY
# ================================================================
cat("\n")
cat("================================================================\n")
cat("RESULTS SUMMARY — what to match against the papers\n")
cat("================================================================\n")
cat("
SECTION 2 (JASA):
  mu0 size  -> all rows ~0.05    (confirms correct null distribution)
  mu1 power -> increases n=20 to n=50  (confirms test has power)
  t3 Table: New column > CQ column throughout
            (paper's core claim: ARE ~ 2.54 under t3 for large p)
  Normal Table: New ≈ CQ  (confirms no efficiency loss under normality)

SECTION 3 (JMVA):
  Clean, t3/t5: FSE(Sig~-PD) > 1  and  FSE(Sig~-PD) > FSE(SCM)
  Clean, Normal: FSE ≈ 1 (small or no gain — expected)
  Contaminated: Sig~-PD maintains FSE advantage over SCM at 10% outliers
  (paper's claim: weighted sign improves efficiency without losing robustness)

SECTION 4 (TRUE MEAN comparison):
  JASA size ~0.05 for all three distributions  [size OK]
  JMVA FSE(PD/Cov) > 1 for t3 and t5          [efficiency gain]
  Ordering: t3 FSE > t5 FSE > Normal FSE       [heavier tail = more gain]
  FSE(PD/SCM) > 1 across all distributions     [weighted beats unweighted]
")
