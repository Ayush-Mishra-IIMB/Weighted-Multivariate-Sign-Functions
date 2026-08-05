# ================================================================
#  UNIFIED COMPARISON v2 — EXACT SAME SETTINGS AS JASA PAPER
#  source("final_unified_v2.R")
#
#  ALL THREE METHODS on IDENTICAL settings:
#    p = 1000,  n = 20 and 50
#    Sigma1, Sigma2, Sigma3  (paper Section 3.1)
#    mu0, mu1, mu2           (paper Section 3.1)
#    Distributions: Normal, t3, Mixture
#
#  STATISTICS:
#    New  = JASA Tn  = sum_{i} sum_{j<i} Zi'Zj
#           Zi = Xi / ||Xi||  (uncentred sign)
#
#    CQ   = Chen-Qin 2010 statistic
#
#    JMVA = weighted sign U-statistic
#           SAME structure as Tn but with CENTRED weighted signs
#           Ri = S(Xi; mu_hat) * W(Xi)
#           Tw = sum_{i} sum_{j<i} Ri'Rj
#           null ~ N(0,1) via same standardisation as JASA
#
#  NO Wald test, NO chi-sq approximation, NO matrix inversion
#  NO assumption on p — works at p=1000 same as JASA
# ================================================================

library(MASS)
set.seed(2025)

B <- 300   # paper uses 1000 — increase for closer match

cat("================================================================\n")
cat("  UNIFIED COMPARISON v2 — exact JASA settings\n")
cat("  p=1000 | n=20,50 | Sig1,2,3 | mu0,1,2\n")
cat("  Normal | t3 | Mixture\n")
cat("  B =", B, "\n")
cat("================================================================\n\n")

# ================================================================
# SECTION 1: FUNCTIONS
# ================================================================

# ------ 1A. JASA Tn ------
# Zi = Xi / ||Xi||  (uncentred)
# Tn = sum_{i} sum_{j<i} Zi'Zj = (||colSums(Z)||^2 - n) / 2

jasa_sign <- function(X) {
  nr <- sqrt(rowSums(X^2))
  Z  <- X / ifelse(nr == 0, 1, nr)
  Z[nr == 0, ] <- 0; Z
}

jasa_Tn <- function(X) {
  Z  <- jasa_sign(X)
  cs <- colSums(Z)
  (sum(cs^2) - nrow(Z)) / 2
}

# Variance via exact eq.(8) of paper
jasa_var <- function(X) {
  Z   <- jasa_sign(X); n <- nrow(Z)
  ZZt <- Z %*% t(Z); ZtZ <- t(Z) %*% Z; Zs <- colMeans(Z)
  t1  <- -n/(n-2)^2
  t2  <- (n-1)/(n*(n-2)^2) * sum(ZZt^2)
  t3  <- (1-2*n)/(n*(n-1)) * as.numeric(t(Zs) %*% ZtZ %*% Zs)
  t4  <- 2/n * sum(Zs^2)
  t5  <- (n-2)^2/(n*(n-1)) * sum(Zs^2)^2
  n*(n-1)/2 * (t1+t2+t3+t4+t5)
}

jasa_pval <- function(X) {
  Tn <- jasa_Tn(X); vn <- jasa_var(X)
  2 * pnorm(-abs(Tn / sqrt(abs(vn))))
}

# ------ 1B. CQ exact (Chen & Qin 2010) ------
cq_pval <- function(X) {
  n   <- nrow(X)
  cq  <- sum(colSums(X)^2) - sum(rowSums(X^2))
  XtX <- t(X) %*% X
  trS2   <- (sum(XtX^2) - sum(diag(XtX)^2)) / (n*(n-1))
  var_cq <- 2 * n * (n-1) * trS2
  2 * pnorm(-abs(cq / sqrt(abs(var_cq))))
}

# ------ 1C. JMVA weighted sign U-statistic ------
# SAME U-statistic structure as JASA Tn
# Ri = S(Xi; mu_hat) * Wi  — centred weighted sign
# Tw = sum_{i} sum_{j<i} Ri'Rj
# Standardised: Tw / sqrt(var_Tw) ~ N(0,1)
# Works at ANY p — no matrix inversion needed

# Centred sign: S(Xi; mu) = (Xi - mu) / ||Xi - mu||
jmva_sign <- function(X, mu) {
  Xc <- sweep(X, 2, mu)
  nr <- sqrt(rowSums(Xc^2))
  S  <- Xc / ifelse(nr == 0, 1, nr)
  S[nr == 0, ] <- 0; S
}

# W_PD: projection depth weight — bounded, stable at any p
w_pd <- function(X, mu) {
  nr    <- sqrt(rowSums(sweep(X, 2, mu)^2))
  mad_r <- median(abs(nr - median(nr)))
  if (mad_r < 1e-8) mad_r <- 1e-8
  w <- nr / (1 + nr/mad_r)
  w / max(w)
}

# W_HSD: half-space depth weight
w_hsd <- function(X, mu) {
  nr <- sqrt(rowSums(sweep(X, 2, mu)^2))
  ecdf(nr)(nr)
}

# JMVA test statistic and p-value
# Key: same fast formula as JASA but Ri have different norms
# Fast formula: Tw = (||colSums(R)||^2 - sum(rowSums(R^2))) / 2
# Variance:     var_Tw = n(n-1)/2 * Tr(Bw^2)
#               Bw = (1/n) R'R  (the weighted analogue of B in JASA)
jmva_pval <- function(X, wtype = "pd") {
  n      <- nrow(X)
  mu_hat <- colMeans(X)                      # estimate centre from data

  # Weighted centred signs
  S <- jmva_sign(X, mu_hat)
  w <- if (wtype == "pd") w_pd(X, mu_hat) else w_hsd(X, mu_hat)
  R <- S * w                                 # Ri = S(Xi; mu_hat) * Wi

  # Test statistic — fast formula
  cs     <- colSums(R)                       # p-vector
  row_sq <- sum(rowSums(R^2))               # sum of ||Ri||^2
  Tw     <- (sum(cs^2) - row_sq) / 2        # U-statistic

  # Variance — Bw = (1/n) R'R, Tr(Bw^2) = ||Bw||_F^2
  Bw    <- t(R) %*% R / n                   # p x p  (same as JASA B but weighted)
  TrBw2 <- sum(Bw^2)                        # Tr(Bw^2)
  var_Tw <- n*(n-1)/2 * TrBw2

  # Standardised: N(0,1) null
  2 * pnorm(-abs(Tw / sqrt(abs(var_Tw))))
}

# ------ 1D. Data generators (exact JASA Section 3.1) ------

make_Sigma <- function(p, type) {
  if (type == 1) {
    # Sigma1: compound symmetry, off-diag = 0.2
    S <- matrix(0.2, p, p); diag(S) <- 1; return(S)
  }
  if (type == 2) {
    # Sigma2: AR(1), rho = 0.8
    return(outer(1:p, 1:p, function(i, j) 0.8^abs(i-j)))
  }
  if (type == 3) {
    # Sigma3: Srivastava, Katayama & Kano (2013)
    d <- 2 + (p - 1:p + 1)/p
    R <- outer(1:p, 1:p, function(i, j)
          ifelse(i==j, 1, (-1)^(i+j) * 0.2^(abs(i-j)/0.1)))
    return(diag(d) %*% R %*% diag(d))
  }
}

make_mu <- function(p, type) {
  if (type == 0) return(rep(0, p))            # mu0: null
  if (type == 1) return(rep(0.25, p))         # mu1: dense alternative
  if (type == 2) {                             # mu2: mixed alternative
    mu <- rep(0, p)
    mu[(floor(p/3)+1):floor(2*p/3)] <-  0.25
    mu[(floor(2*p/3)+1):p]          <- -0.25
    return(mu)
  }
}

gen_data <- function(n, mu, Sigma, dist) {
  p <- length(mu)
  if (dist == "normal")
    return(mvrnorm(n, mu = mu, Sigma = Sigma))
  if (dist == "t3") {
    chi2 <- rchisq(n, 3)
    Z    <- mvrnorm(n, mu = rep(0,p), Sigma = Sigma) / sqrt(chi2/3)
    return(sweep(Z, 2, mu, "+"))
  }
  if (dist == "mix") {
    idx <- sample(1:2, n, replace = TRUE, prob = c(0.9, 0.1))
    X   <- matrix(0, n, p)
    n1  <- sum(idx==1); n2 <- sum(idx==2)
    if (n1>0) X[idx==1,] <- mvrnorm(n1, mu=mu, Sigma=Sigma)
    if (n2>0) X[idx==2,] <- mvrnorm(n2, mu=mu, Sigma=9*Sigma)
    return(X)
  }
}

cat("All functions loaded.\n\n")

# ================================================================
# SECTION 2: SIMULATION CELL
# ================================================================

run_cell <- function(n, p, mu_type, sigma_type, dist) {
  mu    <- make_mu(p, mu_type)       # TRUE mean — analytically known
  Sigma <- make_Sigma(p, sigma_type)
  rej_new <- rej_cq <- rej_pd <- rej_hsd <- 0
  for (b in 1:B) {
    X <- gen_data(n, mu, Sigma, dist)
    if (jasa_pval(X)          < 0.05) rej_new  <- rej_new  + 1
    if (cq_pval(X)            < 0.05) rej_cq   <- rej_cq   + 1
    if (jmva_pval(X, "pd")   < 0.05) rej_pd   <- rej_pd   + 1
    if (jmva_pval(X, "hsd")  < 0.05) rej_hsd  <- rej_hsd  + 1
  }
  c(New = round(rej_new/B, 3),
    CQ  = round(rej_cq /B, 3),
    PD  = round(rej_pd /B, 3),
    HSD = round(rej_hsd/B, 3))
}

# ================================================================
# SECTION 3: PRINT TABLE — exact JASA structure
# ================================================================

# Published paper values at n=20, p=1000 for New and CQ
paper_new <- list(
  normal=list(S1=c(0.066,0.723,0.951),S2=c(0.052,0.795,0.540),S3=c(0.055,0.490,0.242)),
  t3    =list(S1=c(0.083,0.633,0.815),S2=c(0.052,0.682,0.441),S3=c(0.054,0.355,0.198)),
  mix   =list(S1=c(0.063,0.649,0.870),S2=c(0.046,0.678,0.437),S3=c(0.054,0.342,0.178))
)
paper_cq <- list(
  normal=list(S1=c(0.069,0.723,0.826),S2=c(0.051,0.797,0.549),S3=c(0.055,0.438,0.225)),
  t3    =list(S1=c(0.088,0.472,0.371),S2=c(0.053,0.349,0.228),S3=c(0.058,0.174,0.113)),
  mix   =list(S1=c(0.070,0.548,0.449),S2=c(0.063,0.485,0.285),S3=c(0.053,0.207,0.130))
)

print_table <- function(dist_code, dist_label) {

  cat("\n", paste(rep("=",82),collapse=""), "\n", sep="")
  cat(" ", dist_label, "\n")
  cat(" New(paper) CQ(paper) = published JASA values at n=20, p=1000\n")
  cat(" JMVA-PD / JMVA-HSD  = weighted sign U-statistic, same structure as Tn\n")
  cat(paste(rep("=",82),collapse=""), "\n\n", sep="")

  cat(sprintf("%-10s | %3s | %4s | %-14s | %-14s | %-9s | %-9s | %s\n",
              "Setting","n","p",
              "New (paper)",
              "CQ  (paper)",
              "JMVA-PD","JMVA-HSD",
              "Type"))
  cat(paste(rep("-",90),collapse=""),"\n",sep="")

  pn <- paper_new[[dist_code]]
  pc <- paper_cq[[dist_code]]

  for (sig in 1:3) {
    sk <- paste0("S", sig)

    for (mu_t in 0:2) {

      for (n in c(20, 50)) {
        p   <- 1000
        res <- run_cell(n, p, mu_t, sig, dist_code)

        # Paper values only at n=20
        if (n == 20) {
          nw_str <- sprintf("%.3f (%.3f)", res["New"], pn[[sk]][mu_t+1])
          cq_str <- sprintf("%.3f (%.3f)", res["CQ"],  pc[[sk]][mu_t+1])
        } else {
          nw_str <- sprintf("%.3f        ", res["New"])
          cq_str <- sprintf("%.3f        ", res["CQ"])
        }

        # Type and interpretation
        if (mu_t == 0) {
          ok   <- abs(res["New"] - 0.05) < 0.025
          type <- sprintf("SIZE  [%s]", ifelse(ok,"OK","CHECK"))
        } else {
          gap  <- res["New"] - res["CQ"]
          wins <- res["PD"]  > res["CQ"] + 0.03
          type <- sprintf("POWER | New-CQ=%+.3f | JMVA-PD%s CQ",
                          gap, ifelse(wins," >>"," ~"))
        }

        cat(sprintf("Sig%d mu%d    | %3d | %4d | %-14s | %-14s | %9.3f | %9.3f | %s\n",
                    sig, mu_t, n, p,
                    nw_str, cq_str,
                    res["PD"], res["HSD"],
                    type))
      }
    }
    cat("\n")
  }
}

# ================================================================
# SECTION 4: RUN ALL THREE TABLES
# ================================================================

print_table("normal",
  "EXAMPLE 1: Multivariate NORMAL — replicates Table 1")

print_table("t3",
  "EXAMPLE 2: Multivariate t3 (heavy tails) — replicates Table 2")

print_table("mix",
  "EXAMPLE 3: Scale Mixture 0.9N+0.1N(9S) — replicates Table 3")

# ================================================================
# SECTION 5: SUMMARY
# ================================================================

cat("\n", paste(rep("=",82),collapse=""), "\n", sep="")
cat("  SUMMARY — What to tell the professor\n")
cat(paste(rep("=",82),collapse=""), "\n\n", sep="")

cat("JMVA STATISTIC USED:
  R_i = S(X_i; mu_hat) * W_i    (centred weighted sign)
  T_W = sum_{i} sum_{j<i} R_i'R_j    (same U-statistic as JASA Tn)
  Standardised: T_W / sqrt(var_W) ~ N(0,1)
  var_W = n(n-1)/2 * Tr(Bw^2),  Bw = (1/n) R'R
  Works at p=1000 — no matrix inversion, no approximation

DIFFERENCE FROM JASA Tn:
  JASA: Zi = Xi / ||Xi||          (uncentred, ||Zi|| = 1 always)
  JMVA: Ri = S(Xi;mu) * Wi        (centred, ||Ri|| = Wi in (0,1])
  Same U-statistic pair structure — directly comparable

EXPECTED PATTERNS:
  SIZE (mu0):  all four methods ~0.05 across Sig1,2,3
  POWER Normal: New ~ CQ ~ JMVA-PD ~ JMVA-HSD
  POWER t3:     New >> CQ  and  JMVA-PD >> CQ  (sign methods win)
  POWER Mix:    New > CQ   and  JMVA-PD > CQ
  Ordering across Sigma: same pattern in Sig1, Sig2, Sig3
")
cat(paste(rep("=",82),collapse=""),"\n",sep="")
