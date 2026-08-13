#  SOURCE: JMVA paper Section 2, end of page 3
#    "W(X, F) = W(Z, FZ) where Z = Sigma^{-1/2}(X - mu)"
# ================================================================

library(MASS)
set.seed(2025)

B <- 300   # paper uses 1000

cat("================================================================\n")
cat("  UNIFIED COMPARISON v4\n")
cat("  PROPER FIX: full Sigma^{-1/2} via eigendecomposition\n")
cat("  Z = Sigma^{-1/2}(Xi - mu) as per JMVA paper Section 2\n")
cat("  p=1000 | n=20,50 | Sig1,2,3 | mu0,1,2 | B =", B, "\n")
cat("================================================================\n\n")

# ================================================================
# SECTION 1: ALL FUNCTIONS
# ================================================================

# ------ 1A. JASA Tn (unchanged) ------

jasa_sign <- function(X) {
  nr <- sqrt(rowSums(X^2))
  Z  <- X / ifelse(nr == 0, 1, nr)
  Z[nr == 0, ] <- 0; Z
}

jasa_Tn <- function(X) {
  Z  <- jasa_sign(X); cs <- colSums(Z)
  (sum(cs^2) - nrow(Z)) / 2
}

jasa_var <- function(X) {
  Z <- jasa_sign(X); n <- nrow(Z)
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

# ------ 1B. CQ exact (unchanged) ------

cq_pval <- function(X) {
  n   <- nrow(X)
  cq  <- sum(colSums(X)^2) - sum(rowSums(X^2))
  XtX <- t(X) %*% X
  trS2   <- (sum(XtX^2) - sum(diag(XtX)^2)) / (n*(n-1))
  var_cq <- 2 * n * (n-1) * trS2
  2 * pnorm(-abs(cq / sqrt(abs(var_cq))))
}

# ------ 1C. Precompute Sigma^{-1/2} via eigendecomposition ------
#
# WHY eigendecomposition and NOT chol(solve(Sigma)):
#   chol(solve(Sigma)) fails for near-singular Sigma (e.g. high rho AR)
#   Eigendecomposition is numerically stable for any Sigma
#
# Sigma = V D V'  (eigendecomposition)
# Sigma^{-1/2} = V D^{-1/2} V'
#
# Computed ONCE per (Sigma, mu) combination — not inside B loop
# This makes it affordable even at p=1000

compute_Sinvhalf <- function(Sigma) {
  # Eigendecomposition: Sigma = V D V'
  eig  <- eigen(Sigma, symmetric = TRUE)
  V    <- eig$vectors          # p x p eigenvector matrix
  d    <- eig$values           # p eigenvalues
  # Clamp tiny/negative eigenvalues for numerical safety
  d    <- pmax(d, 1e-10)
  # Sigma^{-1/2} = V D^{-1/2} V'
  V %*% diag(1/sqrt(d)) %*% t(V)
}

# Standardise: Z = Sigma^{-1/2}(X - mu)
# Sinvhalf is precomputed outside the B loop
standardise <- function(X, mu, Sinvhalf) {
  Xc <- sweep(X, 2, mu)        # n x p: Xi - mu
  Xc %*% Sinvhalf              # n x p: Z = (Xi - mu) Sigma^{-1/2}
  # Note: Sinvhalf is symmetric so (Xi-mu) Sinvhalf = Sinvhalf (Xi-mu)
}

# ------ 1D. JMVA weight functions — correctly use Z = Sigma^{-1/2}(Xi-mu) ------

# W_PD: projection depth weight
# W(Xi) ∝ ||Zi|| / (1 + ||Zi|| / MAD(||Z||))
# where Zi = Sigma^{-1/2}(Xi - mu)   [JMVA paper p.4]
w_pd <- function(Z_mat) {
  # Z_mat is already standardised: n x p matrix of Zi vectors
  nr    <- sqrt(rowSums(Z_mat^2))     # ||Zi||
  mad_r <- median(abs(nr - median(nr)))
  if (mad_r < 1e-8) mad_r <- 1e-8
  w <- nr / (1 + nr / mad_r)
  w / max(w)
}

# W_HSD: half-space depth weight
# W(Xi) = empirical CDF of ||Zi||
# where Zi = Sigma^{-1/2}(Xi - mu)   [JMVA paper p.3]
w_hsd <- function(Z_mat) {
  nr <- sqrt(rowSums(Z_mat^2))        # ||Zi||
  ecdf(nr)(nr)
}

# ------ 1E. JMVA centred sign (in original space) ------

jmva_sign <- function(X, mu) {
  Xc <- sweep(X, 2, mu)
  nr <- sqrt(rowSums(Xc^2))
  S  <- Xc / ifelse(nr == 0, 1, nr)
  S[nr == 0, ] <- 0; S
}

# ------ 1F. JMVA weighted sign U-statistic ------
#
# KEY POINT:
#   Sign S(Xi; mu_hat) is computed in ORIGINAL space (not standardised)
#   Weight W(Xi) is computed from STANDARDISED Z = Sigma^{-1/2}(Xi - mu_hat)
#   This matches JMVA paper: R(Xi; mu, F) = S(Xi, mu) * W(Xi, F)
#   where W depends on the standardised version as per Section 2
#
# Sinvhalf: precomputed Sigma^{-1/2} passed in from run_cell()

jmva_pval <- function(X, wtype = "pd", Sinvhalf = NULL) {
  n      <- nrow(X); p <- ncol(X)
  mu_hat <- colMeans(X)

  # If no Sinvhalf provided use identity (fallback)
  if (is.null(Sinvhalf)) Sinvhalf <- diag(p)

  # Step 1: centred sign in original space
  S <- jmva_sign(X, mu_hat)

  # Step 2: standardised Z for weight computation
  Z_mat <- standardise(X, mu_hat, Sinvhalf)   # Z = Sigma^{-1/2}(Xi - mu_hat)

  # Step 3: weights based on ||Zi|| (standardised norms)
  w <- switch(wtype,
    "pd"  = w_pd(Z_mat),
    "hsd" = w_hsd(Z_mat)
  )

  # Step 4: weighted sign vectors Ri = S(Xi; mu_hat) * Wi
  R <- S * w

  # Step 5: U-statistic (fast formula)
  cs     <- colSums(R)
  row_sq <- sum(rowSums(R^2))
  Tw     <- (sum(cs^2) - row_sq) / 2

  # Step 6: variance and p-value
  Bw     <- t(R) %*% R / n
  var_Tw <- n*(n-1)/2 * sum(Bw^2)
  2 * pnorm(-abs(Tw / sqrt(abs(var_Tw))))
}

# ------ 1G. Data generators (exact JASA Section 3.1) ------

make_Sigma <- function(p, type) {
  if (type == 1) {
    S <- matrix(0.2, p, p); diag(S) <- 1; return(S)
  }
  if (type == 2) {
    return(outer(1:p, 1:p, function(i,j) 0.8^abs(i-j)))
  }
  if (type == 3) {
    d <- 2 + (p - 1:p + 1)/p
    R <- outer(1:p, 1:p, function(i,j)
          ifelse(i==j, 1, (-1)^(i+j) * 0.2^(abs(i-j)/0.1)))
    return(diag(d) %*% R %*% diag(d))
  }
}

make_mu <- function(p, type) {
  if (type == 0) return(rep(0, p))
  if (type == 1) return(rep(0.25, p))
  if (type == 2) {
    mu <- rep(0, p)
    mu[(floor(p/3)+1):floor(2*p/3)] <-  0.25
    mu[(floor(2*p/3)+1):p]          <- -0.25
    return(mu)
  }
}

gen_data <- function(n, mu, Sigma, dist) {
  p <- length(mu)
  if (dist == "normal")
    return(mvrnorm(n, mu=mu, Sigma=Sigma))
  if (dist == "t3") {
    chi2 <- rchisq(n, 3)
    Z    <- mvrnorm(n, mu=rep(0,p), Sigma=Sigma) / sqrt(chi2/3)
    return(sweep(Z, 2, mu, "+"))
  }
  if (dist == "mix") {
    idx <- sample(1:2, n, replace=TRUE, prob=c(0.9, 0.1))
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
# Sigma^{-1/2} computed ONCE per cell via eigendecomposition
# Then passed into jmva_pval() for all B replications

run_cell <- function(n, p, mu_type, sigma_type, dist) {

  mu    <- make_mu(p, mu_type)
  Sigma <- make_Sigma(p, sigma_type)

  # Precompute Sigma^{-1/2} once — not inside B loop
  cat(sprintf("    Computing Sigma^{-1/2} for Sig%d...", sigma_type))
  Sinvhalf <- compute_Sinvhalf(Sigma)
  cat(" done\n")

  rej_new <- rej_cq <- rej_pd <- rej_hsd <- 0

  for (b in 1:B) {
    X <- gen_data(n, mu, Sigma, dist)

    if (jasa_pval(X)                    < 0.05) rej_new  <- rej_new  + 1
    if (cq_pval(X)                      < 0.05) rej_cq   <- rej_cq   + 1
    if (jmva_pval(X, "pd",  Sinvhalf)  < 0.05) rej_pd   <- rej_pd   + 1
    if (jmva_pval(X, "hsd", Sinvhalf)  < 0.05) rej_hsd  <- rej_hsd  + 1
  }

  c(New = round(rej_new/B, 3),
    CQ  = round(rej_cq /B, 3),
    PD  = round(rej_pd /B, 3),
    HSD = round(rej_hsd/B, 3))
}

# ================================================================
# SECTION 3: PRINT TABLE
# ================================================================

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

  cat("\n", paste(rep("=",86),collapse=""), "\n", sep="")
  cat("  ", dist_label, "\n")
  cat("  New(paper) CQ(paper) = published JASA values at n=20, p=1000\n")
  cat("  JMVA-PD / JMVA-HSD  = Z=Sigma^{-1/2}(Xi-mu) via eigendecomposition\n")
  cat(paste(rep("=",86),collapse=""), "\n\n", sep="")

  cat(sprintf("%-10s | %3s | %4s | %-14s | %-14s | %-9s | %-9s | %s\n",
              "Setting","n","p",
              "New (paper)","CQ  (paper)",
              "JMVA-PD","JMVA-HSD","Type"))
  cat(paste(rep("-",94),collapse=""),"\n",sep="")

  pn <- paper_new[[dist_code]]
  pc <- paper_cq[[dist_code]]

  for (sig in 1:3) {
    sk <- paste0("S",sig)
    for (mu_t in 0:2) {
      for (n in c(20, 50)) {
        p   <- 1000
        cat(sprintf("  Running: Sig%d mu%d n=%d dist=%s\n",
                    sig, mu_t, n, dist_code))
        res <- run_cell(n, p, mu_t, sig, dist_code)

        if (n == 20) {
          nw_str <- sprintf("%.3f (%.3f)", res["New"], pn[[sk]][mu_t+1])
          cq_str <- sprintf("%.3f (%.3f)", res["CQ"],  pc[[sk]][mu_t+1])
        } else {
          nw_str <- sprintf("%.3f        ", res["New"])
          cq_str <- sprintf("%.3f        ", res["CQ"])
        }

        if (mu_t == 0) {
          ok   <- abs(res["New"] - 0.05) < 0.025
          ok_pd <- abs(res["PD"] - 0.05) < 0.025
          type <- sprintf("SIZE | New[%s] PD[%s]",
                          ifelse(ok,"OK","CHECK"),
                          ifelse(ok_pd,"OK","CHECK"))
        } else {
          gap   <- res["New"] - res["CQ"]
          wins  <- res["PD"] > res["CQ"] + 0.03
          type  <- sprintf("POWER | New-CQ=%+.3f | JMVA-PD%s CQ",
                           gap, ifelse(wins," >>"," ~"))
        }

        cat(sprintf("Sig%d mu%d    | %3d | %4d | %-14s | %-14s | %9.3f | %9.3f | %s\n",
                    sig, mu_t, n, p,
                    nw_str, cq_str,
                    res["PD"], res["HSD"], type))
      }
    }
    cat("\n")
  }
}

# ================================================================
# SECTION 4: RUN ALL THREE TABLES
# ================================================================

print_table("normal", "EXAMPLE 1: Normal  — Table 1")
print_table("t3",     "EXAMPLE 2: t3      — Table 2")
print_table("mix",    "EXAMPLE 3: Mixture — Table 3")

# ================================================================
# SECTION 5: SUMMARY
# ================================================================

cat("\n", paste(rep("=",86),collapse=""), "\n", sep="")
cat("  SUMMARY — What changed and why\n")
cat(paste(rep("=",86),collapse=""), "\n\n", sep="")
cat("
WHAT CHANGED FROM v3:
  v3 used diag(diag(Sigma)) as approximation of Sigma^{-1/2}
  For Sigma1 and Sigma2: all diagonal entries = 1
  So diag approx = Identity → Z = Xi - mu → NO standardisation
  This is why v3 showed same results as before for Sigma1, Sigma2

  v4 uses FULL Sigma^{-1/2} via eigendecomposition:
  Sigma = V D V'  →  Sigma^{-1/2} = V D^{-1/2} V'
  Computed once per cell, not inside B loop

HOW WEIGHTS NOW WORK:
  Z_i  = Sigma^{-1/2}(Xi - mu_hat)  [standardised — for weights]
  S_i  = (Xi - mu_hat)/||Xi - mu_hat||  [centred sign — original space]
  W_i  = f(||Z_i||)  [depth weight from standardised norm]
  R_i  = S_i * W_i   [weighted sign]
  T_W  = sum_{i} sum_{j<i} Ri'Rj  [U-statistic]

EXPECTED PATTERNS NOW:
  SIZE (mu0):  JMVA-PD ~ 0.05  [properly calibrated]
  POWER Normal: JMVA-PD ~ New  [no advantage, similar to JASA]
  POWER t3:     JMVA-PD > CQ   [sign methods gain over CQ]
  POWER Mix:    JMVA-PD > CQ   [heavy tail advantage]
  JMVA-HSD:    stable values   [ecdf of ||Zi|| well-behaved]
")
cat(paste(rep("=",86),collapse=""),"\n",sep="")
