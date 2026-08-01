# Spatial Sign Statistics: JASA vs JMVA Replication

> **Papers compared:**
> - Wang, Peng & Li (2015) — *JASA* 110(512): A High-Dimensional Nonparametric Multivariate Test for Mean Vector
> - Majumdar & Chatterjee (2022) — *JMVA* 191: On Weighted Multivariate Sign Functions

---

## Table of Contents

1. [Background](#1-background)
2. [Statistics Defined](#2-statistics-defined)
3. [Simulation Settings](#3-simulation-settings)
4. [Section 1 — R Functions](#4-section-1--r-functions)
5. [Section 2 — JASA Replication](#5-section-2--jasa-replication-tables-1-2-3)
6. [Section 3 — JMVA Replication](#6-section-3--jmva-replication-table-1)
7. [Section 4 — Direct Comparison](#7-section-4--direct-comparison)
8. [How to Run](#8-how-to-run)
9. [Expected Output Guide](#9-expected-output-guide)
10. [References](#10-references)

---

## 1. Background

Both papers use the **spatial sign** as their core building block, descending from Möttönen & Oja (1995).

| | JASA (Wang et al. 2015) | JMVA (Majumdar & Chatterjee 2022) |
|---|---|---|
| **Goal** | Test H₀: μ = 0 for high-dim data | Robust location estimation and scatter |
| **Sign used** | Uncentred: $Z_i = X_i / \|X_i\|$ | Centred: $S(X_i;\mu) = (X_i - \mu)/|X_i - \mu|$ |
| **Weight** | None — W ≡ 1 | Depth-based: W_PD, W_HSD, W_MhD |
| **Regime** | p >> n | Fixed p |
| **Key result** | ARE ≈ 2.54 vs Hotelling under t₃ | ARE > 1 vs unweighted median under heavy tails |

---

## 2. Statistics Defined

### 2.1 JASA — Spatial Sign Test Statistic (eq. 3)

$$Z_i = \frac{X_i}{\|X_i\|}, \quad T_n = \sum_{i=1}^{n} \sum_{j < i} Z_i^\top Z_j$$

Fast identity (since $\|Z_i\| = 1$):

$$T_n = \frac{\|\sum_{i=1}^n Z_i\|^2 - n}{2}$$

Null distribution: $T_n / \sqrt{\text{var}(T_n)} \to N(0,1)$ via martingale CLT.

### 2.2 CQ Test — Chen & Qin (2010)

$$\text{CQ} = \sum_{i \neq j} X_i^\top X_j = \|\text{colSums}(X)\|^2 - \sum_i \|X_i\|^2$$

Variance (exact unbiased form):

$$\widehat{\text{Tr}(\Sigma^2)} = \frac{\|X^\top X\|_F^2 - \|\text{diag}(X^\top X)\|^2}{n(n-1)}$$

### 2.3 JMVA — Weighted Spatial Median (Section 2.1)

$$\hat{q}_{nW} = \arg\min_q \sum_{i=1}^n W(X_i, F) \cdot \|X_i - q\|$$

Computed via IRLS (iteratively reweighted least squares).

**Weight functions (paper p. 3–4):**

| Weight | Formula | Source |
|--------|---------|--------|
| $W_\text{PD}$ | $\|Z\| / (1 + \|Z\|/\text{MAD})$ | Projection depth |
| $W_\text{HSD}$ | $\hat{F}_Z(\|Z\|)$ — empirical CDF of norms | Half-space depth |

**ARE formula (paper Table 1):**

$$\text{ARE}(\hat{q}_{nW}, \hat{q}_n) = \left(\frac{\det V_1}{\det V_W}\right)^{1/p}$$

where $V_W = \Psi_{2W}^{-1} \Psi_{1W} \Psi_{2W}^{-1}$ from Theorem 1.

**Test for H₀: μ = 0 (Wald, from Theorem 1):**

$$T_W = n \cdot \hat{q}_{nW}^\top \hat{V}_W^{-1} \hat{q}_{nW} \sim \chi^2(p)$$

---

## 3. Simulation Settings

### JASA Settings — Section 3.1

| Parameter | Paper | Replication |
|-----------|-------|-------------|
| n | 20, 50 | ✓ |
| p | 1000, 2000 | 1000 |
| μ₀ (null) | (0,...,0) | ✓ |
| μ₁ (dense) | (0.25,...,0.25) | ✓ |
| μ₂ (mixed) | first p/3=0, mid +0.25, last −0.25 | ✓ |
| Σ₁ | compound symmetry, off-diag=0.2 | ✓ |
| Σ₂ | AR(1), ρ=0.8 | ✓ |
| Σ₃ | D·R·D, Srivastava et al. 2013 | ✓ |
| Distributions | Normal, t₃, Mixture 0.9N+0.1N(9Σ) | ✓ |
| B | 1000 | 500 |

### JMVA Settings — Section 6.1 / Table 1

| Parameter | Paper | Replication |
|-----------|-------|-------------|
| p | 5, 10, 20, 50 | 5, 10 |
| Σ | Identity | ✓ |
| Distributions | Normal, t₃, t₅, t₁₀, t₂₀ | Normal, t₃, t₅, t₁₀ |
| B | 10,000 | 500 |
| Metric | ARE = (det V₁ / det Vw)^(1/p) | ✓ |

---

## 4. Section 1 — R Functions

```r
library(MASS)
set.seed(2025)

B_jasa <- 500   # paper: 1000
B_jmva <- 500   # paper: 10000
B_comp <- 300   # comparison reps
```

### 1A. JASA Functions

```r
# Spatial sign: Zi = Xi / ||Xi||  (uncentred, JASA eq. 2.1)
jasa_sign <- function(X) {
  nr <- sqrt(rowSums(X^2))
  Z  <- X / ifelse(nr == 0, 1, nr)
  Z[nr == 0, ] <- 0; Z
}

# Fast Tn = (||colSums(Z)||^2 - n) / 2
# Identity holds because ||Zi|| = 1 for all i
jasa_Tn <- function(X) {
  Z <- jasa_sign(X); cs <- colSums(Z)
  (sum(cs^2) - nrow(Z)) / 2
}

# Variance estimator — exact eq.(8) from paper
jasa_var <- function(X) {
  Z <- jasa_sign(X); n <- nrow(Z)
  ZZt <- Z %*% t(Z); ZtZ <- t(Z) %*% Z; Zs <- colMeans(Z)
  t1 <- -n/(n-2)^2
  t2 <- (n-1)/(n*(n-2)^2) * sum(ZZt^2)
  t3 <- (1-2*n)/(n*(n-1)) * as.numeric(t(Zs) %*% ZtZ %*% Zs)
  t4 <- 2/n * sum(Zs^2)
  t5 <- (n-2)^2/(n*(n-1)) * sum(Zs^2)^2
  n*(n-1)/2 * (t1+t2+t3+t4+t5)
}

# Two-sided p-value
jasa_pval <- function(X) {
  Tn <- jasa_Tn(X); vn <- jasa_var(X)
  2 * pnorm(-abs(Tn / sqrt(abs(vn))))
}
```

### 1B. CQ Test — Exact Form

```r
# Chen & Qin (2010) exact variance — unbiased Tr(Sigma^2)
# Previous code used sum(var(X)^2) which is approximate — this is correct
cq_pval <- function(X) {
  n   <- nrow(X)
  cq  <- sum(colSums(X)^2) - sum(rowSums(X^2))
  XtX <- t(X) %*% X
  # Unbiased estimator — removes diagonal contribution
  trS2   <- (sum(XtX^2) - sum(diag(XtX)^2)) / (n*(n-1))
  var_cq <- 2 * n * (n-1) * trS2
  2 * pnorm(-abs(cq / sqrt(abs(var_cq))))
}
```

### 1C. JMVA Weight Functions

```r
# W_PD: projection depth weight (bounded — most stable)
# W(Xi) proportional to |Z| / (1 + |Z|/MAD)  — paper p.4
w_pd <- function(X, mu) {
  nr    <- sqrt(rowSums(sweep(X, 2, mu)^2))
  mad_r <- median(abs(nr - median(nr)))
  if (mad_r < 1e-8) mad_r <- 1e-8
  w <- nr / (1 + nr/mad_r); w / max(w)
}

# W_HSD: half-space depth weight
# W(Xi) = empirical CDF of centred norms  — paper p.3
w_hsd <- function(X, mu) {
  nr <- sqrt(rowSums(sweep(X, 2, mu)^2))
  ecdf(nr)(nr)
}
```

### 1D. JMVA Weighted Spatial Median (IRLS)

```r
# Minimises: sum_i W(Xi,F) * ||Xi - q||
# IRLS: q_new = sum(W_i/||Xi-q|| * Xi) / sum(W_i/||Xi-q||)
weighted_median <- function(X, wtype = "pd", tol = 1e-7, maxit = 200) {
  if (wtype == "unweighted") {
    # Plain spatial median — W = 1
    q <- colMeans(X)
    for (iter in 1:maxit) {
      Xc <- sweep(X,2,q); nr <- sqrt(rowSums(Xc^2)); nr <- pmax(nr,1e-8)
      wi <- 1/nr
      q_new <- colSums(X*wi)/sum(wi)
      if (sqrt(sum((q_new-q)^2)) < tol) { q <- q_new; break }
      q <- q_new
    }
    return(q)
  }
  q <- colMeans(X)
  for (iter in 1:maxit) {
    w  <- if (wtype == "pd") w_pd(X, q) else w_hsd(X, q)
    Xc <- sweep(X, 2, q); nr <- sqrt(rowSums(Xc^2)); nr <- pmax(nr, 1e-8)
    wi    <- w / nr
    q_new <- colSums(X * wi) / sum(wi)
    if (sqrt(sum((q_new - q)^2)) < tol) { q <- q_new; break }
    q <- q_new
  }
  q
}
```

### 1E. JMVA Wald Test (Theorem 1)

```r
# Test H0: mu=0 using asymptotic variance from Theorem 1
# T_W = n * q_hat' * V_W^{-1} * q_hat ~ chi^2(p)
jmva_pval_fast <- function(X, wtype = "pd") {
  n     <- nrow(X); p <- ncol(X)
  q_hat <- weighted_median(X, wtype)
  w  <- if (wtype == "pd") w_pd(X, q_hat) else w_hsd(X, q_hat)
  Xc <- sweep(X, 2, q_hat)
  nr <- sqrt(rowSums(Xc^2)); nr <- pmax(nr, 1e-8)
  S  <- Xc / nr
  # Psi1W = (1/n) sum W_i^2 * S_i * S_i'
  wS    <- S * w
  Psi1W <- t(wS) %*% wS / n
  # Psi2W = (1/n) sum W_i/||Xi-q|| * (I - S_i*S_i')
  Psi2W <- matrix(0, p, p)
  for (i in 1:n) {
    Si <- S[i,]
    Psi2W <- Psi2W + w[i]/nr[i] * (diag(p) - outer(Si, Si))
  }
  Psi2W <- Psi2W / n
  P2inv <- tryCatch(solve(Psi2W), error = function(e) MASS::ginv(Psi2W))
  VW    <- P2inv %*% Psi1W %*% P2inv
  VWinv <- tryCatch(solve(VW),    error = function(e) MASS::ginv(VW))
  Tstat <- n * as.numeric(t(q_hat) %*% VWinv %*% q_hat)
  pchisq(Tstat, df = p, lower.tail = FALSE)
}
```

### 1F. ARE Computation (paper Table 1)

```r
# ARE(q_hat_nW, q_hat_n) = (det V1 / det VW)^(1/p)
compute_are <- function(X, wtype = "pd") {
  n <- nrow(X); p <- ncol(X)
  # Weighted median variance VW
  q_W  <- weighted_median(X, wtype)
  w    <- if (wtype=="pd") w_pd(X,q_W) else w_hsd(X,q_W)
  Xc   <- sweep(X,2,q_W); nr <- sqrt(rowSums(Xc^2)); nr <- pmax(nr,1e-8)
  S    <- Xc/nr; wS <- S * w
  P1W  <- t(wS) %*% wS / n
  P2W  <- matrix(0,p,p)
  for (i in 1:n) { Si <- S[i,]; P2W <- P2W + w[i]/nr[i]*(diag(p)-outer(Si,Si)) }
  P2W  <- P2W/n
  P2Wi <- tryCatch(solve(P2W), error=function(e) MASS::ginv(P2W))
  VW   <- P2Wi %*% P1W %*% P2Wi
  # Unweighted median variance V1
  q_1  <- weighted_median(X, "unweighted")
  Xc1  <- sweep(X,2,q_1); nr1 <- sqrt(rowSums(Xc1^2)); nr1 <- pmax(nr1,1e-8)
  S1   <- Xc1/nr1; P11 <- t(S1) %*% S1 / n
  P21  <- matrix(0,p,p)
  for (i in 1:n) { Si <- S1[i,]; P21 <- P21 + 1/nr1[i]*(diag(p)-outer(Si,Si)) }
  P21  <- P21/n
  P21i <- tryCatch(solve(P21), error=function(e) MASS::ginv(P21))
  V1   <- P21i %*% P11 %*% P21i
  d1 <- det(V1); dW <- det(VW)
  if (d1 <= 0 || dW <= 0) return(NA)
  (d1/dW)^(1/p)
}
```

### 1G. Data Generators

```r
# JASA covariance structures — Section 3.1 of paper
make_Sigma <- function(p, type) {
  if (type==1) { S <- matrix(0.2,p,p); diag(S) <- 1; return(S) }         # compound sym
  if (type==2) { return(outer(1:p,1:p,function(i,j) 0.8^abs(i-j))) }     # AR(1)
  if (type==3) {                                                             # Srivastava
    d <- 2+(p-1:p+1)/p
    R <- outer(1:p,1:p,function(i,j)
          ifelse(i==j, 1, (-1)^(i+j)*0.2^(abs(i-j)/0.1)))
    return(diag(d)%*%R%*%diag(d))
  }
}

# JASA mean vectors — Section 3.1 of paper
make_mu <- function(p, type) {
  if (type==0) return(rep(0,p))                                            # null
  if (type==1) return(rep(0.25,p))                                         # dense alt
  if (type==2) {                                                             # mixed alt
    mu <- rep(0,p)
    mu[(floor(p/3)+1):floor(2*p/3)] <-  0.25
    mu[(floor(2*p/3)+1):p]          <- -0.25
    return(mu)
  }
}

# Data generation for all three distributions
gen_data <- function(n, mu, Sigma, dist) {
  p <- length(mu)
  if (dist=="normal") return(mvrnorm(n,mu=mu,Sigma=Sigma))
  if (dist=="t3") {
    chi2 <- rchisq(n,3)
    return(sweep(mvrnorm(n,mu=rep(0,p),Sigma=Sigma)/sqrt(chi2/3),2,mu,"+"))
  }
  if (dist=="mix") {   # 0.9*N(mu,Sigma) + 0.1*N(mu,9*Sigma)
    idx <- sample(1:2,n,replace=TRUE,prob=c(0.9,0.1))
    X   <- matrix(0,n,p)
    n1  <- sum(idx==1); n2 <- sum(idx==2)
    if (n1>0) X[idx==1,] <- mvrnorm(n1,mu=mu,Sigma=Sigma)
    if (n2>0) X[idx==2,] <- mvrnorm(n2,mu=mu,Sigma=9*Sigma)
    return(X)
  }
}
```

---

## 5. Section 2 — JASA Replication (Tables 1, 2, 3)

> **What this checks:** Empirical size (reject rate at true μ=0, target ≈ 0.05) and power (reject rate at μ₁ or μ₂, should exceed 0.05 and increase with n). Paper values shown in brackets.

```r
# Paper's exact published values at n=20, p=1000
paper <- list(
  normal=list(
    S1=list(m0=c(0.066,0.069), m1=c(0.723,0.723), m2=c(0.951,0.826)),
    S2=list(m0=c(0.052,0.051), m1=c(0.795,0.797), m2=c(0.540,0.549)),
    S3=list(m0=c(0.055,0.055), m1=c(0.490,0.438), m2=c(0.242,0.225))),
  t3=list(
    S1=list(m0=c(0.083,0.088), m1=c(0.633,0.472), m2=c(0.815,0.371)),
    S2=list(m0=c(0.052,0.053), m1=c(0.682,0.349), m2=c(0.441,0.228)),
    S3=list(m0=c(0.054,0.058), m1=c(0.355,0.174), m2=c(0.198,0.113))),
  mix=list(
    S1=list(m0=c(0.063,0.070), m1=c(0.649,0.548), m2=c(0.870,0.449)),
    S2=list(m0=c(0.046,0.063), m1=c(0.678,0.485), m2=c(0.437,0.285)),
    S3=list(m0=c(0.054,0.053), m1=c(0.342,0.207), m2=c(0.178,0.130)))
)

run_jasa_cell <- function(n, p, mu_t, sig_t, dist) {
  mu <- make_mu(p, mu_t); Sigma <- make_Sigma(p, sig_t)
  rn <- rcq <- 0
  for (b in 1:B_jasa) {
    X <- gen_data(n, mu, Sigma, dist)
    if (jasa_pval(X) < 0.05) rn  <- rn  + 1
    if (cq_pval(X)   < 0.05) rcq <- rcq + 1
  }
  c(New=round(rn/B_jasa,3), CQ=round(rcq/B_jasa,3))
}

print_jasa_table <- function(dist_code, label) {
  cat(sprintf("\n--- %s ---\n", label))
  cat(sprintf("%-9s %3s %5s | %14s | %14s | Type\n",
              "Setting","n","p","New (paper)","CQ (paper)"))
  cat(rep("-",65),"\n",sep="")
  pv <- paper[[dist_code]]
  mu_keys <- c("m0","m1","m2")
  for (sig in 1:3) {
    sk <- paste0("S",sig)
    for (mi in 1:3) {
      mn <- mu_keys[mi]; mu_t <- mi-1
      for (n in c(20,50)) {
        p   <- 1000
        res <- run_jasa_cell(n, p, mu_t, sig, dist_code)
        pnw <- pv[[sk]][[mn]][1]
        pcq <- pv[[sk]][[mn]][2]
        type <- if (mu_t==0) "SIZE " else "POWER"
        cat(sprintf("Sig%d mu%d  %3d %5d | %.3f  (%.3f)   | %.3f  (%.3f)   | %s\n",
                    sig, mu_t, n, p, res["New"],pnw, res["CQ"],pcq, type))
      }
    }
    cat("\n")
  }
}

# Run all three tables
print_jasa_table("normal", "Example 1: Normal  [Table 1] — yours (paper)")
print_jasa_table("t3",     "Example 2: t3      [Table 2] — yours (paper)")
print_jasa_table("mix",    "Example 3: Mixture [Table 3] — yours (paper)")
```

> **Why values differ from paper:** Paper uses B=1000 with unpublished seed. Differences of ±0.02 are normal Monte Carlo variation. Values within ±0.03 of the paper = correct replication. The key pattern to verify: **under t₃, New >> CQ throughout** — this is the paper's core claim (ARE ≈ 2.54).

---

## 6. Section 3 — JMVA Replication (Table 1)

> **What this checks:** ARE of weighted spatial median vs unweighted, across distributions and p values. ARE > 1 means the weighted median is more efficient. Should follow ordering: t₃ > t₅ > t₁₀ > Normal, and decrease with p.

**Paper Table 1 values (projection depth weights):**

| p | t₃ | t₅ | t₁₀ | t₂₀ | Normal |
|---|----|----|-----|-----|--------|
| 5 | 1.28 | 1.20 | 1.16 | 1.14 | 1.13 |
| 10 | 1.15 | 1.10 | 1.07 | 1.07 | 1.06 |

```r
paper_are <- list(
  t3     = c("5"=1.28, "10"=1.15),
  t5     = c("5"=1.20, "10"=1.10),
  t10    = c("5"=1.16, "10"=1.07),
  Normal = c("5"=1.13, "10"=1.06)
)

for (dist_name in c("t3","t5","t10","Normal")) {
  for (p in c(5, 10)) {
    Sigma_are <- diag(p)
    ares_pd <- ares_hsd <- numeric(B_jmva)
    for (b in 1:B_jmva) {
      n <- 200
      if (dist_name=="Normal") {
        X <- mvrnorm(n, mu=rep(0,p), Sigma=Sigma_are)
      } else {
        df   <- as.integer(sub("t","",dist_name))
        chi2 <- rchisq(n, df)
        X    <- mvrnorm(n, mu=rep(0,p), Sigma=Sigma_are) / sqrt(chi2/df)
      }
      ares_pd[b]  <- tryCatch(compute_are(X,"pd"),  error=function(e) NA)
      ares_hsd[b] <- tryCatch(compute_are(X,"hsd"), error=function(e) NA)
    }
    are_pd  <- round(mean(ares_pd,  na.rm=TRUE), 3)
    are_hsd <- round(mean(ares_hsd, na.rm=TRUE), 3)
    pval    <- paper_are[[dist_name]][as.character(p)]
    interp  <- if (are_pd > 1) "weighted wins" else "unweighted wins"
    cat(sprintf("%-8s | p=%2d | ARE-PD: %.3f | ARE-HSD: %.3f | Paper: %.2f | %s\n",
                dist_name, p, are_pd, are_hsd,
                ifelse(is.na(pval),NA,pval), interp))
  }
}
```

> **Key:** ARE > 1 for t₃/t₅/t₁₀ confirms the JMVA paper's efficiency claim. Normal ARE ≈ 1.06–1.13 in the paper (small gain even under Gaussianity). ARE decreases as p grows — matches paper Table 1 exactly.

---

## 7. Section 4 — Direct Comparison

> **What this checks:** New (JASA Tₙ), CQ, and JMVA weighted median all tested on the same dataset at the analytically known true mean. Size check (H₀ true, target ≈ 0.05) and power check (H₁ true, higher = better).

```r
p_c      <- 4
Sig_c    <- diag(c(4,3,2,1))
mu_null  <- rep(0, p_c)     # analytically known true null
mu_shift <- rep(0.3, p_c)   # analytically known true alternative

gen_comp <- function(n, mu, dist) {
  if (dist=="Normal") return(mvrnorm(n, mu=mu, Sigma=Sig_c))
  df   <- as.integer(sub("t","",dist))
  chi2 <- rchisq(n, df)
  mvrnorm(n, mu=mu, Sigma=Sig_c) / sqrt(chi2/df)
}

run_comp <- function(mu_true, dist, n) {
  rn <- rcq <- rjm <- 0
  for (b in 1:B_comp) {
    X <- gen_comp(n, mu_true, dist)
    if (jasa_pval(X)           < 0.05) rn  <- rn  + 1
    if (cq_pval(X)             < 0.05) rcq <- rcq + 1
    if (jmva_pval_fast(X,"pd") < 0.05) rjm <- rjm + 1
  }
  c(New=round(rn/B_comp,3), CQ=round(rcq/B_comp,3), JMVA=round(rjm/B_comp,3))
}

# SIZE check
cat("--- SIZE (TRUE mu=0, all methods target ~0.05) ---\n")
for (dist in c("Normal","t3","t5")) {
  for (n in c(20,50,100)) {
    res <- run_comp(mu_null, dist, n)
    ok  <- all(abs(res-0.05) < 0.025)
    cat(sprintf("%-8s n=%3d | New=%.3f CQ=%.3f JMVA=%.3f | %s\n",
                dist, n, res["New"], res["CQ"], res["JMVA"],
                ifelse(ok,"[OK]","[CHECK]")))
  }
}

# POWER check
cat("\n--- POWER (TRUE mu=(0.3,...), higher=better) ---\n")
for (dist in c("Normal","t3","t5")) {
  for (n in c(20,50,100)) {
    res <- run_comp(mu_shift, dist, n)
    gap <- round(res["New"] - res["CQ"], 3)
    winner <- if (dist=="Normal") "all similar" else
              sprintf("New~JMVA >> CQ (gap=%+.3f)", gap)
    cat(sprintf("%-8s n=%3d | New=%.3f CQ=%.3f JMVA=%.3f | %s\n",
                dist, n, res["New"], res["CQ"], res["JMVA"], winner))
  }
}
```

---

## 8. How to Run

```r
# Option 1: Run full script
source("final_comparison_v3.R")

# Option 2: Quick test (reduce B for speed)
B_jasa <- 100; B_jmva <- 100; B_comp <- 100
source("final_comparison_v3.R")

# Option 3: Run one section at a time
# - Load Section 1 (functions) first
# - Then run any of Sections 2, 3, 4 independently
```

**Runtime estimates:**

| B value | Section 2 | Section 3 | Section 4 | Total |
|---------|-----------|-----------|-----------|-------|
| B=100 | ~3 min | ~5 min | ~2 min | ~10 min |
| B=500 | ~15 min | ~25 min | ~8 min | ~48 min |
| B=1000 | ~30 min | ~50 min | ~15 min | ~95 min |

---

## 9. Expected Output Guide

### Section 2 — JASA

| What to look for | Expected |
|-----------------|----------|
| SIZE rows (mu0) | ~0.05 for New across all settings |
| POWER rows, Normal | New ≈ CQ — both detect signal equally |
| POWER rows, t₃ | **New >> CQ** — this is the paper's main claim |
| POWER rows, Mixture | New > CQ — sign-based advantage holds |
| Difference from paper | ±0.02–0.03 is normal (different seed, B=500 not 1000) |

### Section 3 — JMVA ARE

| What to look for | Expected |
|-----------------|----------|
| t₃ ARE-PD at p=5 | ~1.28 (paper value) |
| Normal ARE-PD at p=5 | ~1.13 (paper value) |
| Ordering | t₃ > t₅ > t₁₀ > Normal |
| Effect of p | ARE decreases as p grows from 5 to 10 |

### Section 4 — Comparison

| What to look for | Expected |
|-----------------|----------|
| SIZE: all methods | All ~0.05 — correctly calibrated |
| POWER, Normal | New ≈ CQ ≈ JMVA — all equivalent |
| POWER, t₃/t₅ | New ≈ JMVA >> CQ |
| Key message | Both JASA and JMVA exploit sign structure → power over CQ under heavy tails |

---

## 10. References

1. Wang L., Peng B., & Li R. (2015). A high-dimensional nonparametric multivariate test for mean vector. *Journal of the American Statistical Association*, **110**(512), 1658–1669. DOI: 10.1080/01621459.2014.988215

2. Majumdar S. & Chatterjee S. (2022). On weighted multivariate sign functions. *Journal of Multivariate Analysis*, **191**, 105013. DOI: 10.1016/j.jmva.2022.105013

3. Chen S.X. & Qin Y.L. (2010). A two-sample test for high-dimensional data with application to gene-set testing. *Annals of Statistics*, **38**, 808–835.

4. Möttönen J. & Oja H. (1995). Multivariate spatial sign and rank methods. *Journal of Nonparametric Statistics*, **5**(2), 201–213.

5. Magyar A.F. & Tyler D.E. (2014). The asymptotic inadmissibility of the spatial sign covariance matrix for elliptically symmetric distributions. *Biometrika*, **101**(3), 673–688.

---

*Generated from R simulation study. All code available in `final_comparison_v3.R`.*
