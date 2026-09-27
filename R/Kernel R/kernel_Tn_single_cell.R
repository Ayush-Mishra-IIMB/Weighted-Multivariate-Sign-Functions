# ================================================================
#  KERNEL Tn PERMUTATION TEST — single-cell correctness check
#  Extends JASA_vs_JMVA_statistic_comparison.R with one new method:
#    T_n^K = 1/C(n,2) * sum_{i<j} K(y_i, y_j),   Y = A X
#    A: n x n, a_ij ~ iid N(0,1), symmetrized + shifted to be PD
#    K: RBF kernel, median-heuristic bandwidth (fixed)
#    p-value via PERMUTATION (not closed-form asymptotics)
#
#  SCOPE FOR THIS RUN (deliberately small, to check correctness first):
#    Sigma: Sig1 only | dist: normal only | mu: mu0 (size), mu1 (power)
#    n = 20, p = 1000 | B = 100 outer reps | B_perm = 100 permutations
#  Once this checks out, scale B/B_perm/Sigma/mu/dist up to match the
#  full grid in the original script.
# ================================================================

library(MASS)
set.seed(2025)

B      <- 100   # outer simulation replications (paper/original script uses 300)
B_perm <- 100   # permutations per replication for the kernel test's p-value

cat("================================================================\n")
cat("  KERNEL Tn — single-cell check (Sig1, normal, mu0 & mu1, n=20,p=1000)\n")
cat("  B =", B, " | B_perm =", B_perm, "\n")
cat("================================================================\n\n")

# ================================================================
# SECTION 1: FUNCTIONS COPIED FROM THE ORIGINAL SCRIPT (unchanged)
# ================================================================

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

cq_pval <- function(X) {
  n   <- nrow(X)
  cq  <- sum(colSums(X)^2) - sum(rowSums(X^2))
  XtX <- t(X) %*% X
  trS2   <- (sum(XtX^2) - sum(diag(XtX)^2)) / (n*(n-1))
  var_cq <- 2 * n * (n-1) * trS2
  2 * pnorm(-abs(cq / sqrt(abs(var_cq))))
}

compute_Sinvhalf <- function(Sigma) {
  eig  <- eigen(Sigma, symmetric = TRUE)
  V    <- eig$vectors
  d    <- eig$values
  d    <- pmax(d, 1e-10)
  V %*% diag(1/sqrt(d)) %*% t(V)
}

standardise <- function(X, mu, Sinvhalf) {
  Xc <- sweep(X, 2, mu)
  Xc %*% Sinvhalf
}

w_pd <- function(Z_mat) {
  nr    <- sqrt(rowSums(Z_mat^2))
  mad_r <- median(abs(nr - median(nr)))
  if (mad_r < 1e-8) mad_r <- 1e-8
  w <- nr / (1 + nr / mad_r)
  w / max(w)
}

w_hsd <- function(Z_mat) {
  nr <- sqrt(rowSums(Z_mat^2))
  ecdf(nr)(nr)
}

jmva_sign <- function(X, mu) {
  Xc <- sweep(X, 2, mu)
  nr <- sqrt(rowSums(Xc^2))
  S  <- Xc / ifelse(nr == 0, 1, nr)
  S[nr == 0, ] <- 0; S
}

jmva_pval <- function(X, wtype = "pd", Sinvhalf = NULL) {
  n      <- nrow(X); p <- ncol(X)
  mu_hat <- colMeans(X)
  if (is.null(Sinvhalf)) Sinvhalf <- diag(p)
  S <- jmva_sign(X, mu_hat)
  Z_mat <- standardise(X, mu_hat, Sinvhalf)
  w <- switch(wtype, "pd" = w_pd(Z_mat), "hsd" = w_hsd(Z_mat))
  R <- S * w
  cs     <- colSums(R)
  row_sq <- sum(rowSums(R^2))
  Tw     <- (sum(cs^2) - row_sq) / 2
  Bw     <- t(R) %*% R / n
  var_Tw <- n*(n-1)/2 * sum(Bw^2)
  2 * pnorm(-abs(Tw / sqrt(abs(var_Tw))))
}

make_Sigma <- function(p, type) {
  if (type == 1) { S <- matrix(0.2, p, p); diag(S) <- 1; return(S) }
  if (type == 2) return(outer(1:p, 1:p, function(i,j) 0.8^abs(i-j)))
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
  if (dist == "normal") return(mvrnorm(n, mu=mu, Sigma=Sigma))
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

cat("Copied JASA/CQ/JMVA functions loaded.\n\n")

# ================================================================
# SECTION 2: NEW — KERNEL Tn AND PERMUTATION P-VALUE
# ================================================================

## ---- A: n x n, iid N(0,1), symmetrized + eigenvalue-shifted to be PD ----
## ASSUMPTION (flagged): PD achieved by symmetrizing + shifting eigenvalues.
## Confirm with advisor if a different PD construction is intended.
generate_A <- function(n) {
  M <- matrix(rnorm(n * n), n, n)
  S <- (M + t(M)) / 2
  eig <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  shift <- max(0, -min(eig)) + 0.01
  S + diag(shift, n)
}

## ---- Fixed RBF kernel, median-heuristic bandwidth ----
rbf_kernel_matrix <- function(Y) {
  D2 <- as.matrix(dist(Y))^2
  sigma2 <- median(D2[upper.tri(D2)])
  if (sigma2 == 0) sigma2 <- 1e-6
  exp(-D2 / (2 * sigma2))
}

## ---- Tn^K = 1/C(n,2) * sum_{i<j} K(yi,yj),  Y = A X (matrix product) ----
kernel_Tn <- function(X, A) {
  n <- nrow(X)
  Y <- A %*% X                    # n x p, matrix product (mixes across samples)
  K <- rbf_kernel_matrix(Y)
  sum(K[upper.tri(K)]) / choose(n, 2)
}

## ---- Permutation p-value: permute ROWS of X, A stays FIXED ----
kernel_perm_pval <- function(X, A, B_perm) {
  n <- nrow(X)
  Tn_obs <- kernel_Tn(X, A)

  Tn_null <- numeric(B_perm)
  for (b in 1:B_perm) {
    X_perm <- X[sample(1:n), , drop = FALSE]
    Tn_null[b] <- kernel_Tn(X_perm, A)
  }

  (1 + sum(Tn_null >= Tn_obs)) / (B_perm + 1)
}

cat("Kernel Tn + permutation p-value functions loaded.\n\n")

# ================================================================
# SECTION 3: SINGLE-CELL SIMULATION (Sig1, normal, mu0 & mu1, n=20,p=1000)
# ================================================================

run_cell_kernel <- function(n, p, mu_type, sigma_type, dist, B, B_perm) {

  mu    <- make_mu(p, mu_type)
  Sigma <- make_Sigma(p, sigma_type)

  cat(sprintf("    Computing Sigma^{-1/2} for Sig%d...", sigma_type))
  Sinvhalf <- compute_Sinvhalf(Sigma)
  cat(" done\n")

  rej_new <- rej_cq <- rej_pd <- rej_hsd <- rej_kernel <- 0

  for (b in 1:B) {
    X <- gen_data(n, mu, Sigma, dist)
    A <- generate_A(n)   # NOTE: fresh A drawn each replication (flagged assumption)

    if (jasa_pval(X)                  < 0.05) rej_new    <- rej_new    + 1
    if (cq_pval(X)                    < 0.05) rej_cq     <- rej_cq     + 1
    if (jmva_pval(X, "pd",  Sinvhalf) < 0.05) rej_pd     <- rej_pd     + 1
    if (jmva_pval(X, "hsd", Sinvhalf) < 0.05) rej_hsd    <- rej_hsd    + 1
    if (kernel_perm_pval(X, A, B_perm) < 0.05) rej_kernel <- rej_kernel + 1

    if (b %% 20 == 0) cat("      rep", b, "/", B, "\n")
  }

  c(New = round(rej_new/B, 3),
    CQ  = round(rej_cq /B, 3),
    PD  = round(rej_pd /B, 3),
    HSD = round(rej_hsd/B, 3),
    Kernel = round(rej_kernel/B, 3))
}

n <- 20; p <- 1000

cat("\n---- SIZE check: Sig1, normal, mu0 (null) ----\n")
res_size <- run_cell_kernel(n, p, mu_type = 0, sigma_type = 1, dist = "normal",
                             B = B, B_perm = B_perm)
cat("\nResults (should all be near 0.05 under H0):\n")
print(res_size)

cat("\n---- POWER check: Sig1, normal, mu1 (alternative) ----\n")
res_power <- run_cell_kernel(n, p, mu_type = 1, sigma_type = 1, dist = "normal",
                              B = B, B_perm = B_perm)
cat("\nResults (Kernel column should rise well above 0.05 if it has power):\n")
print(res_power)

# ================================================================
# SECTION 4: SUMMARY TABLE
# ================================================================

cat("\n", paste(rep("=", 70), collapse=""), "\n", sep="")
cat("  SUMMARY: Sig1, dist=normal, n=20, p=1000\n")
cat(paste(rep("=", 70), collapse=""), "\n\n", sep="")
cat(sprintf("%-6s | %-6s | %-6s | %-6s | %-6s | %-6s\n",
            "Type","New","CQ","PD","HSD","Kernel"))
cat(paste(rep("-", 70), collapse=""), "\n", sep="")
cat(sprintf("%-6s | %6.3f | %6.3f | %6.3f | %6.3f | %6.3f\n",
            "SIZE",  res_size["New"],  res_size["CQ"],
                     res_size["PD"],   res_size["HSD"],  res_size["Kernel"]))
cat(sprintf("%-6s | %6.3f | %6.3f | %6.3f | %6.3f | %6.3f\n",
            "POWER", res_power["New"], res_power["CQ"],
                     res_power["PD"],  res_power["HSD"], res_power["Kernel"]))
