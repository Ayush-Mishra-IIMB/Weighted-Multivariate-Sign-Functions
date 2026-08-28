# ================================================================
#  KERNEL Tn PERMUTATION TEST — simple, single null vs single alternate
#  Null:        X ~ N(0, Sigma1)
#  Alternate:   X ~ N(delta, Sigma1)   <- SAME family, only mean shifts
#  Y = A X,  A: n x n, a_ij ~ iid N(0,1), symmetrized + shifted to be PD
#  Tn = 1/C(n,2) * sum_{i<j} K(yi,yj),  K = RBF, median-heuristic bandwidth
#  p-value: permutation (permute rows of X, A stays fixed)
# ================================================================

library(MASS)
set.seed(2025)

n      <- 20
p      <- 1000
B      <- 100     # outer simulation replications
B_perm <- 100     # permutations per replication
delta  <- 1       # non-zero mean shift for the alternative

# ---- Fixed covariance: Sig1 (compound symmetric), same as before ----
Sigma <- matrix(0.2, p, p)
diag(Sigma) <- 1

# ---- Null and alternate: SAME family, only the mean differs ----
gen_X_null <- function(n, p, Sigma, ...) {
  mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
}

gen_X_alt <- function(n, p, Sigma, delta) {
  mvrnorm(n, mu = rep(delta, p), Sigma = Sigma)
}

# ---- A: n x n, iid N(0,1), symmetrized + eigenvalue-shifted to be PD ----
generate_A <- function(n) {
  M <- matrix(rnorm(n * n), n, n)
  S <- (M + t(M)) / 2
  eig <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  shift <- max(0, -min(eig)) + 0.01
  S + diag(shift, n)
}

# ---- Fixed RBF kernel, median-heuristic bandwidth ----
rbf_kernel_matrix <- function(Y) {
  D2 <- as.matrix(dist(Y))^2
  sigma2 <- median(D2[upper.tri(D2)])
  if (sigma2 == 0) sigma2 <- 1e-6
  exp(-D2 / (2 * sigma2))
}

# ---- Tn = 1/C(n,2) * sum_{i<j} K(yi,yj),  Y = A X ----
kernel_Tn <- function(X, A) {
  n <- nrow(X)
  Y <- A %*% X
  K <- rbf_kernel_matrix(Y)
  sum(K[upper.tri(K)]) / choose(n, 2)
}

# ---- Permutation p-value: permute rows of X, A fixed ----
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

# ---- Run B repetitions, report rejection rate AND store all p-values ----
run_check <- function(gen_fun, n, p, Sigma, B, B_perm, ...) {
  pvals <- numeric(B)
  for (b in 1:B) {
    X <- gen_fun(n, p, Sigma, ...)
    A <- generate_A(n)
    pvals[b] <- kernel_perm_pval(X, A, B_perm)
    if (b %% 20 == 0) cat("  rep", b, "/", B, "\n")
  }
  list(rejection_rate = mean(pvals < 0.05), pvals = pvals)
}

cat("Running NULL (X ~ N(0, Sigma1))...\n")
null_res <- run_check(gen_X_null, n, p, Sigma, B, B_perm)
cat("Size (want ~0.05):", null_res$rejection_rate, "\n")
cat("Null p-values summary:\n")
print(summary(null_res$pvals))

cat("\nRunning ALTERNATE (X ~ N(delta, Sigma1), delta =", delta, ")...\n")
alt_res <- run_check(gen_X_alt, n, p, Sigma, B, B_perm, delta = delta)
cat("Power (want notably > 0.05):", alt_res$rejection_rate, "\n")
cat("Alternate p-values summary:\n")
print(summary(alt_res$pvals))

# ---- Visualize: p-value histograms under null and alternate ----
# Under H0 these should look roughly Uniform(0,1).
# Under H1, if the test has power, they should pile up near 0.
dev.off()   

# then re-run WITHOUT the png() line:
par(mfrow = c(1, 2))
hist(null_res$pvals, breaks = 20, col = "steelblue",
     main = "P-values under NULL ",
     xlab = "p-value", xlim = c(0,1))
hist(alt_res$pvals, breaks = 20, col = "darkorange",
     main = "P-values under ALTERNATE",
     xlab = "p-value", xlim = c(0,1))