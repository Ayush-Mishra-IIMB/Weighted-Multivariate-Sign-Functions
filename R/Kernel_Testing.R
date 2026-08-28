# ================================================================
#  KERNEL Tn PERMUTATION TEST — full grid (3 Sigma types x 3 mu types)
#  Same setting as kernel_Tn_simple.R, extended to match the original
#  JASA/JMVA make_Sigma() / make_mu() definitions exactly.
#
#    Y = A X,  A: n x n, a_ij ~ iid N(0,1), symmetrized + shifted to be PD
#    Tn = 1/C(n,2) * sum_{i<j} K(yi,yj),  K = RBF, median-heuristic bandwidth
#    p-value: permutation (permute rows of X, A stays fixed)
#
#  Grid:
#    mu0 (null)                x Sig1, Sig2, Sig3   -> SIZE checks
#    mu1 (constant shift)      x Sig1, Sig2, Sig3   -> POWER checks
#    mu2 (block +/- pattern)   x Sig1, Sig2, Sig3   -> POWER checks
#  dist = "normal" only, n = 20, p = 1000 (same as kernel_Tn_simple.R)
# ================================================================

library(MASS)
set.seed(2025)

n      <- 20
p      <- 1000
B      <- 100     # outer simulation replications
B_perm <- 100     # permutations per replication

cat("================================================================\n")
cat("  KERNEL Tn — full grid (3 Sigma types x 3 mu types)\n")
cat("  n =", n, " p =", p, " B =", B, " B_perm =", B_perm, "\n")
cat("================================================================\n\n")

# ================================================================
# SECTION 1: DATA-GENERATING SETUP (identical to JASA/JMVA pipeline)
# ================================================================

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

gen_data <- function(n, mu, Sigma) {
  mvrnorm(n, mu = mu, Sigma = Sigma)   # dist = "normal" only, as agreed
}

cat("Data-generating functions loaded (make_Sigma, make_mu, gen_data).\n\n")

# ================================================================
# SECTION 2: KERNEL Tn AND PERMUTATION P-VALUE (unchanged from before)
# ================================================================

generate_A <- function(n) {
  M <- matrix(rnorm(n * n), n, n)
  S <- (M + t(M)) / 2
  eig <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  shift <- max(0, -min(eig)) + 0.01
  S + diag(shift, n)
}

rbf_kernel_matrix <- function(Y) {
  D2 <- as.matrix(dist(Y))^2
  sigma2 <- median(D2[upper.tri(D2)])
  if (sigma2 == 0) sigma2 <- 1e-6
  exp(-D2 / (2 * sigma2))
}

kernel_Tn <- function(X, A) {
  n <- nrow(X)
  Y <- A %*% X
  K <- rbf_kernel_matrix(Y)
  sum(K[upper.tri(K)]) / choose(n, 2)
}

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
# SECTION 3: RUN ONE CELL (mu_type x sigma_type)
# ================================================================

run_cell_kernel <- function(n, p, mu_type, sigma_type, B, B_perm) {
  mu    <- make_mu(p, mu_type)
  Sigma <- make_Sigma(p, sigma_type)

  pvals <- numeric(B)
  for (b in 1:B) {
    X <- gen_data(n, mu, Sigma)
    A <- generate_A(n)   # fresh A drawn each replication
    pvals[b] <- kernel_perm_pval(X, A, B_perm)
  }

  mean(pvals < 0.05)
}

# ================================================================
# SECTION 4: RUN THE FULL 3 x 3 GRID
# ================================================================

mu_types    <- c(0, 1, 2)
sigma_types <- c(1, 2, 3)

results <- matrix(NA, nrow = length(mu_types), ncol = length(sigma_types),
                   dimnames = list(paste0("mu", mu_types), paste0("Sig", sigma_types)))

for (mt in mu_types) {
  for (st in sigma_types) {
    cat("Running mu_type =", mt, " sigma_type =", st, "...\n")
    results[paste0("mu", mt), paste0("Sig", st)] <-
      round(run_cell_kernel(n, p, mu_type = mt, sigma_type = st, B = B, B_perm = B_perm), 3)
  }
}

# ================================================================
# SECTION 5: SUMMARY TABLE
# ================================================================

cat("\n", paste(rep("=", 60), collapse=""), "\n", sep="")
cat("  KERNEL Tn — rejection rates (rows=mu type, cols=Sigma type)\n")
cat("  mu0 row = SIZE (want ~0.05) | mu1, mu2 rows = POWER (want > 0.05)\n")
cat(paste(rep("=", 60), collapse=""), "\n\n", sep="")
print(results)
