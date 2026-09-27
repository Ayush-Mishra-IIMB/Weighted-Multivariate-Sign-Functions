# ================================================================
#  KERNEL T_n — RBF GAMMA (BANDWIDTH) SWEEP
#  Full grid: Sigma1/2/3 x mu0/mu1/mu2, n=20, p=1000, Normal
#
#  Goal: check whether varying the RBF kernel's gamma (equivalently
#  the bandwidth sigma_K^2 = 1/(2*gamma)) changes power in a
#  meaningful way. If power moves across gamma, that's a signal the
#  bandwidth matters and is worth tuning properly; if it's flat,
#  the median-heuristic default is probably fine as-is.
#
#  Kept identical to the earlier calibrated pilot except for one
#  change: instead of always using the median-heuristic bandwidth,
#  sigma_K^2 is set to a MULTIPLE of the median heuristic, i.e.
#      sigma_K^2 = gamma_mult * median(pairwise squared distances)
#  gamma_mult = 1 reproduces the original median-heuristic pilot
#  exactly, so that run is our known-calibrated reference point.
# ================================================================

library(MASS)
set.seed(2025)

# ---- Sweep settings ----
GAMMA_MULT <- c(0.25, 0.5, 1, 2, 4)   # multiples of median heuristic bandwidth
                                       # (gamma_mult=1 == original pilot)
B      <- 200   # replications per cell (reduced from 600 for scan feasibility;
                 # this is a screening scan, not a final-precision run)
N_PERM <- 150   # permutations per replication (reduced from 200 for the same reason)

N      <- 20
P      <- 1000
DIST   <- "normal"

SIGMA_TYPES <- 1:3
MU_TYPES    <- 0:2

cat("================================================================\n")
cat("  KERNEL Tn — RBF GAMMA (BANDWIDTH) SWEEP\n")
cat("  n=20 | p=1000 | Normal | Sigma1/2/3 x mu0/mu1/mu2\n")
cat("  gamma_mult grid:", paste(GAMMA_MULT, collapse=", "), "\n")
cat("  B =", B, "| N_PERM =", N_PERM, " (screening run; not final precision)\n")
cat("================================================================\n\n")

# ================================================================
# SECTION 1: REUSED FUNCTIONS (identical to prior pipeline)
# ================================================================

make_Sigma <- function(p, type) {
  if (type == 1) {
    S <- matrix(0.2, p, p); diag(S) <- 1; return(S)
  }
  if (type == 2) {
    return(outer(1:p, 1:p, function(i, j) 0.8^abs(i - j)))
  }
  if (type == 3) {
    d <- 2 + (p - 1:p + 1) / p
    R <- outer(1:p, 1:p, function(i, j)
      ifelse(i == j, 1, (-1)^(i + j) * 0.2^(abs(i - j) / 0.1)))
    return(diag(d) %*% R %*% diag(d))
  }
}

make_mu <- function(p, type) {
  if (type == 0) return(rep(0, p))
  if (type == 1) return(rep(0.25, p))
  if (type == 2) {
    mu <- rep(0, p)
    mu[(floor(p / 3) + 1):floor(2 * p / 3)] <- 0.25
    mu[(floor(2 * p / 3) + 1):p]            <- -0.25
    return(mu)
  }
}

gen_data <- function(n, mu, Sigma, dist) {
  p <- length(mu)
  if (dist == "normal")
    return(mvrnorm(n, mu = mu, Sigma = Sigma))
}

generate_A <- function(n) {
  matrix(rnorm(n * n), n, n)     # iid N(0,1), n x n
}

make_Y <- function(X, A) {
  A %*% X                        # matrix multiplication, n x p
}

# ================================================================
# SECTION 2: MODIFIED — RBF KERNEL WITH TUNABLE GAMMA MULTIPLIER
# ================================================================
#
# sigma_K^2 = gamma_mult * median(pairwise squared distances)
# gamma_mult = 1 reproduces the original median-heuristic pilot exactly.

rbf_kernel_matrix <- function(Y, gamma_mult = 1) {
  D2 <- as.matrix(dist(Y))^2
  sigma2 <- gamma_mult * median(D2[upper.tri(D2)])
  if (sigma2 == 0) sigma2 <- 1e-6
  exp(-D2 / (2 * sigma2))
}

kernel_Tn <- function(Y, gamma_mult = 1) {
  Kmat <- rbf_kernel_matrix(Y, gamma_mult)
  n    <- nrow(Y)
  s    <- sum(Kmat[upper.tri(Kmat)])
  s / choose(n, 2)
}

# ================================================================
# SECTION 3: PERMUTATION TEST (corrected p-value formula)
# ================================================================

kernel_perm_pval <- function(X, gamma_mult = 1, n_perm = N_PERM) {
  n <- nrow(X)
  A <- generate_A(n)
  Y_obs  <- make_Y(X, A)
  Tn_obs <- kernel_Tn(Y_obs, gamma_mult)

  Tn_perm <- numeric(n_perm)
  for (k in 1:n_perm) {
    perm_idx <- sample(1:n)
    X_perm   <- X[perm_idx, , drop = FALSE]
    Y_perm   <- make_Y(X_perm, A)
    Tn_perm[k] <- kernel_Tn(Y_perm, gamma_mult)
  }

  # corrected p-value formula (avoids exact-zero p-values)
  p_val <- (1 + sum(Tn_perm >= Tn_obs)) / (n_perm + 1)
  p_val
}

# ================================================================
# SECTION 4: SWEEP — ALL SIGMA x MU x GAMMA_MULT COMBINATIONS
# ================================================================

run_cell <- function(n, p, mu_type, sigma_type, dist, gamma_mult, B, n_perm) {
  mu    <- make_mu(p, mu_type)
  Sigma <- make_Sigma(p, sigma_type)

  rej <- 0
  for (b in 1:B) {
    X <- gen_data(n, mu, Sigma, dist)
    p_val <- kernel_perm_pval(X, gamma_mult = gamma_mult, n_perm = n_perm)
    if (p_val < 0.05) rej <- rej + 1
  }
  rej / B
}

results <- data.frame()

total_cells <- length(SIGMA_TYPES) * length(MU_TYPES) * length(GAMMA_MULT)
cell_count  <- 0

for (sigma_type in SIGMA_TYPES) {
  for (mu_type in MU_TYPES) {
    for (gm in GAMMA_MULT) {
      cell_count <- cell_count + 1
      cat(sprintf("[%d/%d] Sigma%d, mu%d, gamma_mult=%.2f ... ",
                   cell_count, total_cells, sigma_type, mu_type, gm))

      rate <- run_cell(N, P, mu_type, sigma_type, DIST, gm, B, N_PERM)

      cat(sprintf("rate = %.3f\n", rate))

      results <- rbind(results, data.frame(
        Sigma = sigma_type,
        mu    = mu_type,
        gamma_mult = gm,
        rejection_rate = rate,
        type  = ifelse(mu_type == 0, "SIZE", "POWER")
      ))
    }
  }
}

cat("\n================================================================\n")
cat("  FULL SWEEP RESULTS\n")
cat("================================================================\n\n")
print(results, row.names = FALSE)

# ================================================================
# SECTION 5: SUMMARY — DOES POWER MOVE WITH GAMMA?
# ================================================================

cat("\n================================================================\n")
cat("  RANGE OF REJECTION RATE ACROSS GAMMA, PER (SIGMA, MU) CELL\n")
cat("  (large range => gamma matters; near-zero range => flat / doesn't matter)\n")
cat("================================================================\n\n")

summary_tab <- aggregate(rejection_rate ~ Sigma + mu + type, data = results,
                          FUN = function(x) max(x) - min(x))
names(summary_tab)[names(summary_tab) == "rejection_rate"] <- "range_across_gamma"
summary_tab <- summary_tab[order(-summary_tab$range_across_gamma), ]
print(summary_tab, row.names = FALSE)

cat("\nNOTE: gamma_mult = 1 reproduces the original median-heuristic\n")
cat("pilot setting exactly (our known-calibrated reference point).\n")
cat("B and N_PERM are reduced here for screening feasibility across\n")
cat("the full 3x3x5 = 45-cell grid; re-run promising cells at full\n")
cat("B=600 / N_PERM=200 to confirm before drawing final conclusions.\n")
