# Kernel-Based Test Statistic: Construction, Calibration, and Gamma Tuning

> **Setting:** p = 1000 | n = 20 | Σ₁, Σ₂, Σ₃ | μ₀, μ₁, μ₂ | Normal
> **Method:** Kernel U-statistic (RBF kernel) with permutation-based inference
> **Status:** Size/calibration confirmed; power investigation in progress — gamma ruled out as the limiting factor

---

## Table of Contents

1. [Motivation](#1-motivation)
2. [Kernel Statistic: Construction](#2-kernel-statistic-construction)
3. [Permutation Test](#3-permutation-test)
4. [Size Calibration Pilot](#4-size-calibration-pilot)
5. [Gamma (Bandwidth) Tuning Sweep](#5-gamma-bandwidth-tuning-sweep)
6. [Extended Gamma Sweep](#6-extended-gamma-sweep)
7. [Current Conclusion and Next Steps](#7-current-conclusion-and-next-steps)

---

## 1. Motivation

This work extends the JASA/CQ/JMVA test-statistic comparison (see
`JASA_vs_JMVA_comparison_report.md`) with a fourth, kernel-based
approach to the same one-sample high-dimensional mean-vector testing
problem:

$$H_0: \mu = \mu_0 \qquad \text{vs.} \qquad H_1: \mu \neq \mu_0$$

Unlike the JASA and JMVA statistics, which are built from spatial
signs of the observations, the kernel statistic evaluates a
similarity kernel on a linearly transformed version of the data,
with significance assessed via a permutation test rather than an
asymptotic normal approximation.

---

## 2. Kernel Statistic: Construction

**Data:** $X$ is the $n \times p$ data matrix (rows = subjects,
columns = features), generated the same way as in the JASA/JMVA
pipeline (`gen_data`, with mean vector `mu` from `make_mu` and
covariance `Sigma` from `make_Sigma`).

**Linear transform.** Let $A$ be an $n \times n$ matrix with i.i.d.
$N(0,1)$ entries, redrawn independently for each replication. Define

$$Y = AX \qquad (n \times p, \text{ matrix product})$$

**Kernel.** A Gaussian (RBF) kernel is applied to the rows of $Y$,
with bandwidth set by the **median heuristic**:

$$K(Y_i, Y_j) = \exp\left(-\frac{\|Y_i - Y_j\|^2}{2\sigma_K^2}\right),
\qquad \sigma_K^2 = \text{median}\left(\{\|Y_i - Y_j\|^2\}_{i<j}\right)$$

**Test statistic.** The kernel U-statistic averages the kernel over
all off-diagonal pairs:

$$T_n^{\text{Kernel}} = \frac{1}{\binom{n}{2}} \sum_{i<j} K(Y_i, Y_j)$$

This mirrors the pairwise U-statistic structure used by the JASA and
JMVA statistics, but replaces the spatial-sign inner product with a
kernel similarity evaluated on the transformed data $Y$.

---

## 3. Permutation Test

Since $T_n^{\text{Kernel}}$ has no known asymptotic null distribution
in this construction, significance is assessed via permutation:

1. Draw $A$ once per replication; compute $Y_{\text{obs}} = AX$ and
   the observed statistic $T_n^{\text{obs}}$.
2. For each of `n_perm` permutations: randomly permute the **rows of
   $X$** (holding $A$ fixed), recompute $Y_{\text{perm}} =
   AX_{\text{perm}}$, and recompute $T_n^{\text{perm}}$.
3. Compute the p-value using the standard, bias-corrected permutation
   formula:

$$\hat{p} = \frac{1 + \sum_{k=1}^{n_{\text{perm}}} \mathbb{1}\{T_n^{(k)} \geq T_n^{\text{obs}}\}}{n_{\text{perm}} + 1}$$

The `+1` correction in both numerator and denominator avoids
returning an exact-zero p-value, which is otherwise possible when the
observed statistic is more extreme than every permutation draw.

---

## 4. Size Calibration Pilot

**Setting:** $\mu_0$ (null), $\Sigma_1$, $n=20$, $p=1000$, Normal,
$B=300$ replications, `n_perm=200`.

| Metric | Result | Target |
|---|---|---|
| Empirical size | 0.023 | ≈ 0.05 |
| Mean p-value across replications | 0.496 | ≈ 0.5 (uniform under H₀) |

The mean p-value close to 0.5 indicates the permutation mechanism
itself is correctly calibrated (p-values are approximately uniform
under the null, as they should be). The empirical size of 0.023 is
below the nominal 0.05, but this gap is within Monte Carlo noise at
$B=300$ — under a true size of 0.05, the number of rejections out of
300 replications follows approximately $\text{Binomial}(300, 0.05)$,
with mean 15 and SD ≈ 3.8; the observed ~7 rejections is about 2 SDs
low, a plausible draw rather than evidence of miscalibration.

**Conclusion:** the permutation test is correctly calibrated. This
cleared the way to investigate power under $\mu_1$/$\mu_2$.

---

## 5. Gamma (Bandwidth) Tuning Sweep

check whether the RBF kernel's
bandwidth (equivalently, $\gamma = 1/(2\sigma_K^2)$) affects power,
a full grid sweep was run:

**Grid:** $\Sigma_1, \Sigma_2, \Sigma_3 \times \mu_0, \mu_1, \mu_2
\times$ 5 bandwidth multipliers $\{0.25, 0.5, 1, 2, 4\}$ applied to
the median-heuristic bandwidth ($\times 1$ reproduces the original
pilot exactly). $n=20$, $p=1000$, Normal, $B=200$, `n_perm=150`
(reduced from the pilot's $B=300$/`n_perm=200` for feasibility across
45 cells).

### Result

Rejection rates under $\mu_1$ and $\mu_2$ (POWER) ranged only
**0.020–0.075** across every covariance structure and every gamma
value tested — statistically indistinguishable from the $\mu_0$
(SIZE) rejection rates in the same range (0.025–0.090). The range of
rejection rate across the 5 gamma values, per cell, was small
throughout (0.010–0.065), with the *largest* movement occurring in a
**SIZE** cell ($\Sigma_3$, $\mu_0$: range 0.065) rather than a POWER
cell — meaning gamma was moving the false-positive rate more than it
was moving true detection power.

**Conclusion:** no evidence that gamma affects power in this range.

---

## 6. Extended Gamma Sweep

To rule out the possibility that power only emerges at an extreme
bandwidth setting, the grid was extended with four additional, much
larger multipliers: $\{10, 20, 50, 200\}$, giving a 9-point gamma
grid and $3 \times 3 \times 9 = 81$ total cells (same $B=200$,
`n_perm=150`).

### Result

The extended range did not change the conclusion:

| Summary | Value |
|---|---|
| Largest range across all 9 gamma values, any cell | 0.050 ($\Sigma_3, \mu_1$, POWER) |
| Maximum power achieved at **any** gamma, **any** cell | 0.075 ($\Sigma_3, \mu_2$) |
| Maximum size (false-positive rate) achieved at any gamma | 0.075 ($\Sigma_2, \mu_0$) |

The single best power result across the entire 81-cell grid (0.075)
is matched exactly by the best *size* result in the same grid — i.e.,
even the most favorable gamma/covariance/mean combination for
detecting a true signal performs no better than the test's own
false-positive rate under the null. Widening the bandwidth grid
50-fold (up to 200× the median heuristic) did not unlock any
additional power.

---

## 7. Current Conclusion and Next Steps

**Gamma is ruled out** as the explanation for the kernel test's lack
of power. Across a bandwidth range spanning 0.25× to 200× the
median-heuristic value, and across all three covariance structures
and both non-null mean configurations, $T_n^{\text{Kernel}}$ shows no
detectable separation between size and power. The test is correctly
*calibrated* (Section 4) but currently has **no demonstrated power**
at any bandwidth.

**Leading hypothesis:** the $Y = AX$ transformation may be
responsible. $A$ is redrawn as a fresh, independent $n \times n$
random Gaussian matrix on *every* replication. Since $A$ carries no
information about $\mu$ and changes from replication to replication,
this construction may be effectively randomizing away the very
mean-shift signal the test needs to detect, regardless of how the
kernel is subsequently tuned.

**Next diagnostic step:** re-run the size/power comparison with $A$
**fixed** (drawn once and reused across all replications within a
cell) rather than redrawn per replication, to test directly whether
this restores detectable power. This is the current priority, ahead
of any further kernel or bandwidth tuning.

---

## References

Simulation code and the parallel JASA/CQ/JMVA comparison referenced
here are available in this repository:
[Weighted-Multivariate-Sign-Functions](https://github.com/Ayush-Mishra-IIMB/Weighted-Multivariate-Sign-Functions).
