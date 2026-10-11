# LLM Benchmark-to-Target Distribution-Shift Study

## Research log, statistical design, results, and interpretation

**Status:** Research log through 12 October 2026\
**Benchmark:** GSM8K\
**Models:** GPT-OSS 120B, DeepSeek V4 Pro, Llama 3.3 70B Instruct Turbo

> **Reporting convention:** Empirical model evaluations and the separate
> known-truth simulation are reported separately. Simulation assumptions
> are not presented as measured model accuracies.

## 1. Executive summary

This project studies how estimates of LLM accuracy behave when the
benchmark question distribution differs from a target distribution. Work
to date includes a 4,000-question GSM8K benchmark/source set, a
900-question held-out target pool, a 100-question direct-target sample,
validation of three callable models, empirical benchmark-versus-target
comparisons, and a controlled repeated-sampling simulation.

Frozen model IDs:

  -------------------------------------------------------------------------------
  Model                               Together AI model ID
  ----------------------------------- -------------------------------------------
  GPT-OSS 120B                        `openai/gpt-oss-120b`

  DeepSeek V4 Pro                     `deepseek-ai/DeepSeek-V4-Pro-0813`

  Llama 3.3 70B Instruct Turbo        `meta-llama/Llama-3.3-70B-Instruct-Turbo`
  -------------------------------------------------------------------------------

**Direct-target sample (100 questions):** GPT-OSS 99/100 (99%), DeepSeek
97/100 (97%), Llama 99/100 (99%).\
**Benchmark comparison sample (100 questions):** GPT-OSS 90/100 (90%),
DeepSeek 96/100 (96%), Llama 96/100 (96%).\
All 600 model-question evaluations across these two samples were
reported as API-successful and parseable.

A separate simulation used assumed stratum-specific correctness rates to
create a known-truth scenario. In that simulation, a benchmark-only
estimator was biased for target accuracy, while direct target sampling
was nearly unbiased and had near-nominal 95% coverage but wider
intervals. These are simulation results, not measured accuracies of the
three models.

## 2. Research objective and estimands

The objective is to evaluate statistical inference for LLM accuracy
under benchmark-to-target distribution shift.

Let: - P_B = benchmark question distribution. - P_T = target question
distribution. - Y_m(x) = operational correctness outcome for model m on
question x.

The primary estimand is target accuracy:

**Target accuracy:** θ_m\^T = E\[Y_m(X) \| X drawn from P_T\].

For a fully evaluated finite target frame T containing N_T questions,
the reference value is:

**Finite-frame target accuracy:** θ_m,T = (1 / N_T) × Σ from j = 1 to
N_T of Y_m(X_j\^T).

The first controlled shift uses question-length strata. Length is an
observable covariate and a convenient shift variable, but it is not a
complete measure of difficulty or a guarantee of real-world distribution
shift. Conclusions must be scoped to the constructed population.

### Research questions

1.  How biased is a benchmark-only estimate when the target question
    distribution differs from the benchmark?
2.  How does a small direct-target sample affect bias, RMSE, 95%
    confidence-interval coverage, and interval width?
3.  How do results change with target sample size and shift severity?

## 3. Dataset construction and roles

  ------------------------------------------------------------------------
  Dataset / symbol                              Size Role
  --------------------- ---------------------------- ---------------------
  Benchmark/source set               4,000 questions Benchmark population
  B                                                  

  Held-out target pool                 900 questions Target evaluation
  T                                                  pool/reference set in
                                                     the design

  Direct-target sample                 100 questions Sample used for
  S_T                                                reported
                                                     direct-target
                                                     accuracy

  Benchmark comparison                 100 questions Sample drawn from the
  sample                                             4,000-question
                                                     benchmark
  ------------------------------------------------------------------------

The 900-question held-out pool and 100-question direct-target sample are
not interchangeable. The model accuracy tables below are based on 100
evaluated questions per model per sample. Do not describe the
900-question pool as 900 completed model evaluations unless it was
separately scored in full.

### 3.1 Direct-target sample composition

  Stratum         Count      Share
  ----------- --------- ----------
  Short              15        15%
  Medium             35        35%
  Long               50        50%
  **Total**     **100**   **100%**

Recorded sample seed: `20261012`. The recorded SHA-256 fingerprint
begins `c28a88e` and ends `805e0f1`; use the full fingerprint in the
project manifest for formal verification.

### 3.2 Benchmark comparison sample composition

  Stratum         Count      Share
  ----------- --------- ----------
  Short              28        28%
  Medium             48        48%
  Long               24        24%
  **Total**     **100**   **100%**

The recorded SHA-256 fingerprint begins `a14be3` and ends `e2000f`.
These are custom length groups, not official GSM8K labels. Because the
benchmark and target samples have different stratum allocations, raw
accuracy differences can reflect both model performance and sample
composition.

### 3.3 Data-integrity requirements

Preserve question IDs, seeds, strata, sampling probabilities, and
hashes. Check duplicates and intended disjointness among benchmark,
direct-target, and held-out target sets. Validate numerical ground truth
before scoring. Do not silently remove or manually repair outcomes after
seeing model performance.

## 4. Model selection and technical validation

All three models were found in the Together AI catalog and returned
nonempty responses in a basic API preflight. Each model was then tested
on three known-answer arithmetic prompts, with two repetitions per
prompt:

-   17 + 25 = 42
-   144 ÷ 12 = 12
-   23 × 17 = 391

  ------------------------------------------------------------------------
  Model             Calls     Nonempty      Correct       Errors Status
                             responses      answers              
  ---------- ------------ ------------ ------------ ------------ ---------
  GPT-OSS               6            6            6            0 Passed
  120B                                                           

  DeepSeek              6            6            6            0 Passed
  V4 Pro                                                         

  Llama 3.3             6            6            6            0 Passed
  70B                                                            
  ------------------------------------------------------------------------

This establishes basic callability and elementary answer formatting, not
benchmark-level reliability.

### Primary-run settings recorded

-   Frozen model IDs as listed above.
-   Temperature: `0`.
-   Direct-target sample: 100 shared questions; 300 model-question
    evaluations.
-   API retries: up to 2 according to the run log.
-   GPT-OSS and DeepSeek used the Responses API; Llama used Chat
    Completions.
-   Saved outcomes included response, API status, parser status,
    predicted numerical answer, ground truth, and operational
    correctness.
-   A separate benchmark comparison used 100 questions per model,
    another 300 evaluations.
-   All 600 evaluations across the two runs were reported as
    API-successful and parseable.

## 5. Scoring definition

For model m and question x, operational correctness is:

-   Y_m(x) = 1 if the model returns a usable, parseable final answer
    that matches the numerical ground truth.
-   Y_m(x) = 0 if the answer is incorrect, the API call fails, or
    parsing fails.

This keeps every assigned question in the primary denominator. API
failures and unparseable answers count as operational failures rather
than being silently dropped.

Secondary conditional accuracy is:

**Conditional accuracy** = (number correct and parseable) / (number of
successful, parseable outcomes).

Always report the conditional denominator. Conditional accuracy is not
the same as operational accuracy.

## 6. Accuracy and Wilson confidence intervals

For n assigned questions and k operationally correct answers:

**Sample accuracy:** p̂ = k / n.

The 95% Wilson interval uses z = 1.96. Define:

-   Denominator: D = 1 + z²/n
-   Center: C = (p̂ + z²/(2n)) / D
-   Half-width: H = \[z / D\] × √\[p̂(1 − p̂)/n + z²/(4n²)\]

Then:

**Wilson 95% confidence interval:** \[C − H, C + H\].

Equivalently, the interval endpoints are:

Lower = {p̂ + z²/(2n) − z√\[p̂(1 − p̂)/n + z²/(4n²)\]} / {1 + z²/n}

Upper = {p̂ + z²/(2n) + z√\[p̂(1 − p̂)/n + z²/(4n²)\]} / {1 + z²/n}.

Wilson intervals are used for the descriptive 100-question results
below. For stratified samples or finite-population sampling without
replacement, uncertainty procedures should reflect the sampling design
where relevant.

## 7. Empirical results: direct-target sample

  -------------------------------------------------------------------------------------
  Model         Correct          n   Accuracy Wilson 95% CI           API   Unparseable
                                                                 failures 
  ---------- ---------- ---------- ---------- ---------------- ---------- -------------
  GPT-OSS            99        100        99% 94.55%--99.82%            0             0
  120B                                                                    

  DeepSeek           97        100        97% 91.55%--98.97%            0             0
  V4 Pro                                                                  

  Llama 3.3          99        100        99% 94.55%--99.82%            0             0
  70B                                                                     
  -------------------------------------------------------------------------------------

All 300 direct-target evaluations were reported as successful and
parseable.

**Interpretation:** GPT-OSS and Llama both scored 99/100; DeepSeek
scored 97/100. The intervals overlap substantially, so these descriptive
results alone do not establish a population-level difference. Since the
models answered the same questions, formal comparisons should preserve
question-level pairing.

## 8. Empirical results: benchmark comparison sample

  -------------------------------------------------------------------------------------
  Model         Correct          n   Accuracy Wilson 95% CI           API   Unparseable
                                                                 failures 
  ---------- ---------- ---------- ---------- ---------------- ---------- -------------
  GPT-OSS            90        100        90% 82.56%--94.48%            0             0
  120B                                                                    

  DeepSeek           96        100        96% 90.16%--98.43%            0             0
  V4 Pro                                                                  

  Llama 3.3          96        100        96% 90.16%--98.43%            0             0
  70B                                                                     
  -------------------------------------------------------------------------------------

All 300 benchmark-comparison evaluations were reported as successful and
parseable.

### 8.1 Raw target-minus-benchmark differences

Define Δ_m = p̂\_m,T − p̂\_m,B.

  Model               Benchmark sample   Direct-target sample             Difference
  ----------------- ------------------ ---------------------- ----------------------
  GPT-OSS 120B                     90%                    99%   +9 percentage points
  DeepSeek V4 Pro                  96%                    97%    +1 percentage point
  Llama 3.3 70B                    96%                    99%   +3 percentage points

These are descriptive differences between two samples with different
length-group allocations. They should not be interpreted as causal
effects of distribution shift without accounting for sample design,
target construction, and uncertainty.

## 9. Estimators for the distribution-shift study

Let p_T,h be the known target share in stratum h, and θ̂\_m,h\^B the
benchmark estimate of accuracy in that stratum.

### Method 1 --- Naive benchmark-only

**Estimator:** θ̂\_m,1 = (1/n_B) × Σ from i = 1 to n_B of Y_m(X_i\^B).

This is unweighted benchmark accuracy. Under shift, it need not estimate
target accuracy.

### Method 2 --- Benchmark post-stratification

**Estimator:** θ̂\_m,2 = Σ from h = 1 to H of p_T,h × θ̂\_m,h\^B.

This reweights benchmark stratum means to target stratum shares. It
relies on conditional transportability: within each stratum, benchmark
outcomes must be informative about target outcomes.

### Method 3 --- Direct target sample

**Estimator:** θ̂\_m,3 = (1/n_T) × Σ from i = 1 to n_T of Y_m(X_i\^T).

Under probability sampling from the target population, this directly
estimates target accuracy. A small target sample can have high variance.

### Method 4 --- Combined benchmark and target estimator

Within each stratum, combine benchmark and target observations using
prespecified weights that respect the sampling design, then weight the
stratum means by target shares:

**Estimator:** θ̂\_m,4 = Σ from h = 1 to H of p_T,h × θ̂\_m,h\^(combined).

This may improve precision by using the larger benchmark sample while
targeting the target distribution. The pooling weights and variance
formula must be fixed before confirmatory results are inspected.
Post-stratified and combined simulation results must be computed, not
inferred from theoretical expectations.

## 10. Controlled-shift simulation settings

This simulation is separate from the empirical model runs. It uses
assumed correctness rates to create a known-truth setting, then
repeatedly samples from it to evaluate statistical procedures.

  -----------------------------------------------------------------------
  Setting                             Value
  ----------------------------------- -----------------------------------
  Benchmark population size N_B       4,000

  Direct target sample size n_T       100

  Repeated sampling replicates R      2,000

  Simulation seed                     `20261010`

  Confidence level                    95%

  Simple-proportion interval          Wilson

  Benchmark weights                   0.27650 / 0.48525 / 0.23825
  (Short/Medium/Long)                 

  Target weights (Short/Medium/Long)  0.15 / 0.35 / 0.50

  Assumed correctness by stratum      0.90 / 0.75 / 0.55
  (Short/Medium/Long)                 

  Main metrics                        Coverage, bias, RMSE, average
                                      interval width
  -----------------------------------------------------------------------

The correctness rates 0.90, 0.75, and 0.55 are simulation assumptions,
not empirical results for any LLM.

### 10.1 Known benchmark and target accuracy

Benchmark accuracy:

θ_B = 0.27650 × 0.90 + 0.48525 × 0.75 + 0.23825 × 0.55 = 0.743825.

Target accuracy:

θ_T = 0.15 × 0.90 + 0.35 × 0.75 + 0.50 × 0.55 = 0.6725.

Target-minus-benchmark difference:

δ = θ_T − θ_B = 0.6725 − 0.743825 = −0.071325.

Thus, in this simulated scenario, target accuracy is about 7.13
percentage points lower than benchmark accuracy because the target
places more weight on the longest-question stratum, whose assumed
correctness rate is lowest.

### 10.2 Repeated-sampling metrics

For replicate r = 1, ..., R, let θ̂\_r be an estimate and θ_T the known
target reference.

**Bias:**

Bias_hat = (1/R) × Σ from r = 1 to R of (θ̂\_r − θ_T).

**Root mean squared error (RMSE):**

RMSE_hat = √{(1/R) × Σ from r = 1 to R of (θ̂\_r − θ_T)²}.

**Empirical 95% coverage:**

Coverage_hat = (1/R) × Σ from r = 1 to R of I(θ_T lies between L_r and
U_r).

Here, I(condition) is 1 when the condition is true and 0 otherwise.

**Average interval width:**

W_bar = (1/R) × Σ from r = 1 to R of (U_r − L_r).

Coverage is a repeated-sampling property. It cannot be established from
one confidence interval in one empirical run.

### 10.3 Recorded simulation results

  Metric                     Naive benchmark-only   Direct-target sample
  ------------------------ ---------------------- ----------------------
  Bias                                    +0.0716                +0.0004
  RMSE                                    0.07197                0.04607
  Empirical 95% coverage                    0.000                  0.947
  Average interval width                  0.02703                0.17997

**Interpretation:**

-   The naive benchmark-only estimator has positive bias of about 7.16
    percentage points because it estimates the higher benchmark accuracy
    rather than target accuracy. Its intervals are narrow but had zero
    coverage in the recorded simulation.
-   The direct-target estimator has bias close to zero and coverage
    about 94.7%, near the nominal 95%, but much wider intervals because
    it uses only 100 target observations.
-   This illustrates the bias--precision trade-off: a narrow interval
    around the wrong estimand is not a good estimate of target accuracy.
-   The post-stratified and combined methods still need to be evaluated
    and reported under the same simulation design.

With R = 2,000 and nominal coverage 0.95, the Monte Carlo standard error
of coverage is approximately:

SE_MC = √\[0.95 × (1 − 0.95) / 2,000\] ≈ 0.0049,

or about 0.49 percentage points.

## 11. Stratified confidence intervals

For a stratified estimator θ̂\_strat = Σ from h = 1 to H of p_T,h × θ̂\_h,
a common variance approximation for independent stratum samples is:

Var_hat(θ̂\_strat) = Σ from h = 1 to H of p_T,h² × θ̂\_h(1 − θ̂\_h) / n_h.

The standard error is SE_hat = √Var_hat. A basic 95% normal interval is:

θ̂\_strat ± 1.96 × SE_hat.

Because this interval can extend outside \[0, 1\], the implementation
should use a bounded interval or a prespecified stratified bootstrap
where appropriate. For sampling without replacement from a finite
population, include the finite-population correction 1 − n_h/N_h when
justified by the sampling design. The final interval method should be
fixed before confirmatory analysis.

## 12. Paired model comparisons

The three models are evaluated on the same questions within each sample.
Preserve this pairing. For models a and b, define the question-level
difference:

D_i = Y_a(X_i) − Y_b(X_i), so D_i is −1, 0, or 1.

The paired accuracy difference is:

Δ̂\_a,b = (1/n) × Σ from i = 1 to n of D_i.

An exact McNemar test can compare paired binary outcomes using
discordant pairs. With three models, there are three pairwise
comparisons; use a prespecified multiplicity correction such as Holm's
procedure if all are tested. Overlapping Wilson intervals are not a
replacement for a paired test, and a nonsignificant test does not prove
equivalence.

## 13. What the results support

### Supported by completed empirical runs

-   The three selected models passed basic API and arithmetic
    validation.
-   Direct-target sample accuracy was 99% for GPT-OSS, 97% for DeepSeek,
    and 99% for Llama.
-   Benchmark comparison sample accuracy was 90% for GPT-OSS, 96% for
    DeepSeek, and 96% for Llama.
-   All 600 evaluations across the two empirical samples were reported
    successful and parseable.
-   The two samples had different length-stratum allocations.

### Supported by the recorded simulation

-   Benchmark and target accuracies can differ under a controlled change
    in stratum weights.
-   A naive benchmark estimator can be biased for target accuracy.
-   Direct target sampling can provide near-nominal coverage but wider
    intervals at small sample sizes.

### Not yet established

-   That any model is universally superior.
-   That length shift captures all real-world distribution shift.
-   That empirical accuracy differences are caused solely by
    distribution shift.
-   That post-stratification or combined estimation is superior until
    their simulation results are computed.
-   That the simulation's assumed correctness probabilities describe the
    three LLMs.

## 14. Reproducibility and artifact log

Known result artifacts:

  ------------------------------------------------------------------------------
  Artifact                                   Purpose
  ------------------------------------------ -----------------------------------
  `gsm8k_primary_results.csv`                300 outcomes: 100 direct-target
                                             questions × 3 models

  `gsm8k_benchmark_comparison_results.csv`   300 outcomes: 100 benchmark
                                             questions × 3 models

  Target sample manifest                     Target question IDs, strata, seed,
                                             and fingerprint

  Benchmark sample manifest                  Benchmark comparison question IDs,
                                             strata, and fingerprint

  Simulation configuration/results           Assumptions, seed, repetitions, and
                                             estimator metrics
  ------------------------------------------------------------------------------

Confirm exact on-disk paths in the active project before archiving or
moving files. Preserve the earlier exploratory run with DeepSeek
failures as a historical artifact; do not combine it with the newer
validated results.

Suggested structure:

``` text
LLM STATS AI/
├── data/
│   ├── raw/benchmark_B.csv
│   ├── processed/
│   │   ├── gsm8k_primary_results.csv
│   │   └── gsm8k_benchmark_comparison_results.csv
│   └── manifests/
├── notebooks/
│   ├── model_preflight_validation.ipynb
│   ├── empirical_accuracy_comparison.ipynb
│   └── distribution_shift_simulation.ipynb
└── reports/
    └── llm_benchmark_target_shift_research_log.md
```

## 15. Overall interpretation and next steps

The empirical runs show high observed accuracy on the sampled GSM8K
questions, but the model ranking and accuracy levels differ between the
benchmark comparison sample and the direct-target sample. Benchmark
accuracy should not automatically be treated as target-population
accuracy.

The known-truth simulation illustrates the statistical problem directly:
in the specified scenario, the naive benchmark estimator is biased for
target accuracy and has zero coverage, while the direct-target estimator
is nearly unbiased with approximately nominal coverage but wider
intervals.

The next step is to implement and validate the repeated-sampling
comparison of benchmark post-stratification and the combined estimator,
using the same known-truth target population and all four evaluation
criteria. More complex methods such as importance weighting, doubly
robust estimation, or an LLM-generated target population should be added
only after the baseline methods are validated.

## Appendix A. Frozen model configuration

``` python
STUDY_MODEL_IDS = {
    "GPT-OSS": "openai/gpt-oss-120b",
    "DeepSeek": "deepseek-ai/DeepSeek-V4-Pro-0813",
    "Llama": "meta-llama/Llama-3.3-70B-Instruct-Turbo",
}
```

## Appendix B. Key terms

  -----------------------------------------------------------------------
  Term                                Meaning
  ----------------------------------- -----------------------------------
  Benchmark accuracy                  Expected operational correctness
                                      under P_B

  Target accuracy                     Expected operational correctness
                                      under P_T

  Operational accuracy                Correct usable answers divided by
                                      all assigned questions

  Conditional accuracy                Correct answers among successful
                                      parseable outcomes

  Bias                                Average estimate minus the
                                      reference value

  RMSE                                Square root of average squared
                                      estimation error

  95% coverage                        Fraction of repeated-sampling
                                      intervals containing the reference
                                      value

  Interval width                      Upper endpoint minus lower endpoint

  Covariate shift                     Change in observed covariate
                                      distribution under a conditional
                                      transportability assumption

  Direct-target sample                Probability sample drawn from the
                                      target population

  Finite-frame reference              Accuracy over every item in the
                                      defined finite target frame, when
                                      outcomes are available
  -----------------------------------------------------------------------
