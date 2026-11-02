# Statistical Inference for LLM Performance Under Benchmark-to-Target Distribution Shift

## Consolidated Research Record --- Stages 1--9

**Project:** IS T2 AI STATS\
**Empirical pilot:** GSM8K; DeepSeek-V4-Pro-0813 and GPT-OSS-120b\
**Consolidation date:** 26 October 2026

## 1. Executive summary

This project studies inference for LLM accuracy when the benchmark
distribution differs from the target population. The empirical pilot
uses Short and Long GSM8K question strata. The benchmark mixture is 62%
Short / 38% Long. The target mixture used in the pilot is 30% Short /
70% Long, a designed scenario rather than a verified real-world
population proportion.

The project retains exactly three estimators:

1.  **Benchmark-only accuracy:** benchmark stratum accuracies weighted
    by the benchmark mixture.
2.  **Reweighted benchmark accuracy:** benchmark stratum accuracies
    weighted by the assumed target mixture.
3.  **Direct target accuracy:** target stratum accuracies weighted by
    the assumed target mixture.

The central finding from simulation is that reweighting can correct
composition differences when conditional accuracy is stable across
benchmark and target populations within strata. It cannot, by itself,
correct conditional distribution shift. Direct target estimation was
approximately unbiased in the simulated designs because the target
sample was generated from the target distribution, but this relies on
target-sample representativeness.

In the empirical pilot, both models had high observed accuracy.
DeepSeek's direct-target accuracy was 100% among scored target
responses, but three target Long responses were unscored. Counting these
as incorrect gives a 97% target estimate under the assumed 30%/70%
mixture. Paired point estimates favored DeepSeek in complete pairs, but
intervals included or touched zero; statistical superiority is not
established.

The simulation results are not empirical claims about LLMs. The
simulation used synthetic data-generating mechanisms, not LLM-generated
synthetic questions. No new LLM calls were made in Stages 2--9.

## 2. Research question and motivation

**Research question:** How should LLM accuracy and its uncertainty be
estimated when the benchmark and target population differ in
question-type composition and potentially in conditional accuracy?

Overall accuracy is a weighted average of stratum accuracies. A
different target mix can change accuracy even if within-stratum
performance remains stable. Reweighting benchmark stratum accuracies to
the target mix can address this composition shift. If conditional
accuracy changes between benchmark and target within a stratum,
reweighting alone does not recover target accuracy.

## 3. Estimands and notation

Let (g `\in `{=tex}{S,L}) represent Short and Long questions. Let
(p\_{B,g}) be benchmark accuracy in stratum (g), and (p\_{T,g}) target
accuracy. Benchmark weights are (w\_{B,S}=0.62) and (w\_{B,L}=0.38).
Target weights are (w\_{T,S}=q) and (w\_{T,L}=1-q), with (q=0.30) in the
pilot.

### 3.1 Benchmark-only accuracy

\[
`\widehat`{=tex}`\theta`{=tex}*B=0.62`\widehat `{=tex}p*{B,S}+0.38`\widehat `{=tex}p\_{B,L}.
\]

This describes the benchmark mixture. It is not automatically the
target-population accuracy.

### 3.2 Reweighted benchmark accuracy

\[
`\widehat`{=tex}`\theta`{=tex}*{RW}=q`\widehat `{=tex}p*{B,S}+(1-q)`\widehat `{=tex}p\_{B,L}.
\]

This applies target weights to benchmark stratum accuracies. It
addresses composition shift under conditional stability, not arbitrary
distribution shift.

### 3.3 Direct target accuracy

\[
`\widehat`{=tex}`\theta`{=tex}*T=q`\widehat `{=tex}p*{T,S}+(1-q)`\widehat `{=tex}p\_{T,L}.
\]

This uses target-sample stratum accuracies. Its validity depends on the
target sample representing the target population and on the specified
target weights being appropriate.

**Important:** the 30% Short / 70% Long target mix is a design
assumption, not a validated real-world target composition.

## 4. Empirical pilot and data-quality audit

The existing pilot file is `analysis_outputs/llm_two_model_results.csv`.
It contains 400 rows: 100 benchmark questions and 100 target questions
evaluated by each of two models. The binary correctness outcome is
`correct`.

Three DeepSeek target responses were unscored because generation ended
with `finish_reason = length`. A data-quality discrepancy was also
observed: `evaluation_status` was `pending` for rows even when `correct`
was populated. This discrepancy should be audited before formal
publication.

### Observed stratum counts

  ---------------------------------------------------------------------------------------
  Sample      Model      Stratum         Asked      Scored     Correct   Missing/unscored
                                                                 among 
                                                                scored 
  ----------- ---------- --------- ----------- ----------- ----------- ------------------
  Benchmark   DeepSeek   Short              62          62          59                  0

  Benchmark   DeepSeek   Long               38          38          38                  0

  Benchmark   GPT-OSS    Short              62          62          58                  0

  Benchmark   GPT-OSS    Long               38          38          38                  0

  Target      DeepSeek   Short              30          30          30                  0

  Target      DeepSeek   Long               70          67          67                  3

  Target      GPT-OSS    Short              30          30          29                  0

  Target      GPT-OSS    Long               70          70          68                  0
  ---------------------------------------------------------------------------------------

## 5. Stage 1 --- Pilot estimates and bootstrap intervals

Stage 1 used 10,000 stratified percentile-bootstrap replicates, seed
`20261025`, resampling outcomes within strata with fixed stratum sizes
and fixed weights.

### Accuracy estimates

  Model      Estimator                          Estimate
  ---------- -------------------------------- ----------
  DeepSeek   Benchmark-only                       0.9700
  DeepSeek   Reweighted benchmark                 0.9855
  DeepSeek   Direct target, scored outcomes       1.0000
  GPT-OSS    Benchmark-only                       0.9600
  GPT-OSS    Reweighted benchmark                 0.9806
  GPT-OSS    Direct target                        0.9700

### Bootstrap results

  ------------------------------------------------------------------------
  Model            Estimator                 Bootstrap SE 95% percentile
                                                          CI
  ---------------- ---------------- --------------------- ----------------
  DeepSeek         Benchmark-only                0.016852 \[0.930, 1.000\]

  DeepSeek         Reweighted                    0.008113 \[0.9661,
                   benchmark                              1.000\]

  DeepSeek         Direct target,                0.000000 \[1.000, 1.000\]
                   scored outcomes                        

  GPT-OSS          Benchmark-only                0.019249 \[0.920, 0.990\]

  GPT-OSS          Reweighted                    0.009487 \[0.9613,
                   benchmark                              0.9952\]

  GPT-OSS          Direct target                 0.017073 \[0.930, 1.000\]
  ------------------------------------------------------------------------

DeepSeek's degenerate direct-target bootstrap interval results from all
scored target responses being correct. It does not establish true
accuracy of 100% and does not account for the three unscored responses.

### Initial paired contrasts

DeepSeek minus GPT-OSS:

  Estimator                                  Difference 95% bootstrap CI
  ---------------------------------------- ------------ ---------------------
  Benchmark-only                                +0.0100 \[0.0000, 0.0300\]
  Reweighted benchmark                         +0.00484 \[0.0000, 0.01452\]
  Direct target, initial paired analysis        +0.0300 \[0.0000, 0.0700\]

These intervals touch zero. The results do not establish statistical
superiority or equivalence.

## 6. Stage 2 --- Monte Carlo estimator study

### Design

-   Seed: `20261026`.
-   Initial run: 2,000 Monte Carlo replicates per scenario cell; 500
    bootstrap replicates per sample.
-   Benchmark weights: Short 0.62, Long 0.38.
-   Target Short weights: 0.30 and 0.50.
-   Sample sizes: 100, 250, 500.
-   Benchmark correctness probabilities: Short 0.93, Long 0.97.

Target scenarios:

1.  **No conditional shift:** target Short 0.93, Long 0.97.
2.  **Uniform −5 percentage-point shift:** target Short 0.88, Long 0.92.
3.  **Differential shift:** target Short 0.85, Long 0.95.

Under no shift and target Short weight 0.30, target truth is
(0.30(0.93)+0.70(0.97)=0.958). Under target Short weight 0.50, it is
0.950.

### Main findings

-   Under conditional stability, reweighting reduces the composition
    bias of benchmark-only accuracy.
-   Direct target accuracy is approximately unbiased in the simulation
    because target observations are generated from the target
    distribution.
-   Under conditional shift, benchmark-based estimators can be biased
    for target accuracy even at large sample sizes.
-   A narrow interval around a biased estimator can have poor coverage
    of the target truth.

Illustrative initial-run results at target Short weight 0.30:

  ----------------------------------------------------------------------------------
  Scenario              n Estimator                   Bias         RMSE       95% CI
                                                                            coverage
  ---------- ------------ ---------------- --------------- ------------ ------------
  No shift            100 Benchmark-only          −0.01329      0.02617        0.941

  No shift            100 Reweighted             −0.000889      0.02190        0.826
                          benchmark                                     

  No shift            100 Direct target          +0.000145      0.02002       0.9195

  No shift            500 Reweighted             +0.000088      0.00958        0.937
                          benchmark                                     

  No shift            500 Direct target          +0.000170      0.00874        0.955

  Uniform −5          500 Benchmark-only          +0.03719          ---          ---
  pp shift                                                              

  Uniform −5          500 Reweighted              +0.05008          ---        0.014
  pp shift                benchmark                                     

  Uniform −5          500 Direct target    approximately 0          ---        0.954
  pp shift                                                              
  ----------------------------------------------------------------------------------

Dashes denote values not preserved in the recorded summary; they should
not be inferred.

## 7. Stages 2--3 and 7 --- Confidence-interval comparisons

### Methods compared

1.  Stratified percentile bootstrap.
2.  Wilson score intervals within strata combined using Bonferroni
    adjustment.
3.  Clopper--Pearson exact intervals within strata combined using
    Bonferroni adjustment.

The Wilson and Clopper--Pearson combinations are conservative reference
constructions for the weighted estimate, not automatically the preferred
methods in every setting. The weighted combination itself does not
guarantee nominal coverage under transport bias.

### Stage 2 CI validation

The initial implementation was optimized after a slow version was
interrupted. The optimized version used seed `20261027`, 3,000 Monte
Carlo replicates per cell, and 500 bootstrap replicates. A later reduced
user run used 1,000 Monte Carlo replicates and 200 bootstrap replicates.
Results differ slightly across runs because settings and random seeds
differ.

Reduced validation results, target Short weight 0.30:

  -------------------------------------------------------------------------------
  Scenario                     n Estimator              Bootstrap   CP-Bonferroni
                                                         coverage        coverage
  -------------- --------------- ---------------- --------------- ---------------
  No shift                   100 Benchmark-only             0.923           0.981

  No shift                   100 Reweighted                 0.825           0.998
                                 benchmark                        

  No shift                   100 Direct target              0.905           1.000

  No shift                   500 Benchmark-only             0.796           0.958

  No shift                   500 Reweighted                 0.939           1.000
                                 benchmark                        

  No shift                   500 Direct target              0.959           1.000

  Uniform −5 pp              500 Benchmark-only             0.075           0.550
  shift                                                           

  Uniform −5 pp              500 Reweighted                 0.009           0.206
  shift                          benchmark                        

  Uniform −5 pp              500 Direct target              0.944           1.000
  shift                                                           

  Differential               500 Benchmark-only             0.353           0.883
  shift                                                           

  Differential               500 Reweighted                 0.082           0.548
  shift                          benchmark                        

  Differential               500 Direct target              0.948           1.000
  shift                                                           
  -------------------------------------------------------------------------------

### Stage 3 CI-method comparison

Stage 3 used seed `20261028`, 1,000 Monte Carlo replicates and 200
bootstrap replicates. Under no shift at target Short weight 0.30 and
n=500:

  Estimator              Method                          Coverage   Mean interval width
  ---------------------- ----------------------------- ---------- ---------------------
  Benchmark-only         Percentile bootstrap               0.792                0.0384
  Benchmark-only         Wilson-Bonferroni                  0.910                0.0628
  Benchmark-only         Clopper--Pearson-Bonferroni        0.952                0.0650
  Reweighted benchmark   Percentile bootstrap               0.933                0.0364
  Reweighted benchmark   Wilson-Bonferroni                  0.996                0.0608
  Reweighted benchmark   Clopper--Pearson-Bonferroni        1.000                0.0627
  Direct target          Percentile bootstrap               0.953                0.0338
  Direct target          Wilson-Bonferroni                  1.000                0.0579
  Direct target          Clopper--Pearson-Bonferroni        1.000                0.0602

Under uniform −5 pp shift at n=500, Stage 3 bootstrap coverages were
0.103 for benchmark-only, 0.015 for reweighted benchmark, and 0.946 for
direct target. Under differential shift, they were 0.331, 0.078, and
0.953, respectively.

### Interpretation

-   Wilson and Clopper--Pearson constructions often increased coverage
    in tested scenarios, but were wider.
-   Conservative stratum intervals do not guarantee nominal coverage for
    a weighted estimate under systematic bias.
-   No confidence-interval method corrects transport bias.
-   With 1,000 Monte Carlo replicates, coverage near 0.95 has Monte
    Carlo standard error around 0.7 percentage points.
-   The evidence does not establish a universally best CI method.

## 8. Stage 4 --- Empirical intervals and missingness

Stage 4 compared bootstrap, Wilson-Bonferroni, and
Clopper--Pearson-Bonferroni intervals and calculated a
missing-as-incorrect sensitivity estimate.

  -------------------------------------------------------------------------------------------
  Model       Estimator        Estimate among Bootstrap   Wilson-Bonferroni   CP-Bonferroni
                                       scored 95% CI      95% CI              95% CI
                                     outcomes                                 
  ----------- ---------------- -------------- ----------- ------------------- ---------------
  DeepSeek    Benchmark-only           0.9700 \[0.9300,   \[0.8626, 0.9911\]  \[0.8664,
                                              1.0000\]                        0.9952\]

  DeepSeek    Reweighted               0.9855 \[0.9661,   \[0.8732, 0.9957\]  \[0.8791,
              benchmark                       1.0000\]                        0.9977\]

  DeepSeek    Direct target            1.0000 \[1.0000,   \[0.9081, 1.0000\]  \[0.9149,
                                              1.0000\]                        1.0000\]

  GPT-OSS     Benchmark-only           0.9600 \[0.9200,   \[0.8490, 0.9861\]  \[0.8523,
                                              0.9900\]                        0.9910\]

  GPT-OSS     Reweighted               0.9806 \[0.9613,   \[0.8667, 0.9933\]  \[0.8723,
              benchmark                       0.9952\]                        0.9957\]

  GPT-OSS     Direct target            0.9700 \[0.9300,   \[0.8618, 0.9938\]  \[0.8639,
                                              1.0000\]                        0.9982\]
  -------------------------------------------------------------------------------------------

For DeepSeek target Long questions, 70 responses were requested, 67 were
scored, and all 67 scored responses were correct. The scored-only
accuracy is 1.00. Counting all three unscored responses as incorrect
gives (67/70 `\approx 0.9571`{=tex}) in Long and target-weighted
accuracy (0.30(1.00)+0.70(67/70)=0.97).

The missing-as-incorrect value is a sensitivity rule, not proof that the
missing answers would be wrong. The bootstrap interval \[1,1\] is
degenerate because all scored outcomes were successes; it does not
establish true population accuracy of 100%.

## 9. Stage 5 --- Paired model comparison

Paired contrasts use shared question IDs and complete outcomes for both
models in the primary analysis.

  -----------------------------------------------------------------------------
  Estimator              DeepSeek minus 95%                         Total pairs
                                GPT-OSS paired-bootstrap   
                                        CI                 
  ---------------- -------------------- ------------------ --------------------
  Benchmark-only                +0.0100 \[0.0000, 0.0300\]                  100

  Reweighted                   +0.00484 \[0.0000,                           100
  benchmark                             0.01452\]          

  Direct target,               +0.03090 \[0.0000,                            97
  complete pairs                        0.07179\]          
  -----------------------------------------------------------------------------

Pairing audit:

-   Benchmark Short: 62 complete pairs; one DeepSeek-only correct
    outcome and no GPT-OSS-only correct outcomes.
-   Benchmark Long: 38 complete pairs; all both correct.
-   Target Short: 30 complete pairs; one DeepSeek-only correct outcome
    and no GPT-OSS-only correct outcomes.
-   Target Long: 67 complete pairs; three DeepSeek outcomes missing;
    among complete pairs, two DeepSeek-only correct, zero GPT-OSS-only
    correct, 65 both correct, and zero both wrong.

All observed complete-pair discordances favored DeepSeek, but
discordances were sparse and intervals included or touched zero. The
direct-target complete-pair estimate differs from Stage 1 because three
target Long pairs are excluded. When missing DeepSeek target outcomes
are counted as incorrect, the direct-target contrast is approximately
zero with an interval around \[−0.05, +0.05\] under the 30%/70% mixture.
Report the complete-pair and missing-as-incorrect results separately.

## 10. Stage 6 --- Target-composition sensitivity

Stage 6 varied target Short weight from 0 to 1 in increments of 0.05,
with Long weight as its complement. It used 5,000 bootstrap replicates
and seed `20261031`.

At target Short weight 0.30:

  ------------------------------------------------------------------------
  Analysis                    DeepSeek minus GPT-OSS 95% interval
  --------------------- ---------------------------- ---------------------
  Benchmark-only                             +0.0100 \[0.0000, 0.0300\]

  Reweighted benchmark                      +0.00484 \[0.0000, 0.01452\]

  Direct target,                            +0.03090 \[0.0000, 0.07135\]
  complete pairs                                     

  Direct target,                     approximately 0 \[−0.05, +0.05\]
  missing counted                                    
  incorrect                                          
  ------------------------------------------------------------------------

The reweighted contrast varies with the Short weight because the
observed benchmark contrast is concentrated in the Short stratum. The
complete-pair direct-target contrast remains positive and relatively
stable over tested mixtures because both observed stratum contrasts are
positive. The missing-as-incorrect contrast varies more with target
composition and can change sign. This is sensitivity analysis, not
evidence about the true target mix.

## 11. Stage 7 --- Cross-stage synthesis

At target Short weight 0.30 and n=500, representative coverage results
were:

  --------------------------------------------------------------------------
  Scenario        Estimator          Bootstrap coverage        CP-Bonferroni
                                                                    coverage
  --------------- ---------------- -------------------- --------------------
  No shift        Benchmark-only            about 0.796          about 0.958

  No shift        Reweighted                about 0.939                1.000
                  benchmark                             

  No shift        Direct target             about 0.959                1.000

  Uniform −5 pp   Benchmark-only            about 0.075          about 0.550
  shift                                                 

  Uniform −5 pp   Reweighted                about 0.009          about 0.206
  shift           benchmark                             

  Uniform −5 pp   Direct target             about 0.944                1.000
  shift                                                 

  Differential    Benchmark-only            about 0.353          about 0.883
  shift                                                 

  Differential    Reweighted                about 0.082          about 0.548
  shift           benchmark                             

  Differential    Direct target             about 0.948                1.000
  shift                                                 
  --------------------------------------------------------------------------

Coverage estimates vary slightly across runs because of differing
settings and seeds. The key distinction is between uncertainty around an
estimator and estimator bias for the target estimand. Wider intervals
may improve coverage in some designs but do not remove transport bias.

## 12. Stage 8 --- Empirical robustness consolidation

Stage 8 read five existing CSVs and produced 226 consolidated rows:

  Input                                         Rows
  ---------------------------------------- ---------
  Stage 4 empirical CI results                    18
  Stage 4 missingness audit                        8
  Stage 5 paired model comparison                  4
  Stage 5 pairing audit                            4
  Stage 6 target-composition sensitivity         192
  **Total**                                  **226**

Outputs: -
`analysis_outputs/gsm8k_stage8_empirical_robustness_synthesis.csv` -
`analysis_outputs/gsm8k_stage8_input_inventory.csv` -
`analysis_outputs/gsm8k_stage8_research_log.txt`

The consolidated file combines heterogeneous tables with source-stage
labels. Blank cells are expected where a source table did not contain a
column; they should not be interpreted as zero. Stage 8 consolidated the
records and did not recalculate the estimates.

## 13. Assumptions and limitations

### Empirical limitations

1.  The pilot is small and is not sufficient for a definitive model
    ranking.
2.  The target mix is designed, not verified against a real deployment
    population.
3.  Three DeepSeek target responses are missing; scored-only and
    missing-as-incorrect analyses answer different questions.
4.  Boundary estimates produce degenerate percentile-bootstrap
    intervals.
5.  Paired discordances are sparse.
6.  The `evaluation_status` field was `pending` even for some rows with
    correctness values and should be audited.
7.  The available data do not establish that the target sample
    represents a real deployment population.
8.  Results refer to the specific model identifiers used; they should
    not be generalized to later versions without evidence.

### Simulation limitations

1.  Only two strata and specified correctness probabilities were
    simulated.
2.  The simulated target distribution is known by construction, unlike
    most real target populations.
3.  Coverage estimates have Monte Carlo error.
4.  Reduced runs used 1,000 Monte Carlo and 200 bootstrap replicates, so
    small differences may reflect Monte Carlo variability.
5.  Simulations do not establish a universal estimator or interval
    ranking.
6.  Simulations are not empirical evidence of either LLM's real-world
    target accuracy.

## 14. Defensible conclusions and unsupported claims

### Supported by current evidence

-   Both models had high observed accuracy in this GSM8K pilot under the
    three specified estimators.
-   DeepSeek's 100% direct-target estimate applies to scored responses;
    three target Long responses were unscored.
-   Paired point estimates favored DeepSeek, but the intervals include
    or touch zero and do not establish superiority.
-   Reweighting substantially reduced composition bias under the
    simulated conditional-stability scenario.
-   Under simulated conditional shift, benchmark-based estimates were
    biased for target accuracy and intervals could severely undercover.
-   Direct target estimation was approximately unbiased under the
    simulated target-sampling design.
-   Confidence-interval construction affects coverage and width, but
    does not correct transport bias.
-   Results depend on target weights, target sample design, and
    missing-response handling.

### Not supported

-   DeepSeek is statistically superior to GPT-OSS.
-   DeepSeek's true target accuracy is 100%.
-   The 30%/70% mixture is the actual target-population composition.
-   Reweighting always corrects distribution shift.
-   Direct target estimation is universally best regardless of sampling
    design.
-   Wilson or Clopper--Pearson weighted intervals guarantee nominal
    coverage for every target estimand.
-   Simulation results prove real-world performance of either LLM.
-   The approach is validated for all benchmarks, target domains, or
    distribution-shift mechanisms.

## 15. Recommended reporting protocol

1.  Define the target estimand and mixture before presenting accuracy.
2.  Report benchmark-only, reweighted benchmark, and direct target
    estimates separately.
3.  Label the target mix as assumed unless independently measured.
4.  Report asked, scored, correct, and missing counts within each
    stratum.
5.  Present scored-only and missing-as-incorrect sensitivity results
    separately.
6.  Identify paired question IDs and missing-response handling in model
    contrasts.
7.  Report point estimates, interval method, interval bounds, and
    interval width.
8.  Separate empirical pilot results from simulation results.
9.  Record seeds, Monte Carlo replicates, and bootstrap replicates.
10. Do not select an interval method solely because it gives the
    narrowest interval.
11. Audit the `evaluation_status` discrepancy and preserve the raw CSV
    unchanged.
12. Do not rerun LLM inference unless new empirical data are explicitly
    approved.

## 16. Reproducibility record

The user ran the project in Python from VS Code on macOS in a project
virtual environment. The project folder is
`/Users/ayushmishra/LLM STATS AI`; derived files are stored in
`analysis_outputs/`.

  ----------------------------------------------------------------------------------------------
  Stage                   Script / output family                         Purpose
  ----------------------- ---------------------------------------------- -----------------------
  Pilot                   `scripts/05_gsm8k_pilot.py`;                   Existing two-model
                          `analysis_outputs/llm_two_model_results.csv`   pilot; do not rerun
                                                                         without approval

  Stage 1                 `06_gsm8k_stage1_bootstrap.py` and Stage 1     Data audit, estimates,
                          outputs                                        intervals, initial
                                                                         contrasts

  Stage 2                 `07_gsm8k_stage2_monte_carlo.py`               Estimator behavior
                                                                         under known simulated
                                                                         truth

  Stage 2 validation      `08_gsm8k_stage2_ci_validation.py`             CI coverage validation

  Stage 3                 `09_gsm8k_stage3_ci_method_comparison.py`      Bootstrap, Wilson and
                                                                         CP comparison

  Stage 4                 `10_gsm8k_stage4_empirical_ci.py`              Empirical intervals and
                                                                         missingness sensitivity

  Stage 5                 `11_gsm8k_stage5_paired_model_comparison.py`   Paired contrasts and
                                                                         pairing audit

  Stage 6                 `12_gsm8k_stage6_target_composition.py`        Target-mixture
                                                                         sensitivity

  Stage 7                 `13_gsm8k_stage7_ci_synthesis.py`              Cross-stage CI summary

  Stage 8                 `14_gsm8k_stage8_empirical_robustness.py`      Empirical robustness
                                                                         consolidation
  ----------------------------------------------------------------------------------------------

Recorded seeds: Stage 1 `20261025`; initial Stage 2 `20261026`;
optimized CI validation `20261027`; Stage 3 `20261028`; Stage 6
`20261031`. Distinguish planned and reduced validation settings in any
paper.

Retain the original pilot CSV unchanged, along with derived tables,
scripts, and stage logs.

## 17. Final conclusion

Target-population inference requires explicit attention to the estimand,
distribution shift, sampling design, missingness, and uncertainty
quantification. Benchmark-only accuracy describes the benchmark mixture.
Reweighted benchmark accuracy changes the weights and can correct
composition shift under conditional stability, but cannot by itself
correct conditional accuracy changes. Direct target estimation is
informative when the target sample is representative, but remains
vulnerable to sampling error, missing responses, and uncertain target
composition.

The GSM8K pilot illustrates these issues but is too small, with sparse
discordances and missing outcomes, to support a definitive model
ranking. The simulation study supplies controlled evidence about
estimator and interval behavior under the stated data-generating
mechanisms; it must not be presented as empirical LLM evidence.

The defensible contribution at this stage is a reproducible comparison
of the three agreed estimators, paired with explicit assumptions and
simulation evidence showing how conditional shift can undermine
benchmark-based inference. Further claims require a validated target
sampling design, an audit of the evaluation-status discrepancy, and
additional empirical evidence.

## Appendix A. Key output files to retain

-   `analysis_outputs/llm_two_model_results.csv`
-   `analysis_outputs/gsm8k_stage1_data_audit.csv`
-   `analysis_outputs/gsm8k_stage1_estimates.csv`
-   `analysis_outputs/gsm8k_stage1_model_differences.csv`
-   `analysis_outputs/gsm8k_stage2_monte_carlo_summary.csv`
-   `analysis_outputs/gsm8k_stage2_ci_validation_summary.csv`
-   `analysis_outputs/gsm8k_stage2_truth_and_bias_checks.csv`
-   `analysis_outputs/gsm8k_stage3_ci_method_comparison.csv`
-   `analysis_outputs/gsm8k_stage4_empirical_ci_results.csv`
-   `analysis_outputs/gsm8k_stage4_missingness_audit.csv`
-   `analysis_outputs/gsm8k_stage5_paired_model_comparison.csv`
-   `analysis_outputs/gsm8k_stage5_pairing_audit.csv`
-   `analysis_outputs/gsm8k_stage6_target_composition_sensitivity.csv`
-   `analysis_outputs/gsm8k_stage7_ci_synthesis_all.csv`
-   `analysis_outputs/gsm8k_stage7_ci_synthesis_q030.csv`
-   `analysis_outputs/gsm8k_stage8_empirical_robustness_synthesis.csv`
-   `analysis_outputs/gsm8k_stage8_input_inventory.csv`

## Appendix B. Final audit checklist

-   [ ] Three estimators remain unchanged.
-   [ ] Target mixture is described as assumed, not known truth.
-   [ ] Empirical and simulation evidence are separated.
-   [ ] Missing-response handling is explicit.
-   [ ] Paired comparisons disclose exclusions.
-   [ ] CI method and settings are reported.
-   [ ] Monte Carlo uncertainty is acknowledged.
-   [ ] `evaluation_status` discrepancy is investigated before
    publication.
-   [ ] Original pilot data are preserved.
-   [ ] No new LLM calls were made during Stages 2--9.



# Stage 10 Addendum — Whole-Sample Accuracy, Question-Length Distribution, and Exploratory Tests

**Stage status:** Completed using the existing pilot CSV. No new LLM calls were made and the original data were not modified.

## 18. Stage 10 objective

Stage 10 adds a descriptive whole-sample accuracy table, question-word-count distribution summaries, length-stratified correctness tests, and paired model comparisons by question-length group. It supplements the three agreed estimators; it does not introduce a fourth estimator.

Input: `analysis_outputs/llm_two_model_results.csv`  
Rows read: 400  
Unique questions: 200 (100 benchmark and 100 target questions)

## 19. Whole-sample accuracy

The table below reports accuracy among scored responses and Wilson 95% intervals.

| Sample | Model | Asked | Scored | Correct | Scored-only accuracy | Wilson 95% CI | Missing/unscored |
|---|---|---:|---:|---:|---:|---|---:|
| Benchmark | DeepSeek | 100 | 100 | 97 | 97% | [91.55%, 98.97%] | 0 |
| Benchmark | GPT-OSS | 100 | 100 | 96 | 96% | [90.16%, 98.43%] | 0 |
| Target | DeepSeek | 100 | 97 | 97 | 100% | [96.19%, 100%] | 3 |
| Target | GPT-OSS | 100 | 100 | 97 | 97% | [91.55%, 98.97%] | 0 |

**Interpretation:** DeepSeek's benchmark accuracy is 97% and GPT-OSS's is 96%, a descriptive difference of one percentage point. In the target sample, DeepSeek's scored-only accuracy is 100%, but three responses were unscored. Counting these as incorrect gives 97/100 = 97% across all requested target responses. The scored-only target estimate and the missing-as-incorrect sensitivity result must be reported separately.

The whole-benchmark accuracy is the ordinary overall accuracy for the benchmark's existing mixture of question lengths. It is the benchmark-only estimand, not a new estimator.

## 20. Question-length distribution

Each unique question was counted once for this distribution analysis.

| Statistic | Benchmark | Target |
|---|---:|---:|
| Questions with word counts | 100 | 100 |
| Mean word count | 42.44 | 52.24 |
| Median word count | 39 | 49 |

Tests comparing word counts between the observed benchmark and target sets:

| Test | Test statistic | p-value |
|---|---:|---:|
| Welch two-sample t-test | −3.81497 | 0.000183441 |
| Mann–Whitney U | 3254.5 | 0.0000199709 |
| Two-sample Kolmogorov–Smirnov | 0.33 | 0.0000321143 |

The three tests indicate that the observed question-word-count distributions differ. Welch's test targets a difference in means; Mann–Whitney compares ranks/distributions and is not strictly a test of medians in general; the KS test is sensitive to differences between the empirical distributions.

**Interpretation caveat:** the target question set was deliberately designed, including a 30% Short / 70% Long mix. These tests describe the difference between the two observed sets. They do not establish that a naturally occurring deployment population has this distribution, nor do they establish the cause of the difference.

## 21. Correctness by question-length group

Fisher's exact test was used as the primary exploratory association test because correctness is binary and some cells have zero errors. One-way ANOVA on the 0/1 correctness indicator was also computed as a diagnostic only.

| Sample | Model | Short correct/scored | Long correct/scored | Fisher exact p-value | Exploratory ANOVA F | Exploratory ANOVA p-value |
|---|---|---:|---:|---:|---:|---:|
| Benchmark | DeepSeek | 59/62 | 38/38 | 0.286061 | 1.893559 | 0.171937 |
| Benchmark | GPT-OSS | 58/62 | 38/38 | 0.294472 | 2.568276 | 0.112245 |
| Target | DeepSeek | 30/30 | 67/67 | 1.000000 | Not defined | Not defined |
| Target | GPT-OSS | 29/30 | 68/70 | 1.000000 | 0.016039 | 0.899480 |

The target DeepSeek ANOVA was not defined in the script output because all scored correctness outcomes were 1, leaving no within-group variation. The target Fisher test returned p=1 with no observed errors in either group; the odds ratio is undefined.

At the conventional 5% level, none of the reported Fisher tests or defined exploratory ANOVA tests rejects equal correctness rates across Short and Long groups. This is not evidence that the true rates are identical. The sample is small, observed accuracy is close to 100%, and the number of errors is very low, so the tests have limited ability to detect differences.

**Methodological decision:** do not make ordinary ANOVA the primary inferential test for this binary outcome. Retain it as an exploratory diagnostic if useful, and report the binary counts and Fisher's exact results. The existing estimator and transport analysis remains the central analysis.

## 22. Paired model comparisons by length group

The exact paired discordance test was calculated on question IDs with scored outcomes for both models. It is an exact binomial test on the discordant pairs, equivalent to the exact form of McNemar's test for a 2×2 paired binary outcome.

| Sample | Length | Complete pairs | DeepSeek-only correct | GPT-OSS-only correct | Discordant pairs | Exact paired p-value | Paired accuracy difference |
|---|---|---:|---:|---:|---:|---:|---:|
| Benchmark | Long | 38 | 0 | 0 | 0 | 1.000 | 0.0000 |
| Benchmark | Short | 62 | 1 | 0 | 1 | 1.000 | +0.01613 |
| Target | Long | 67 | 2 | 0 | 2 | 0.500 | +0.02985 |
| Target | Short | 30 | 1 | 0 | 1 | 1.000 | +0.03333 |

All observed discordant complete pairs favored DeepSeek, but there were only one or two discordances in each nonempty group. The exact tests therefore do not provide evidence of a statistically significant model difference. A p-value of 1.0 with only one discordant pair does not mean the models are equivalent; it reflects the extremely limited information in that stratum.

These length-specific results are consistent with the earlier whole-benchmark paired analysis: the observed benchmark accuracy difference is +1 percentage point, but the paired bootstrap interval touched zero. Missing DeepSeek target outcomes mean that target complete-pair results and missing-as-incorrect sensitivity results answer different questions.

## 23. Stage 10 files and reproducibility

Stage 10 created:

- `analysis_outputs/gsm8k_stage10_accuracy_metrics_by_length.csv`
- `analysis_outputs/gsm8k_stage10_question_length_distribution.csv`
- `analysis_outputs/gsm8k_stage10_length_accuracy_tests.csv`
- `analysis_outputs/gsm8k_stage10_paired_model_tests_by_length.csv`
- `analysis_outputs/gsm8k_stage10_length_metrics_research_log.txt`

Script: `scripts/15_gsm8k_stage10_length_metrics_anova.py`

The script read the existing pilot CSV only. No LLM calls were made and no simulations were performed.

## 24. Integrated interpretation after Stage 10

Stage 10 adds three useful empirical facts:

1. **Whole-benchmark baseline:** DeepSeek scored 97/100 and GPT-OSS scored 96/100. The observed difference is small and does not establish statistical superiority.
2. **Question-length distribution:** the target set has a higher mean word count than the benchmark set (52.24 versus 42.44), and the observed word-count distributions differ in the reported tests. The designed Short/Long mixture and question-word-count differences should be described transparently; neither by itself establishes a real deployment population shift.
3. **Length and correctness:** the pilot does not provide statistically significant evidence of a Short-versus-Long accuracy difference in the reported exploratory tests. Near-perfect accuracy and sparse errors limit power, so the result should not be interpreted as proof that question length has no relationship with correctness.

These descriptive and exploratory analyses complement the main transport-inference study. The central methodological conclusion remains unchanged: reweighting can address composition shift under conditional stability, but cannot correct conditional accuracy shift on its own. Direct target inference depends on a representative target sample, valid target weights, and explicit handling of missing responses.

## 25. Updated final audit checklist

- [ ] Three agreed estimators remain unchanged.
- [ ] Whole-benchmark accuracy is clearly labelled as benchmark-only accuracy.
- [ ] Scored-only target accuracy and missing-as-incorrect sensitivity are separated.
- [ ] Target mixture is described as designed/assumed rather than known population truth.
- [ ] Empirical results are separated from simulation results.
- [ ] Question word-count summaries and their tests are reported as observed-set comparisons.
- [ ] Binary correctness is not analyzed with ANOVA as the sole or primary inferential test.
- [ ] Paired model comparisons report discordant counts and missing-pair handling.
- [ ] The `evaluation_status` discrepancy is investigated before publication.
- [ ] Original pilot data are preserved.
- [ ] No new LLM calls were made in Stages 2–10.

3. Models and evaluation output
Model identifiers appearing in the evaluation results:
- deepseek-ai/DeepSeek-V4-Pro-0813
- openai/gpt-oss-120b
Evaluation file:

analysis_outputs/gsm8k_5000/gsm8k_5000_model_results.csv

The reported columns include sample type, question ID, source split, question text and word count, length group, model, prompt version, gold and predicted answers, correctness, scoring/evaluation status, model response, finish reason, token counts, error, and timestamps.

Other files reported in analysis_outputs/gsm8k_5000/:
- gsm8k_5000_B.csv
- gsm8k_5000_T.csv
- gsm8k_5000_model_results.csv
- gsm8k_5000_run_metadata.json
- gsm8k_5000_sample_audit.json
- gsm8k_5000_sample_research_log.txt
- 09_estimator_bootstrap_monte_carlo.py
The evaluation CSV was reported to contain 10,251 rows. The intended design would normally imply 10,000 model-question evaluations if both models were evaluated once on each of the 5,000 selected questions (5,000 × 2). The additional 251 rows must be reconciled before final analysis; possible explanations must be checked from the data rather than assumed.
4. Preliminary accuracy results
These are scored-response-only descriptive estimates from the current CSV summary. They are not yet final inferential results.
Sample	Model	Rows in file	Scored	Correct	Accuracy among scored
Benchmark	DeepSeek V4 Pro	4,197	3,863	3,785	97.9808%
Benchmark	GPT-OSS 120B	4,008	3,998	3,845	96.1731%
Target	DeepSeek V4 Pro	1,043	967	951	98.3454%
Target	GPT-OSS 120B	1,003	999	961	96.1962%


Accuracy by question-length group
Sample	Model	Length group	Rows	Scored	Correct	Accuracy among scored
Benchmark	DeepSeek V4 Pro	Long	1,845	1,630	1,593	97.7301%
Benchmark	DeepSeek V4 Pro	Short	2,352	2,233	2,192	98.1639%
Benchmark	GPT-OSS 120B	Long	1,725	1,718	1,640	95.4598%
Benchmark	GPT-OSS 120B	Short	2,283	2,280	2,205	96.7105%
Target	DeepSeek V4 Pro	Long	460	405	394	97.2840%
Target	DeepSeek V4 Pro	Short	583	562	557	99.1103%
Target	GPT-OSS 120B	Long	433	429	405	94.4056%
Target	GPT-OSS 120B	Short	570	570	556	97.5439%


Preliminary descriptive interpretation
- Among responses that were scored, DeepSeek's observed accuracy is higher than GPT-OSS's in both the benchmark and target samples.
- Both models have lower observed scored-only accuracy on long than short questions in the reported length breakdown.
- These patterns are descriptive only. They do not establish statistically significant model superiority or a causal effect of question length.
- The scored-only denominator may introduce bias if scoring failures are related to model, question difficulty, response format, or correctness. The missingness mechanism must be investigated.
## 1. 

This log starts at Stage 12 and continues through the full-scale Direction B pairwise-preference evaluation and its offline inference. The project has two distinct analyses:

- **Direction A — correctness-based inference:** model accuracy on benchmark and target samples, including paired comparisons, standardization, and sensitivity to unscored responses.
- **Direction B — judge-based pairwise preference:** an independent LLM judge compares two models' answers to the same question and chooses A, B, or a tie.

These outcomes must not be conflated. Judge preference is not gold-answer accuracy, and neither is by itself a universal measure of model quality.

## 2. Models and data

- DeepSeek: `deepseek-ai/DeepSeek-V4-Pro-0813`
- GPT-OSS: `openai/gpt-oss-120b`
- Direction B judge: `zai-org/GLM-5.3`
- Dataset: GSM8K mathematical word problems
- Direction B prompt version: `direction_b_gsm8k_v1`

The judge rubric covered correctness, validity of reasoning, completeness, and clarity; verbosity was not rewarded; ties were allowed. Outputs included preference, correctness labels, reasoning labels, confidence, and rationale. The two candidate answers were matched on question ID. Answer identities were randomized as A/B during judging and decoded afterward with a separate private key.

Existing outputs were reused for offline analysis. The inference script did not call any LLM API. Private mapping files should not be published with public research artifacts.

# Part I — Stage 12 corrected correctness analysis

## 3. Objective

Stage 12 continued the correctness-based comparison of DeepSeek and GPT-OSS using matched GSM8K questions. The corrected analysis examined benchmark and target differences, standardized differences, and sensitivity to assumptions about unscored responses.

## 4. Stage 12 outputs

Output directory: `analysis_outputs/gsm8k_pairwise/stage12_corrected/`

| Reference comparison | DeepSeek − GPT-OSS |
|---|---:|
| Benchmark | +2.127 percentage points |
| Target | +2.131 percentage points |
| Standardized result 1 | +2.391 percentage points |
| Standardized result 2 | +2.378 percentage points |

These are reference point differences from the corrected Stage 12 analysis. They do not, by themselves, establish statistical significance.

### Sensitivity analysis

- Zero was included in 12 of 36 expanded sensitivity scenarios.
- Worst-case bounds included zero.

### Interpretation

The reference point estimates were positive and similar across the reported comparisons, favoring DeepSeek descriptively. However, the sensitivity results show that the conclusion is not invariant to every scenario or assumption. The defensible interpretation is that the reference estimates favor DeepSeek while robustness depends on assumptions. Exact scenario-level estimates and interval endpoints are not reproduced here because they are not available in the consolidated result excerpt; no missing numbers are inferred.

# Part II — Direction B pairwise preference evaluation

## 5. Research questions

1. How often does the judge prefer either answer versus returning a tie?
2. Among decisive judgments, is the winner split balanced?
3. When exactly one model answer is gold-correct, does the judge prefer the correct answer?
4. Are mapped outcomes stable when answer order is reversed?

These questions concern the judge's behavior under a specified protocol, not a direct estimate of model accuracy.

## 6. Full-scale sample design

The sample contains 1,500 unique matched questions, selected with seed `20261031`:

| Sample type | Questions |
|---|---:|
| Benchmark | 1,200 |
| Target | 300 |
| **Total** | **1,500** |

Primary input artifacts in `analysis_outputs/gsm8k_direction_b/`:

- `direction_b_1500_pair_sample.csv` — question, gold answer, model responses, and saved correctness indicators
- `judge_input_1500.csv` — judge-facing input
- `judge_key_private_1500.csv` — private A/B identity and correctness mapping
- `judge_results_glm53_1500.jsonl` — saved judge outputs
- `judge_run_log_glm53_1500.json` — run metadata and usage log

## 7. Full-scale judging run

Runner: `scripts/24_gsm8k_direction_b_run_judge_1500.py`

Reported run summary:

- 1,500 primary judgments
- 150 reversed-order judgments
- 1,650 unique run IDs
- Zero non-success records
- Input tokens: 902,245
- Output tokens: 122,437
- Total tokens: 1,024,682
- Estimated cost: USD 1.8019

These are the completed run's reported counts. Later inference used saved results and did not make API calls.

## 8. Mapping and audit

Mapping script: `scripts/25_gsm8k_direction_b_analyze_judgments.py`

Outputs:

- `direction_b_primary_judgments_mapped.csv`
- `direction_b_reversed_order_mapped.csv`
- `direction_b_order_consistency.csv`

The initial mapping attempt did not extract a nested judge result because the parser did not inspect the `judge_result` key. The parser was patched to inspect that key, after which the saved judgments mapped successfully. This was a parsing fix, not a new judging run.

# Part III — Primary preference results

## 9. Overall results

| Primary judge outcome | Count | Percentage |
|---|---:|---:|
| DeepSeek preferred | 64 | 4.27% |
| GPT-OSS preferred | 14 | 0.93% |
| Tie | 1,422 | 94.80% |
| **Total** | **1,500** | **100%** |

There were 78 decisive judgments (5.20% of all comparisons). Among decisive judgments, DeepSeek received 64/78 = **82.05%** and GPT-OSS 14/78 = **17.95%**.

The 82.05% is conditional on the judge making a non-tie choice; it is not DeepSeek's share of all comparisons. The tie rate is a central result, not a nuisance to discard.

## 10. Benchmark and target breakdown

### Benchmark (n = 1,200)

| Outcome | Count | Percentage |
|---|---:|---:|
| DeepSeek preferred | 57 | 4.75% |
| GPT-OSS preferred | 12 | 1.00% |
| Tie | 1,131 | 94.25% |

There were 69 decisive benchmark judgments. DeepSeek received 57/69 = **82.61%** of those decisive preferences.

### Target (n = 300)

| Outcome | Count | Percentage |
|---|---:|---:|
| DeepSeek preferred | 7 | 2.33% |
| GPT-OSS preferred | 2 | 0.67% |
| Tie | 291 | 97.00% |

There were only 9 decisive target judgments. DeepSeek received 7/9 = **77.78%** of those decisive preferences.

### Interpretation

Decisive preferences descriptively favored DeepSeek in both samples, but the target estimate is based on only nine decisive judgments and is highly uncertain. These observed proportions do not establish a benchmark-versus-target difference or identify a causal effect of distribution shift.

## 11. Exact binomial test

The offline inference script tested the 78 decisive judgments against a 50:50 winner split:

- Decisive judgments: 78
- DeepSeek wins: 64
- GPT-OSS wins: 14
- Exact two-sided p-value: **8.58218 × 10⁻⁹**

Under the binomial model and its assumptions, the decisive winner split is inconsistent with 50:50. The test is conditional on the judge not returning a tie. It does not test equality of model accuracy and does not validate the judge as an unbiased measure of human preference. For reporting, show the p-value together with the 94.8% tie rate and absolute counts.

# Part IV — Gold correctness and alignment

## 12. Correctness parsing issue and fix

The first run of `scripts/26_gsm8k_direction_b_inference.py` reported all 1,500 gold-correctness labels missing. Inspection showed that saved correctness values were numeric representations such as `1.0` and `0.0`, but the parser only recognized string values `"1"` and `"0"`.

The boolean conversion function was corrected to recognize `true`, `1`, `1.0`, `yes`, `correct` as true, and `false`, `0`, `0.0`, `no`, `incorrect` as false. The corrected run parsed the saved labels. This changed parsing only; it did not change model responses, gold answers, judge outputs, or sample data.

## 13. Gold-correctness states

| Gold correctness status | Pairs | Percentage |
|---|---:|---:|
| Both correct | 1,472 | 98.13% |
| Neither correct | 10 | 0.67% |
| DeepSeek only correct | 12 | 0.80% |
| GPT-OSS only correct | 6 | 0.40% |
| **Total** | **1,500** | **100%** |

Exactly one answer was correct in 18 pairs (1.2%). This small discordant subset limits the precision of correctness-alignment analysis. These are saved gold-correctness indicators, not the judge's own correctness labels.

## 14. Judge alignment among exactly-one-correct pairs

| Judge outcome relative to gold correctness | Count | Percentage of 18 |
|---|---:|---:|
| Aligned — preferred the correct answer | 15 | 83.33% |
| Tie | 1 | 5.56% |
| Misaligned — preferred the incorrect answer | 2 | 11.11% |
| **Total** | **18** | **100%** |

The judge preferred the gold-correct answer in 15 of 18 discordant pairs. This is a descriptive, small-sample diagnostic and is not enough to claim strong general judge reliability. A few additional cases could materially change the proportions. It says nothing definitive about judge behavior when both answers are correct or both are incorrect.

# Part V — Reversed-order robustness

## 15. Results

A subset of 150 pairs was judged again with answer order reversed:

- Reversed-order pairs: 150
- Same mapped outcome, including ties: **148/150 = 98.67%**
- Decisive in both orderings: 4 pairs
- Same winner among those four: 4/4

### Interpretation

Mapped outcomes were stable for 148 of 150 repeated pairs, but agreement was dominated by ties. Only four pairs were decisive in both orderings, so the evidence about order robustness when a winner is selected is very limited. This does not justify claiming that order bias is absent. Repeated judgments should not be added to the primary n = 1,500 denominator because they are repeated evaluations of a subset, not new questions.

# Part VI — Scripts, artifacts, and reproducibility

## 16. Inference script

Script: `scripts/26_gsm8k_direction_b_inference.py`

Inputs:

- `analysis_outputs/gsm8k_direction_b/judge_results_glm53_1500.jsonl`
- `analysis_outputs/gsm8k_direction_b/judge_key_private_1500.csv`
- `analysis_outputs/gsm8k_direction_b/direction_b_1500_pair_sample.csv`

Functions performed by the script:

- Map primary judge outputs to question IDs and model identities
- Summarize overall, benchmark, and target preference outcomes
- Calculate Wilson confidence intervals for proportions
- Perform an exact two-sided binomial test of decisive wins
- Compare benchmark and target decisive-winner composition using Fisher's exact test
- Compare benchmark and target tie/decisive rates using Fisher's exact test
- Assess preference alignment for pairs where exactly one answer is gold-correct
- Save a summary CSV and Markdown research log

Outputs:

- `analysis_outputs/gsm8k_direction_b/direction_b_inference_summary.csv`
- `analysis_outputs/gsm8k_direction_b/gsm8k_direction_b_research_log.md`

The terminal output confirmed 1,500 mapped primary judgments, 78 decisive outcomes (64 DeepSeek, 14 GPT-OSS), the exact p-value, the four correctness-state counts, and alignment counts of 15/1/2. It confirmed both files were saved and that no API calls were made. Exact Wilson interval and Fisher-test values were not included in the shared terminal excerpt; retrieve them from the saved summary CSV rather than guessing.

## 17. Reproduction commands

From the project root:

```bash
python scripts/25_gsm8k_direction_b_analyze_judgments.py
python scripts/26_gsm8k_direction_b_inference.py
```

These commands map and analyze saved results offline. Do not rerun `scripts/24_gsm8k_direction_b_run_judge_1500.py` just to regenerate summary statistics; the saved JSONL and CSV files suffice.

## 18. Artifact ledger

| Component | Artifact | Purpose |
|---|---|---|
| Stage 12 | `analysis_outputs/gsm8k_pairwise/stage12_corrected/` | Corrected correctness and sensitivity outputs |
| Direction B sample | `direction_b_1500_pair_sample.csv` | Paired responses and gold correctness |
| Judge input | `judge_input_1500.csv` | Blinded judge inputs |
| Private mapping | `judge_key_private_1500.csv` | A/B model identity and correctness key |
| Raw judge results | `judge_results_glm53_1500.jsonl` | Saved primary/reversed judgments |
| Run log | `judge_run_log_glm53_1500.json` | Usage and run metadata |
| Mapped primary | `direction_b_primary_judgments_mapped.csv` | Decoded primary outcomes |
| Mapped reversed | `direction_b_reversed_order_mapped.csv` | Decoded reversed-order outcomes |
| Order consistency | `direction_b_order_consistency.csv` | Reversal agreement |
| Inference summary | `direction_b_inference_summary.csv` | CI and test summaries |
| Generated log | `gsm8k_direction_b_research_log.md` | Direction B component log |

Direction B files are in `analysis_outputs/gsm8k_direction_b/`, except the Stage 12 folder.

# Part VII — Assumptions and interpretation boundaries

## 19. Estimands

The main Direction B estimands are:

1. Overall proportions of DeepSeek wins, GPT-OSS wins, and ties under this pair-sampling and judge protocol.
2. Conditional probability that DeepSeek is preferred given a decisive judgment.
3. The same outcome proportions within benchmark and target samples.
4. Correctness alignment among pairs where exactly one answer is gold-correct.
5. Mapped-outcome agreement between original and reversed answer ordering in the 150 repeated pairs.

These are sample- and protocol-specific, not universal model properties.

## 20. Assumptions and limitations

- Inference beyond the sampled questions depends on the sample design and target population definition.
- Model outcomes are paired because both models answer the same question; analyses should preserve this dependence.
- GLM-5.3 is an automated judge, not a human gold-standard preference panel.
- The 94.8% tie rate leaves only 78 decisive comparisons for the overall conditional-winner test.
- The target subset has only nine decisive judgments.
- Correctness alignment has only 18 discordant-correctness pairs.
- Reversal agreement is tie-dominated; only four pairs were decisive in both orderings.
- Benchmark/target differences do not, by themselves, identify a causal effect of distribution shift.
- Judge preference must not be substituted for gold-answer accuracy.
- Secondary benchmark-target tests should be called exploratory unless a testing plan was specified in advance.
- Direction A missing/unscored responses are a separate issue and must be handled in the correctness analysis rather than combined silently with Direction B outcomes.

## 21. Relationship to the wider 5,000-question correctness evaluation

The wider GSM8K correctness project used 4,000 benchmark questions and 1,000 target questions, sampled without overlap. Previously reported scored/unscored counts were:

| Sample/model | Scored | Unscored |
|---|---:|---:|
| Benchmark DeepSeek | 3,803 | 197 |
| Benchmark GPT-OSS | 3,992 | 8 |
| Target DeepSeek | 957 | 43 |
| Target GPT-OSS | 997 | 3 |

Previously reported accuracy among scored responses:

| Sample/model | Accuracy among scored |
|---|---:|
| Benchmark DeepSeek | 0.983697 |
| Benchmark GPT-OSS | 0.962425 |
| Target DeepSeek | 0.985371 |
| Target GPT-OSS | 0.963892 |

These belong to Direction A, not the 1,500-pair Direction B analysis. Keep them in separate tables and explicitly state that the accuracy figures are among scored answers; missingness handling remains important.

# Part VIII — Conclusions and next steps

## 22. Conclusions supported by the shared results

1. Stage 12 corrected reference differences were positive: +2.127 pp benchmark, +2.131 pp target, and +2.391/+2.378 pp for the two standardized results. Sensitivity cautions against claiming robustness across all assumptions: zero was included in 12/36 expanded scenarios and in worst-case bounds.
2. Direction B full-scale judging completed with 1,500 primary pairs and 150 reversed-order repeats; the run log reported zero non-success records.
3. Ties dominate: 1,422/1,500 primary comparisons (94.8%).
4. Among 78 decisive judgments, DeepSeek won 64 (82.05%); the exact two-sided binomial p-value was 8.58218 × 10⁻⁹ under the specified 50:50 decisive-win null.
5. The target decisive comparison is weakly informed: 7 DeepSeek wins versus 2 GPT-OSS wins among nine decisive outcomes.
6. Only 18 pairs had exactly one gold-correct answer; the judge aligned in 15, tied in one, and misaligned in two.
7. Reversed-order outcomes agreed in 148/150 cases, but only four were decisive in both orderings.
8. Pairwise preference and gold-answer accuracy are distinct outcomes and must remain separate in interpretation.

## 23. Next steps

1. Preserve raw judge outputs and the private key unchanged.
2. Inspect `direction_b_inference_summary.csv` for exact Wilson intervals and Fisher-test results; report those only after verifying their denominators and interpretation.
3. Report unconditional preference proportions alongside conditional decisive-winner results.
4. Treat correctness alignment as exploratory because n = 18.
5. Treat order robustness as preliminary because only four pairs were decisive in both orderings.
6. Present Stage 12 correctness inference and Direction B judge preference as distinct analyses in the paper.
7. Do not make additional judge API calls solely to reproduce offline statistics.

## Appendix — Compact result ledger

| Item | Result |
|---|---|
| Stage 12 output directory | `analysis_outputs/gsm8k_pairwise/stage12_corrected/` |
| Stage 12 differences | +2.127 pp; +2.131 pp; +2.391 pp; +2.378 pp |
| Sensitivity | Zero in 12/36 scenarios; worst-case bounds include zero |
| Direction B sample | 1,500 = 1,200 benchmark + 300 target |
| Judge/prompt | `zai-org/GLM-5.3`; `direction_b_gsm8k_v1` |
| Primary results | 64 DeepSeek / 14 GPT-OSS / 1,422 ties |
| Decisive split | 82.05% DeepSeek / 17.95% GPT-OSS |
| Exact binomial p-value | `8.58218e-09` |
| Gold correctness | 1,472 both; 10 neither; 12 DeepSeek only; 6 GPT-OSS only |
| Correctness alignment | 15 aligned / 1 tie / 2 misaligned, n = 18 |
| Order reversal | 148/150 same outcome; 4/4 same winner among decisive in both |
| API calls in inference script | None |

