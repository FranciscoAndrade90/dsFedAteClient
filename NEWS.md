# dsRctFedClient 0.0.0.9000

## New features

* `ds.EstimateAte()` now supports variance estimation for binary outcomes
  (`outcome_type = "binary"`, `compute_variance = TRUE`). The delta-method
  variance has two components:
  - *Parameter estimation*: **g**₁ᵀ **H**₁⁻¹ **g**₁ + **g**₀ᵀ **H**₀⁻¹ **g**₀,
    where **H**_w⁻¹ = `ds.glm(viewVarCov = TRUE)$VarCovMatrix` (the per-arm
    pooled variance–covariance matrix, equal to H_w⁻¹ for binomial dispersion = 1)
    and **g**_w is the gradient of τ̂ w.r.t. θ_w, assembled from unnormalised
    gradient sums returned by `localAteDS`.
  - *Sampling*: Var(δᵢ)/n where δᵢ = μ̂₁(xᵢ) − μ̂₀(xᵢ), pooled from per-site
    scatter scalars via the parallel-axis theorem.
  The per-arm Fisher information matrices are obtained from `ds.glm` rather
  than assembled from per-site disclosures, minimising server-side output.
  Verified at machine precision (ATE diff < 1e-10, variance diff < 1e-10)
  against the manual pooled formula using the same federated coefficients across
  K = 3 heterogeneous sites. The `heterogeneous` argument has no effect for
  binary outcomes — the formula is automatically valid under covariate shift.

* `ds.EstimateAte()` gains a `heterogeneous` argument (default `TRUE`)
  implementing the leverage-corrected asymptotic variance under full
  heterogeneity (covariate shift `H → X` *and* heterogeneous treatment
  assignment `H → W`).  When `TRUE` and `compute_variance = TRUE`:
  - Both `ds.glm()` calls are made with `viewVarCov = TRUE` to obtain the
    pooled per-arm variance–covariance matrices `VarCovMatrix^(w)`.
  - Each server returns its covariate column sums `n_k * X̄_k` via the new
    `heterogeneous` argument of `localAteDS()`.  The client sums these to form
    the global covariate mean `X̄`, then augments it to `X̄' = (1, X̄ᵀ)ᵀ`.
  - Term (I) of the variance is computed via eq. (15) of the companion variance
    paper rather than the homogeneous closed form:
    `Î = σ̂² * X̄'ᵀ (VarCovMatrix¹ + VarCovMatrix⁰) X̄'`,
    where `VarCovMatrix^(w)` is the unscaled Gram inverse `(G^(w))⁻¹` returned
    by `ds.glm(viewVarCov = TRUE)`.  This automatically
    carries the leverage correction `δ^(w) ≥ 0` (Lemma 2 of the companion
    variance paper) and is algebraically identical to the centralised plug-in
    variance (Proposition 2).
  - Term (II) and `σ̂²` are unchanged.
  - When `compute_variance = FALSE`, the argument has no effect.
  - Works uniformly when combined with `study_effects = TRUE`: the augmented
    covariate mean `X̄ = (X̄_X, ρ̂)` (including the study-proportion D-block)
    is obtained from the same `n_k_X_bar` scalars without any additional
    communication round.

---

## Breaking changes

* `ds.EstimateAte()` now returns a named list instead of a single numeric
  scalar. The list contains: `ate`, `variance`, `ci_lower`, `ci_upper`,
  `n_total`, `n_k`, `ate_k`. Any code that treated the return value as a
  number must be updated to use `result$ate`.

## New features

* `ds.EstimateAte()` gains a `study_effects` argument (default `FALSE`)
  implementing the adjusted GD estimator (Definition 9 of Khellaf, Bellet &
  Josse 2025) for multi-site settings with confounded study effects:
  - When `TRUE` and K > 1, a call to the new `addStudyIdDS()` assign function
    creates a new server-side object `<base_name>_study_adj` containing the
    original data frame augmented with K-1 integer study-dummy columns. The
    caller's original object is never modified. The federated GLMs are then
    fitted on the augmented frame with the study dummies in the formula, removing
    confounding between heterogeneous treatment probabilities and study baseline
    shifts without requiring a separate communication round.
  - The degrees of freedom in the variance estimator adjust automatically
    from n - 2(d+1) to n - 2(d+K) because `covariate_names` is extended with
    the K-1 study dummies before being passed to the rest of the pipeline.
  - When `study_effects = TRUE` but K = 1, a warning is emitted and the
    unadjusted estimator is used as a fallback.
  - Simulation shows that with h = c(3, 0, -3) and p = c(0.8, 0.5, 0.2), the
    unadjusted ATE is biased by ≈ 2.3 while the adjusted estimator has bias
    < 0.1 and its 95% CI covers the true ATE of 0.5.

* `ds.EstimateAte()` now computes and returns the full delta-method variance
  of the GD estimator. Two components are estimated:
  - OLS noise: sigma-hat-squared / (n * p-hat * (1 - p-hat)), using pooled
    df n - 2(d+1) from per-site residual sums of squares.
  - Covariate-mean uncertainty: Delta-beta-hat^T Sigma-hat Delta-beta-hat / n,
    pooled from per-site projected scatter scalars via the parallel-axis
    decomposition.
  Both components are assembled from a single call to `localAteDS()` with no
  additional communication round.
* 95% confidence intervals (`ci_lower`, `ci_upper`) are included in the return
  value.
* `.check_sample_sizes()` warns when either federated arm size falls below
  5*(d+1), signalling that the asymptotic approximation may be unreliable.

## Documentation

* Vignette revised: notation unified throughout (`alpha_w` / `beta_w` for
  true parameters, `alpha_hat_w` / `beta_hat_w` for estimates, matching the
  code); variance section rewritten with pointer to Khellaf et al. Table 1,
  plug-in formula, and RSS_k definition; Disclosure paragraph added explaining
  the two-step information flow (IRLS matrices + per-server scalars).
* Example 2 added to vignette illustrating how extreme treatment probabilities
  trigger `nfilter.glm` and how permissive mode enables estimation.

## Bug fixes

* `dev/SettingC_1.R`: corrected argument name from `connections` to
  `datasources`.

## Earlier changes

* Renamed `connections` argument to `datasources` in `ds.EstimateAte()`.
* Added `datashield.connections_find()` fallback when `datasources = NULL`.
* Added `.check_datasources()`, `.check_args()`, `.check_string()` for input
  validation.
* Renamed internal helper to `.localAte()` and marked it internal.
* Added `cli` and `methods` under `Imports` in `DESCRIPTION`.
* Added unit tests for `ds.EstimateAte()` and all input validation helpers.
* `tests/testthat/setupDSLite.R` now uses `devtools::load_all()` instead of
  `library(dsRctFed)` so client tests always run against current source
  without requiring a prior install.
