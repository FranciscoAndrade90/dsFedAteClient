# dsRctFedClient 0.0.0.9000

## Breaking changes

* `ds.EstimateAte()` now returns a named list instead of a single numeric
  scalar. The list contains: `ate`, `variance`, `ci_lower`, `ci_upper`,
  `n_total`, `n_k`, `ate_k`. Any code that treated the return value as a
  number must be updated to use `result$ate`.

## New features

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
