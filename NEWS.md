# dsFedAteClient 0.1.0

Initial CRAN release.

## Features

* `ds.EstimateAte()` orchestrates federated ATE estimation across K
  independent studies using arm-specific OLS models fitted via `ds.glm()`.

* `compute_variance` argument (default `TRUE`) returns a delta-method
  variance estimate and 95% confidence interval assembled from per-server
  residual sums of squares and projected scatter scalars, with no
  additional communication round. Set to `FALSE` to limit server disclosure
  to the three scalars strictly needed for the point estimate (`n_k`,
  `n_k_1`, `tau_k`).

* `heteroscedastic_arms` argument (default `FALSE`) estimates arm-specific
  noise variances \eqn{\sigma_1^2} and \eqn{\sigma_0^2} rather than a pooled
  \eqn{\sigma^2}, adjusting the variance formula in both the delta-method
  and homogeneous closed-form branches.

* `heterogeneous` argument (default `TRUE`) uses the leverage-corrected
  delta-method variance (valid under covariate shift and heterogeneous
  treatment probabilities across sites). Set to `FALSE` for the simpler
  homogeneous closed form.

* `study_effects` argument (default `FALSE`) enables the adjusted
  generalised difference estimator for settings where study membership
  directly shifts the outcome. Appends K-1 integer study-dummy columns
  server-side via `addStudyIdDS()` and adjusts degrees of freedom
  accordingly.

* Input validation helpers `.check_datasources()`, `.check_args()`, and
  `.check_string()` give informative errors on bad arguments.

* `.check_sample_sizes()` warns when either federated arm size falls
  below 5*(d+1), indicating the asymptotic approximation may be unreliable.
