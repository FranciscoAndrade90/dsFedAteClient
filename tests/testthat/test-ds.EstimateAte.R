library(testthat)

# Skip the whole file when dsFedAte is neither loadable from source nor installed.
# This happens in devtools::check() when Suggests are not force-installed.
testthat::skip_if(
  !dir.exists("../../../dsFedAte") && !requireNamespace("dsFedAte", quietly = TRUE),
  "dsFedAte is not installed and source path is unavailable"
)

source("setupDSLite.R")
conns <- setupDSLite()

# Happy path: return structure -------------------------------------------------

test_that("ds.EstimateAte returns a named list with seven elements", {
  result <- ds.EstimateAte(
    datasources     = conns,
    base_name       = "study_data",
    outcome_name    = "Y",
    treatment_name  = "W",
    covariate_names = c("X1", "X2")
  )
  expect_type(result, "list")
  expect_named(result, c("ate", "variance", "ci_lower", "ci_upper",
                         "n_total", "n_k", "ate_k"))
})

test_that("ds.EstimateAte ATE is finite and near the true value of 0.5", {
  result <- ds.EstimateAte(
    datasources     = conns,
    base_name       = "study_data",
    outcome_name    = "Y",
    treatment_name  = "W",
    covariate_names = c("X1", "X2")
  )
  expect_false(is.na(result$ate))
  expect_false(is.infinite(result$ate))
  expect_equal(result$ate, 0.5, tolerance = 0.15)
})

test_that("ds.EstimateAte variance is positive and finite", {
  result <- ds.EstimateAte(
    datasources     = conns,
    base_name       = "study_data",
    outcome_name    = "Y",
    treatment_name  = "W",
    covariate_names = c("X1", "X2")
  )
  expect_gt(result$variance, 0)
  expect_false(is.na(result$variance))
  expect_false(is.infinite(result$variance))
})

test_that("ds.EstimateAte CI contains the true ATE", {
  result <- ds.EstimateAte(
    datasources     = conns,
    base_name       = "study_data",
    outcome_name    = "Y",
    treatment_name  = "W",
    covariate_names = c("X1", "X2")
  )
  expect_lt(result$ci_lower, 0.5)
  expect_gt(result$ci_upper, 0.5)
})

test_that("ds.EstimateAte CI is symmetric around the ATE", {
  result <- ds.EstimateAte(
    datasources     = conns,
    base_name       = "study_data",
    outcome_name    = "Y",
    treatment_name  = "W",
    covariate_names = c("X1", "X2")
  )
  half_width_lower <- result$ate - result$ci_lower
  half_width_upper <- result$ci_upper - result$ate
  expect_equal(half_width_lower, half_width_upper, tolerance = 1e-10)
})

test_that("ds.EstimateAte CI half-width equals 1.96 * sqrt(variance)", {
  result <- ds.EstimateAte(
    datasources     = conns,
    base_name       = "study_data",
    outcome_name    = "Y",
    treatment_name  = "W",
    covariate_names = c("X1", "X2")
  )
  expect_equal(result$ate - result$ci_lower, 1.96 * sqrt(result$variance),
               tolerance = 1e-10)
})

test_that("ds.EstimateAte n_total matches data size", {
  result <- ds.EstimateAte(
    datasources     = conns,
    base_name       = "study_data",
    outcome_name    = "Y",
    treatment_name  = "W",
    covariate_names = c("X1", "X2")
  )
  expect_equal(result$n_total, 200L)
  expect_equal(sum(result$n_k), result$n_total)
})

# compute_variance = FALSE -----------------------------------------------------

test_that("ds.EstimateAte with compute_variance=FALSE returns four elements", {
  result <- ds.EstimateAte(
    datasources      = conns,
    base_name        = "study_data",
    outcome_name     = "Y",
    treatment_name   = "W",
    covariate_names  = c("X1", "X2"),
    compute_variance = FALSE
  )
  expect_named(result, c("ate", "n_total", "n_k", "ate_k"))
})

test_that("ds.EstimateAte compute_variance=FALSE ATE matches full result", {
  full     <- ds.EstimateAte(
    datasources     = conns, base_name = "study_data",
    outcome_name = "Y", treatment_name = "W", covariate_names = c("X1", "X2")
  )
  ate_only <- ds.EstimateAte(
    datasources = conns, base_name = "study_data",
    outcome_name = "Y", treatment_name = "W", covariate_names = c("X1", "X2"),
    compute_variance = FALSE
  )
  expect_equal(ate_only$ate,   full$ate)
  expect_equal(ate_only$n_k,   full$n_k)
  expect_equal(ate_only$ate_k, full$ate_k)
})

# .compute_variance ------------------------------------------------------------

test_that(".compute_variance returns a positive finite scalar", {
  v <- .compute_variance(
    n = 100L, n_k = c(60L, 40L), tau_k = c(0.54, 0.56),
    ate_hat = (60L * 0.54 + 40L * 0.56) / 100L, s_k = c(3.0, 2.5),
    n_1 = 60L, heterogeneous = TRUE, heteroscedastic_arms = FALSE,
    g = c(1, 0.1, -0.05), vcov_treated = (1/60) * diag(3L), vcov_control = (1/40) * diag(3L),
    sigma2_hat = 1.0
  )
  expect_length(v, 1L)
  expect_gt(v, 0)
  expect_false(is.na(v))
  expect_false(is.infinite(v))
})

test_that(".compute_variance parametric term is zero when g = 0", {
  n     <- 100L
  n_k   <- c(60L, 40L)
  tau_k <- c(0.54, 0.46)
  ate_hat <- sum(n_k * tau_k) / n
  s_k   <- c(3.0, 2.5)

  v         <- .compute_variance(
    n = n, n_k = n_k, tau_k = tau_k, ate_hat = ate_hat, s_k = s_k,
    n_1 = 60L, heterogeneous = TRUE, heteroscedastic_arms = FALSE,
    g = rep(0, 3L), vcov_treated = diag(3L), vcov_control = diag(3L),
    sigma2_hat = 1.0
  )
  tss       <- sum(s_k + n_k * (tau_k - ate_hat)^2)
  samp_only <- tss / (n * (n - 1L))
  expect_equal(v, samp_only, tolerance = 1e-10)
})

test_that(".compute_variance sampling term is zero when tau_k = ate_hat and s_k = 0", {
  n   <- 100L
  n_1 <- 60L
  g   <- c(1, 0.1)
  V1  <- (1 / n_1) * diag(2L)
  V0  <- (1 / (n - n_1)) * diag(2L)

  v          <- .compute_variance(
    n = n, n_k = c(n_1, n - n_1), tau_k = c(0.5, 0.5), ate_hat = 0.5, s_k = c(0, 0),
    n_1 = n_1, heterogeneous = TRUE, heteroscedastic_arms = FALSE,
    g = g, vcov_treated = V1, vcov_control = V0,
    sigma2_hat = 1.0
  )
  param_only <- as.numeric(g %*% V1 %*% g) + as.numeric(g %*% V0 %*% g)
  expect_equal(v, param_only, tolerance = 1e-10)
})

test_that(".compute_variance het and hom formulas agree when X_bar = 0 and balanced Gram", {
  # With X_bar = 0, g = c(1, 0), and (G^w)^{-1} = (1/n_w) * I (balanced, centred X),
  # param_var_het = sigma2*(1/n_1 + 1/n_0) = sigma2/(n*p*(1-p)) = param_var_hom.
  # Both sampling terms are zero (s_k = 0, tau_k = ate_hat).
  n          <- 200L
  n_1        <- 100L
  d          <- 1L
  sigma2_hat <- 40.0 / (n - 2L * (d + 1L))
  vcov_1     <- (1 / n_1)       * diag(d + 1L)
  vcov_0     <- (1 / (n - n_1)) * diag(d + 1L)

  v_het <- .compute_variance(
    n = n, n_k = n, tau_k = 0.5, ate_hat = 0.5, s_k = 0.0,
    n_1 = n_1, heterogeneous = TRUE, heteroscedastic_arms = FALSE,
    g = c(1, 0.0), vcov_treated = vcov_1, vcov_control = vcov_0,
    sigma2_hat = sigma2_hat
  )
  v_hom <- .compute_variance(
    n = n, n_k = n, tau_k = 0.5, ate_hat = 0.5, s_k = 0.0,
    n_1 = n_1, heterogeneous = FALSE, heteroscedastic_arms = FALSE,
    sigma2_hat = sigma2_hat
  )
  expect_equal(v_het, v_hom, tolerance = 1e-10)
})

test_that(".compute_variance under !heterogeneous uses arm-specific sigmas when heteroscedastic_arms = TRUE", {
  # Expected formula (!heterogeneous, heteroscedastic_arms):
  #   V_param = sigma2_1 / (n * p) + sigma2_0 / (n * (1 - p)).
  # It must differ from the pooled formula sigma2 / (n * p * (1 - p)) when sigma2_1 != sigma2_0.
  # Both calls share the same (zero) sampling term so the comparison is clean.
  n            <- 200L
  n_1          <- 80L
  p_hat        <- n_1 / n
  sigma2_hat   <- 2.0
  sigma2_hat_1 <- 3.0
  sigma2_hat_0 <- 1.0

  v_het <- .compute_variance(
    n = n, n_k = n, tau_k = 0.5, ate_hat = 0.5, s_k = 0.0,
    n_1 = n_1, heterogeneous = FALSE, heteroscedastic_arms = TRUE,
    sigma2_hat = sigma2_hat,
    sigma2_hat_1 = sigma2_hat_1,
    sigma2_hat_0 = sigma2_hat_0
  )
  v_hom <- .compute_variance(
    n = n, n_k = n, tau_k = 0.5, ate_hat = 0.5, s_k = 0.0,
    n_1 = n_1, heterogeneous = FALSE, heteroscedastic_arms = FALSE,
    sigma2_hat = sigma2_hat
  )

  expected_het <- sigma2_hat_1 / (n * p_hat) + sigma2_hat_0 / (n * (1 - p_hat))
  expected_hom <- sigma2_hat / (n * p_hat * (1 - p_hat))

  expect_equal(v_het, expected_het, tolerance = 1e-10)
  expect_equal(v_hom, expected_hom, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(v_het, v_hom)))
})

# .check_sample_sizes ----------------------------------------------------------

test_that(".check_sample_sizes warns when an arm is below 5*(d+1)", {
  expect_warning(
    .check_sample_sizes(n_k_0 = c(3L, 3L), n_k_1 = c(3L, 3L), d = 2L),
    "Asymptotic"
  )
})

test_that(".check_sample_sizes is silent when both arms are large enough", {
  expect_silent(
    .check_sample_sizes(n_k_0 = c(50L, 50L), n_k_1 = c(50L, 50L), d = 2L)
  )
})

# Sad path: .check_datasources -------------------------------------------------

test_that(".check_datasources rejects a non-list", {
  expect_error(.check_datasources("not_a_list"), "DSConnection")
})

test_that(".check_datasources rejects a list of non-DSConnection objects", {
  expect_error(.check_datasources(list(a = 1, b = 2)), "DSConnection")
})

# Sad path: .check_args --------------------------------------------------------

test_that(".check_args rejects NULL base_name", {
  expect_error(.check_args(NULL, "Y", "W", c("X1")), "base_name")
})

test_that(".check_args rejects non-character outcome_name", {
  expect_error(.check_args("study_data", 123, "W", c("X1")), "outcome_name")
})

test_that(".check_args rejects NULL covariate_names", {
  expect_error(.check_args("study_data", "Y", "W", NULL), "covariate_names")
})

test_that(".check_args rejects an empty covariate_names vector", {
  expect_error(.check_args("study_data", "Y", "W", character(0)), "covariate_names")
})

# Sad path: .check_glm_validity ------------------------------------------------

test_that("ds.EstimateAte errors with informative message when disclosure controls block the GLM", {
  set.seed(1)
  n   <- 50
  imb <- data.frame(W = sample(c(rep(1, 45), rep(0, 5))),
                    X1 = rnorm(n), X2 = rnorm(n))
  imb$Y <- rnorm(n)

  dslite.server.imb <- DSLite::newDSLiteServer(tables = list(study_data = imb))
  imb_env <- new.env(parent = emptyenv())
  options(datashield.env = imb_env)
  assign("dslite.server.imb", dslite.server.imb, envir = imb_env)
  dslite.server.imb$config(defaultDSConfiguration(include = c("dsBase", "dsFedAte")))
  dslite.server.imb$aggregateMethod("localAteDS",                        "localAteDS")
  dslite.server.imb$aggregateMethod("listDisclosureSettingsDS",           "listDisclosureSettingsDS")
  dslite.server.imb$aggregateMethod("checkPermissivePrivacyControlLevel", "checkPermissivePrivacyControlLevel")

  builder <- DSI::newDSLoginBuilder()
  builder$append(server = "server_imb", url = "dslite.server.imb",
                 table = "study_data", driver = "DSLiteDriver")
  conns_imb <- DSI::datashield.login(logins = builder$build(), assign = FALSE)
  DSI::datashield.assign.table(conns = conns_imb, symbol = "study_data",
                               table = "study_data")

  expect_error(
    ds.EstimateAte(
      datasources     = conns_imb,
      base_name       = "study_data",
      outcome_name    = "Y",
      treatment_name  = "W",
      covariate_names = c("X1", "X2")
    ),
    "DataSHIELD blocked the GLM"
  )

  DSI::datashield.logout(conns_imb)
})

# Sad path: .check_string ------------------------------------------------------

test_that(".check_string rejects NULL input", {
  expect_error(.check_string(NULL, "my_arg"), "my_arg")
})

test_that(".check_string rejects numeric input", {
  expect_error(.check_string(42, "my_arg"), "my_arg")
})

# study_effects ----------------------------------------------------------------

conns_se <- setupDSLiteStudyEffects()
TRUE_ATE  <- 0.5   # c(1) - c(0) with beta(1) == beta(0)

test_that("ds.EstimateAte study_effects=TRUE returns seven-element list", {
  result <- ds.EstimateAte(
    datasources     = conns_se,
    base_name       = "study_data",
    outcome_name    = "Y",
    treatment_name  = "W",
    covariate_names = c("X1", "X2"),
    study_effects   = TRUE
  )
  expect_named(result, c("ate", "variance", "ci_lower", "ci_upper",
                         "n_total", "n_k", "ate_k"))
})

test_that("ds.EstimateAte study_effects=TRUE ATE is close to the true value", {
  result <- ds.EstimateAte(
    datasources     = conns_se,
    base_name       = "study_data",
    outcome_name    = "Y",
    treatment_name  = "W",
    covariate_names = c("X1", "X2"),
    study_effects   = TRUE
  )
  expect_equal(result$ate, TRUE_ATE, tolerance = 0.25)
})

test_that("ds.EstimateAte study_effects=TRUE CI covers the true ATE", {
  result <- ds.EstimateAte(
    datasources     = conns_se,
    base_name       = "study_data",
    outcome_name    = "Y",
    treatment_name  = "W",
    covariate_names = c("X1", "X2"),
    study_effects   = TRUE
  )
  expect_lt(result$ci_lower, TRUE_ATE)
  expect_gt(result$ci_upper, TRUE_ATE)
})

test_that("ds.EstimateAte study_effects=FALSE is substantially biased", {
  # With h = c(3, 0, -3) and p = c(0.8, 0.5, 0.2) the expected bias is ~2.4
  result <- ds.EstimateAte(
    datasources     = conns_se,
    base_name       = "study_data",
    outcome_name    = "Y",
    treatment_name  = "W",
    covariate_names = c("X1", "X2"),
    study_effects   = FALSE
  )
  expect_gt(abs(result$ate - TRUE_ATE), 1.0)
})

test_that("ds.EstimateAte study_effects=TRUE variance is positive and finite", {
  result <- ds.EstimateAte(
    datasources     = conns_se,
    base_name       = "study_data",
    outcome_name    = "Y",
    treatment_name  = "W",
    covariate_names = c("X1", "X2"),
    study_effects   = TRUE
  )
  expect_gt(result$variance, 0)
  expect_false(is.na(result$variance))
  expect_false(is.infinite(result$variance))
})

test_that("ds.EstimateAte study_effects=TRUE with compute_variance=FALSE returns four elements", {
  result <- ds.EstimateAte(
    datasources      = conns_se,
    base_name        = "study_data",
    outcome_name     = "Y",
    treatment_name   = "W",
    covariate_names  = c("X1", "X2"),
    study_effects    = TRUE,
    compute_variance = FALSE
  )
  expect_named(result, c("ate", "n_total", "n_k", "ate_k"))
})

test_that("ds.EstimateAte study_effects=TRUE with K=1 warns and falls back", {
  expect_warning(
    result <- ds.EstimateAte(
      datasources     = conns,      # single-server fixture from setupDSLite()
      base_name       = "study_data",
      outcome_name    = "Y",
      treatment_name  = "W",
      covariate_names = c("X1", "X2"),
      study_effects   = TRUE
    ),
    "no effect"
  )
  # Falls back to unadjusted — should still return a valid result
  expect_named(result, c("ate", "variance", "ci_lower", "ci_upper",
                         "n_total", "n_k", "ate_k"))
})

# .add_study_id ----------------------------------------------------------------

# heterogeneous = TRUE --------------------------------------------------------

test_that("ds.EstimateAte heterogeneous=TRUE returns a seven-element list", {
  result <- ds.EstimateAte(
    datasources     = conns,
    base_name       = "study_data",
    outcome_name    = "Y",
    treatment_name  = "W",
    covariate_names = c("X1", "X2"),
    heterogeneous   = TRUE
  )
  expect_named(result, c("ate", "variance", "ci_lower", "ci_upper",
                         "n_total", "n_k", "ate_k"))
})

test_that("ds.EstimateAte heterogeneous=TRUE ATE is close to the true value of 0.5", {
  result <- ds.EstimateAte(
    datasources     = conns,
    base_name       = "study_data",
    outcome_name    = "Y",
    treatment_name  = "W",
    covariate_names = c("X1", "X2"),
    heterogeneous   = TRUE
  )
  expect_equal(result$ate, 0.5, tolerance = 0.15)
})

test_that("ds.EstimateAte heterogeneous=TRUE variance is positive and finite", {
  result <- ds.EstimateAte(
    datasources     = conns,
    base_name       = "study_data",
    outcome_name    = "Y",
    treatment_name  = "W",
    covariate_names = c("X1", "X2"),
    heterogeneous   = TRUE
  )
  expect_gt(result$variance, 0)
  expect_false(is.na(result$variance))
  expect_false(is.infinite(result$variance))
})

test_that("ds.EstimateAte heterogeneous=TRUE and FALSE give the same ATE", {
  r_het <- ds.EstimateAte(
    datasources     = conns,
    base_name       = "study_data",
    outcome_name    = "Y",
    treatment_name  = "W",
    covariate_names = c("X1", "X2"),
    heterogeneous   = TRUE
  )
  r_hom <- ds.EstimateAte(
    datasources     = conns,
    base_name       = "study_data",
    outcome_name    = "Y",
    treatment_name  = "W",
    covariate_names = c("X1", "X2"),
    heterogeneous   = FALSE
  )
  expect_equal(r_het$ate, r_hom$ate)
})

test_that("ds.EstimateAte heterogeneous=TRUE with compute_variance=FALSE returns four elements", {
  result <- ds.EstimateAte(
    datasources      = conns,
    base_name        = "study_data",
    outcome_name     = "Y",
    treatment_name   = "W",
    covariate_names  = c("X1", "X2"),
    heterogeneous    = TRUE,
    compute_variance = FALSE
  )
  expect_named(result, c("ate", "n_total", "n_k", "ate_k"))
})

test_that("ds.EstimateAte heterogeneous=TRUE CI contains the true ATE", {
  result <- ds.EstimateAte(
    datasources     = conns,
    base_name       = "study_data",
    outcome_name    = "Y",
    treatment_name  = "W",
    covariate_names = c("X1", "X2"),
    heterogeneous   = TRUE
  )
  expect_lt(result$ci_lower, 0.5)
  expect_gt(result$ci_upper, 0.5)
})

# heteroscedastic_arms = TRUE ---------------------------------------------------

test_that("ds.EstimateAte heteroscedastic_arms=TRUE returns a seven-element list", {
  result <- ds.EstimateAte(
    datasources        = conns,
    base_name          = "study_data",
    outcome_name       = "Y",
    treatment_name     = "W",
    covariate_names    = c("X1", "X2"),
    heteroscedastic_arms = TRUE
  )
  expect_named(result, c("ate", "variance", "ci_lower", "ci_upper",
                         "n_total", "n_k", "ate_k"))
})

test_that("ds.EstimateAte heteroscedastic_arms=TRUE ATE equals heteroscedastic_arms=FALSE ATE", {
  r_pooled <- ds.EstimateAte(
    datasources        = conns,
    base_name          = "study_data",
    outcome_name       = "Y",
    treatment_name     = "W",
    covariate_names    = c("X1", "X2"),
    heteroscedastic_arms = FALSE
  )
  r_split <- ds.EstimateAte(
    datasources        = conns,
    base_name          = "study_data",
    outcome_name       = "Y",
    treatment_name     = "W",
    covariate_names    = c("X1", "X2"),
    heteroscedastic_arms = TRUE
  )
  expect_equal(r_split$ate, r_pooled$ate)
})

test_that("ds.EstimateAte heteroscedastic_arms=TRUE variance is positive and finite", {
  result <- ds.EstimateAte(
    datasources        = conns,
    base_name          = "study_data",
    outcome_name       = "Y",
    treatment_name     = "W",
    covariate_names    = c("X1", "X2"),
    heteroscedastic_arms = TRUE
  )
  expect_gt(result$variance, 0)
  expect_false(is.na(result$variance))
  expect_false(is.infinite(result$variance))
})

test_that("ds.EstimateAte heterogeneous=FALSE heteroscedastic_arms=TRUE variance is positive and finite", {
  result <- ds.EstimateAte(
    datasources        = conns,
    base_name          = "study_data",
    outcome_name       = "Y",
    treatment_name     = "W",
    covariate_names    = c("X1", "X2"),
    heterogeneous      = FALSE,
    heteroscedastic_arms = TRUE
  )
  expect_gt(result$variance, 0)
  expect_false(is.na(result$variance))
  expect_false(is.infinite(result$variance))
})

# quiet = TRUE -----------------------------------------------------------------

test_that("ds.EstimateAte quiet=TRUE returns the same result as quiet=FALSE", {
  r_noisy <- ds.EstimateAte(
    datasources     = conns,
    base_name       = "study_data",
    outcome_name    = "Y",
    treatment_name  = "W",
    covariate_names = c("X1", "X2"),
    quiet           = FALSE
  )
  r_quiet <- ds.EstimateAte(
    datasources     = conns,
    base_name       = "study_data",
    outcome_name    = "Y",
    treatment_name  = "W",
    covariate_names = c("X1", "X2"),
    quiet           = TRUE
  )
  expect_equal(r_quiet$ate,      r_noisy$ate)
  expect_equal(r_quiet$variance, r_noisy$variance)
  expect_equal(r_quiet$n_total,  r_noisy$n_total)
})

test_that("ds.EstimateAte quiet=TRUE suppresses messages", {
  expect_no_message(
    ds.EstimateAte(
      datasources     = conns,
      base_name       = "study_data",
      outcome_name    = "Y",
      treatment_name  = "W",
      covariate_names = c("X1", "X2"),
      quiet           = TRUE
    )
  )
})

# .add_study_id ----------------------------------------------------------------

test_that(".add_study_id creates a new object with study_id2 and study_id3, leaving original unchanged", {
  .add_study_id(conns_se, "study_data", length(conns_se), "study_data_adj_test")
  adj_cols  <- dsBaseClient::ds.colnames("study_data_adj_test", datasources = conns_se)
  orig_cols <- dsBaseClient::ds.colnames("study_data",          datasources = conns_se)
  for (srv in names(adj_cols)) {
    expect_true("study_id2"  %in% adj_cols[[srv]])
    expect_true("study_id3"  %in% adj_cols[[srv]])
    expect_false("study_id"  %in% adj_cols[[srv]])    # no factor column
    expect_false("study_id2" %in% orig_cols[[srv]])   # original untouched
  }
})
