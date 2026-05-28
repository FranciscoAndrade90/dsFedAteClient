library(testthat)

# Skip the whole file when dsRctFed is neither loadable from source nor installed.
# This happens in devtools::check() when Suggests are not force-installed.
testthat::skip_if(
  !dir.exists("../../../dsRctFed") && !requireNamespace("dsRctFed", quietly = TRUE),
  "dsRctFed is not installed and source path is unavailable"
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

test_that(".compute_variance returns a positive scalar", {
  n       <- 100L
  n_1     <- 50L
  d       <- 2L
  rss_k   <- c(20.0, 15.0)
  s_k     <- c(3.0, 2.5)
  n_k     <- c(60L, 40L)
  tau_k   <- c(0.54, 0.56)
  ate_hat <- sum(n_k * tau_k) / n

  v <- .compute_variance(n, n_1, rss_k, s_k, n_k, tau_k, ate_hat, d)
  expect_length(v, 1L)
  expect_gt(v, 0)
  expect_false(is.na(v))
})

test_that(".compute_variance first term equals sigma2 / (n * p * (1-p))", {
  n       <- 100L
  n_1     <- 50L
  d       <- 1L
  rss_k   <- c(19.6)   # sigma2_hat = 19.6 / (100 - 4)
  s_k     <- 0.0       # zero within-site scatter
  n_k     <- 100L
  tau_k   <- 0.5       # τ̂_k = τ̂ → between-site term = 0
  ate_hat <- 0.5

  v          <- .compute_variance(n, n_1, rss_k, s_k, n_k, tau_k, ate_hat, d)
  sigma2_hat <- rss_k / (n - 2L * (d + 1L))
  p_hat      <- n_1 / n
  expected   <- sigma2_hat / (n * p_hat * (1 - p_hat))
  expect_equal(v, expected, tolerance = 1e-10)
})

test_that(".compute_variance second term is zero when all tau_k equal ate_hat", {
  n       <- 100L
  n_1     <- 50L
  d       <- 2L
  rss_k   <- c(10.0, 8.0)
  n_k     <- c(60L, 40L)
  tau_k   <- c(0.5, 0.5)   # τ̂_k − τ̂ = 0; s_k_scalar also zero
  ate_hat <- 0.5

  v_second_zero <- .compute_variance(n, n_1, rss_k, c(0, 0), n_k, tau_k, ate_hat, d)
  sigma2_hat    <- sum(rss_k) / (n - 2L * (d + 1L))
  p_hat         <- n_1 / n
  first_only    <- sigma2_hat / (n * p_hat * (1 - p_hat))
  expect_equal(v_second_zero, first_only, tolerance = 1e-10)
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
  assign("dslite.server.imb", dslite.server.imb, envir = .GlobalEnv)
  dslite.server.imb$config(defaultDSConfiguration(include = c("dsBase", "dsRctFed")))
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
