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

test_that(".compute_variance heterogeneous=TRUE returns a positive finite scalar", {
  n       <- 100L
  n_1     <- 50L
  d       <- 2L
  rss_k   <- c(10.0, 8.0)
  s_k     <- c(3.0, 2.5)
  n_k     <- c(60L, 40L)
  tau_k   <- c(0.54, 0.46)
  ate_hat <- sum(n_k * tau_k) / n
  # vcov_w = (G^w)^{-1} (unscaled), as returned by ds.glm's VarCovMatrix
  vcov_1  <- (1 / 60) * diag(d + 1L)
  vcov_0  <- (1 / 40) * diag(d + 1L)
  X_bar   <- c(0.1, -0.05)

  v <- .compute_variance(n, n_1, rss_k, s_k, n_k, tau_k, ate_hat, d,
         heterogeneous = TRUE,
         vcov_treated  = vcov_1,
         vcov_control  = vcov_0,
         X_bar         = X_bar)
  expect_length(v, 1L)
  expect_gt(v, 0)
  expect_false(is.na(v))
  expect_false(is.infinite(v))
})

test_that(".compute_variance heterogeneous=TRUE Term(I) reduces to homogeneous when delta=0", {
  # When vcov_w = sigma2_w * (G^w)^{-1} and X_bar = arm means = global mean,
  # Term(I) = sigma2_hat * X_bar' (G1_inv + G0_inv) X_bar'.
  # For a simple balanced case with identity Gram (so G^{-1} = I/n^w for a
  # standardised design), the quadratic form gives sigma2_hat/(n*p*(1-p)).
  # Here we just verify numerical agreement against the closed-form path.
  n       <- 200L
  n_1     <- 100L
  d       <- 1L
  rss_k   <- 40.0
  s_k     <- 0.0
  n_k     <- 200L
  tau_k   <- 0.5
  ate_hat <- 0.5

  sigma2_hat <- rss_k / (n - 2L * (d + 1L))
  p_hat      <- n_1 / n

  # When G^(w) = n^(w) * I_{d+1} (identity Gram, as for centred orthogonal X),
  # (G^(w))^{-1} = (1/n^(w)) * I_{d+1} — this is what ds.glm returns as VarCovMatrix.
  # With X_bar = 0, Lemma 2 gives X_bar' (G1_inv + G0_inv) X_bar' = 1/n_1 + 1/n_0 = 1/(n*p*(1-p)).
  # So Term I = sigma2_hat / (n * p * (1-p)), recovering the homogeneous closed form.
  vcov_1  <- (1 / n_1)      * diag(d + 1L)
  vcov_0  <- (1 / (n - n_1)) * diag(d + 1L)
  X_bar   <- 0.0

  v_het  <- .compute_variance(n, n_1, rss_k, s_k, n_k, tau_k, ate_hat, d,
               heterogeneous = TRUE,
               vcov_treated  = vcov_1,
               vcov_control  = vcov_0,
               X_bar         = X_bar)
  v_hom  <- .compute_variance(n, n_1, rss_k, s_k, n_k, tau_k, ate_hat, d)

  expect_equal(v_het, v_hom, tolerance = 1e-10)
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

# .add_study_id ----------------------------------------------------------------

# outcome_type = "binary" ------------------------------------------------------

conns_bin <- setupDSLiteBinary()

test_that("ds.EstimateAte binary returns four-element list", {
  result <- ds.EstimateAte(
    datasources      = conns_bin,
    base_name        = "study_data",
    outcome_name     = "Y",
    treatment_name   = "W",
    covariate_names  = c("X1", "X2"),
    outcome_type     = "binary",
    compute_variance = FALSE
  )
  expect_named(result, c("ate", "n_total", "n_k", "ate_k"))
})

test_that("ds.EstimateAte binary ATE is a finite value in (-1, 1)", {
  result <- ds.EstimateAte(
    datasources      = conns_bin,
    base_name        = "study_data",
    outcome_name     = "Y",
    treatment_name   = "W",
    covariate_names  = c("X1", "X2"),
    outcome_type     = "binary",
    compute_variance = FALSE
  )
  expect_false(is.na(result$ate))
  expect_gt(result$ate, -1)
  expect_lt(result$ate,  1)
})

test_that("ds.EstimateAte binary with compute_variance=TRUE returns seven elements", {
  result <- ds.EstimateAte(
    datasources      = conns_bin,
    base_name        = "study_data",
    outcome_name     = "Y",
    treatment_name   = "W",
    covariate_names  = c("X1", "X2"),
    outcome_type     = "binary",
    compute_variance = TRUE
  )
  expect_named(result, c("ate", "variance", "ci_lower", "ci_upper",
                         "n_total", "n_k", "ate_k"))
})

test_that("ds.EstimateAte binary variance is positive and finite", {
  result <- ds.EstimateAte(
    datasources      = conns_bin,
    base_name        = "study_data",
    outcome_name     = "Y",
    treatment_name   = "W",
    covariate_names  = c("X1", "X2"),
    outcome_type     = "binary",
    compute_variance = TRUE
  )
  expect_gt(result$variance, 0)
  expect_false(is.na(result$variance))
  expect_false(is.infinite(result$variance))
})

test_that("ds.EstimateAte binary CI half-width equals 1.96 * sqrt(variance)", {
  result <- ds.EstimateAte(
    datasources      = conns_bin,
    base_name        = "study_data",
    outcome_name     = "Y",
    treatment_name   = "W",
    covariate_names  = c("X1", "X2"),
    outcome_type     = "binary",
    compute_variance = TRUE
  )
  expect_equal(result$ate - result$ci_lower, 1.96 * sqrt(result$variance),
               tolerance = 1e-10)
})

test_that("ds.EstimateAte binary n_total matches data size", {
  result <- ds.EstimateAte(
    datasources      = conns_bin,
    base_name        = "study_data",
    outcome_name     = "Y",
    treatment_name   = "W",
    covariate_names  = c("X1", "X2"),
    outcome_type     = "binary",
    compute_variance = FALSE
  )
  expect_equal(result$n_total, 200L)
  expect_equal(sum(result$n_k), result$n_total)
})

# .compute_variance_binary -----------------------------------------------------

make_binary_inputs <- function() {
  set.seed(99)
  n      <- 100L
  n_k    <- c(60L, 40L)
  alpha1 <- 0.5;  alpha0 <- 0.0
  beta1  <- c(0.4, -0.2);  beta0 <- c(0.3, -0.1)

  make_site_data <- function(nk) {
    X   <- cbind(rnorm(nk), rnorm(nk))
    W   <- rbinom(nk, 1, 0.5)
    mu1 <- as.numeric(plogis(alpha1 + X %*% beta1))
    mu0 <- as.numeric(plogis(alpha0 + X %*% beta0))
    list(X = X, W = W, mu1 = mu1, mu0 = mu0)
  }
  sites <- lapply(n_k, make_site_data)

  # Aggregate Hessians for vcov (sum over sites, then invert)
  H1_total <- H0_total <- matrix(0, 3L, 3L)
  g1_raw_k <- g0_raw_k <- vector("list", 2L)
  s_k_delta <- numeric(2L)
  tau_k     <- numeric(2L)

  for (k in seq_along(sites)) {
    s    <- sites[[k]]
    Xt   <- cbind(1, s$X)
    v1   <- s$mu1 * (1 - s$mu1)
    v0   <- s$mu0 * (1 - s$mu0)
    idx1 <- s$W == 1L
    H1_total      <- H1_total + crossprod(Xt[ idx1, ], v1[ idx1] * Xt[ idx1, ])
    H0_total      <- H0_total + crossprod(Xt[!idx1, ], v0[!idx1] * Xt[!idx1, ])
    g1_raw_k[[k]] <- colSums(v1 * Xt)
    g0_raw_k[[k]] <- colSums(v0 * Xt)
    delta_i       <- s$mu1 - s$mu0
    tau_k[k]      <- mean(delta_i)
    s_k_delta[k]  <- sum((delta_i - tau_k[k])^2)
  }
  ate_hat      <- sum(n_k * tau_k) / n
  vcov_treated <- solve(H1_total)
  vcov_control <- solve(H0_total)

  list(n = n, n_k = n_k, tau_k = tau_k, ate_hat = ate_hat,
       vcov_treated = vcov_treated, vcov_control = vcov_control,
       g1_raw_k = g1_raw_k, g0_raw_k = g0_raw_k,
       s_k_delta = s_k_delta)
}

test_that(".compute_variance_binary returns a positive finite scalar", {
  inp <- make_binary_inputs()
  v   <- .compute_variance_binary(
    n            = inp$n,
    n_k          = inp$n_k,
    tau_k        = inp$tau_k,
    ate_hat      = inp$ate_hat,
    vcov_treated = inp$vcov_treated,
    vcov_control = inp$vcov_control,
    g1_raw_k     = inp$g1_raw_k,
    g0_raw_k     = inp$g0_raw_k,
    s_k_delta    = inp$s_k_delta
  )
  expect_length(v, 1L)
  expect_gt(v, 0)
  expect_false(is.na(v))
  expect_false(is.infinite(v))
})

test_that(".compute_variance_binary matches manual pooled formula (single site)", {
  # Single site: vcov comes from glm(), which equals H^{-1} for binomial.
  # This mirrors the ds.glm(viewVarCov=TRUE) path used in production.
  set.seed(7)
  n  <- 100L
  X  <- cbind(rnorm(n), rnorm(n))
  W  <- rbinom(n, 1, 0.5)
  df <- data.frame(
    Y  = rbinom(n, 1, plogis(-0.2 + X %*% c(0.4, -0.2))),
    W  = W, X1 = X[, 1], X2 = X[, 2]
  )
  fit1 <- glm(Y ~ X1 + X2, data = df[df$W == 1, ], family = binomial)
  fit0 <- glm(Y ~ X1 + X2, data = df[df$W == 0, ], family = binomial)

  Xt  <- cbind(1, X)
  mu1 <- as.numeric(plogis(Xt %*% coef(fit1)))
  mu0 <- as.numeric(plogis(Xt %*% coef(fit0)))
  v1  <- mu1 * (1 - mu1)
  v0  <- mu0 * (1 - mu0)
  delta <- mu1 - mu0
  ate   <- mean(delta)

  g1 <-  colMeans(v1 * Xt)
  g0 <- -colMeans(v0 * Xt)

  v_manual <- as.numeric(
    crossprod(g1, vcov(fit1) %*% g1) + crossprod(g0, vcov(fit0) %*% g0)
  ) + var(delta) / n

  v_helper <- .compute_variance_binary(
    n            = n,
    n_k          = n,
    tau_k        = ate,
    ate_hat      = ate,
    vcov_treated = vcov(fit1),
    vcov_control = vcov(fit0),
    g1_raw_k     = list(colSums(v1 * Xt)),
    g0_raw_k     = list(colSums(v0 * Xt)),
    s_k_delta    = sum((delta - ate)^2)
  )
  expect_equal(v_helper, v_manual, tolerance = 1e-10)
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
