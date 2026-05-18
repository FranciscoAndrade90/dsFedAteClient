library(testthat)

source("setupDSLite.R")
conns <- setupDSLite()

# Happy path ------------------------------------------------------------------

test_that("ds.EstimateAte returns a single finite numeric ATE", {
  result <- ds.EstimateAte(
    datasources     = conns,
    base_name       = "study_data",
    outcome_name    = "Y",
    treatment_name  = "W",
    covariate_names = c("X1", "X2")
  )
  expect_type(result, "double")
  expect_length(result, 1L)
  expect_false(is.na(result))
  expect_false(is.infinite(result))
})

# Sad path: .check_datasources ------------------------------------------------

test_that(".check_datasources rejects a non-list", {
  expect_error(.check_datasources("not_a_list"), "DSConnection")
})

test_that(".check_datasources rejects a list of non-DSConnection objects", {
  expect_error(.check_datasources(list(a = 1, b = 2)), "DSConnection")
})

# Sad path: .check_args -------------------------------------------------------

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

# Sad path: .check_glm_validity -----------------------------------------------

test_that("ds.EstimateAte errors with informative message when disclosure controls block the GLM", {
  # 45 treated / 5 control: control arm passes nfilter.subset (5 >= 3) but
  # fails nfilter.glm (3 params / 0.33 ≈ 10 obs required, 5 < 10).
  set.seed(1)
  n <- 50
  imb <- data.frame(W = sample(c(rep(1, 45), rep(0, 5))), X1 = rnorm(n), X2 = rnorm(n))
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
  DSI::datashield.assign.table(conns = conns_imb, symbol = "study_data", table = "study_data")

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

# Sad path: .check_string -----------------------------------------------------

test_that(".check_string rejects NULL input", {
  expect_error(.check_string(NULL, "my_arg"), "my_arg")
})

test_that(".check_string rejects numeric input", {
  expect_error(.check_string(42, "my_arg"), "my_arg")
})
