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

# Sad path: .check_string -----------------------------------------------------

test_that(".check_string rejects NULL input", {
  expect_error(.check_string(NULL, "my_arg"), "my_arg")
})

test_that(".check_string rejects numeric input", {
  expect_error(.check_string(42, "my_arg"), "my_arg")
})
