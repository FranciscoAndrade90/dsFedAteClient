library(DSLite)
# devtools::test() runs from tests/testthat/ so the relative path resolves.
# devtools::check() copies to a temp dir where it does not — fall back to the
# installed package (if available) or skip all tests.
server_src <- "../../../dsRctFed"
if (dir.exists(server_src)) {
  devtools::load_all(server_src, quiet = TRUE)
} else if (requireNamespace("dsRctFed", quietly = TRUE)) {
  library(dsRctFed)
} else {
  setupDSLite <- function() {
    testthat::skip("dsRctFed is not installed and source path is unavailable")
  }
}

setupDSLite <- function() {
  set.seed(42)
  n          <- 200
  study_data <- data.frame(
    W  = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )
  study_data$Y <- 0.5 * study_data$W +
    0.3 * study_data$X1 -
    0.2 * study_data$X2 +
    rnorm(n, sd = 0.5)

  dslite.server <- DSLite::newDSLiteServer(
    tables = list(study_data = study_data)
  )
  # Assign to the global environment so DSLite can resolve the URL
  # and sessions remain accessible after this function returns.
  assign("dslite.server", dslite.server, envir = .GlobalEnv)
  options(datashield.env = .GlobalEnv)

  dslite.server$config(defaultDSConfiguration(include = c("dsBase", "dsRctFed")))
  dslite.server$aggregateMethod("localAteDS",                      "localAteDS")
  dslite.server$aggregateMethod("listDisclosureSettingsDS",         "listDisclosureSettingsDS")
  dslite.server$aggregateMethod("checkPermissivePrivacyControlLevel", "checkPermissivePrivacyControlLevel")

  options(
    datashield.privacyControlLevel = "banana",
    nfilter.tab                    = 3,
    nfilter.subset                 = 3,
    nfilter.glm                    = 0.33,
    nfilter.string                 = 80,
    nfilter.stringShort            = 20,
    nfilter.kNN                    = 3,
    nfilter.levels.density         = 0.33,
    nfilter.levels.max             = 40,
    nfilter.noise                  = 0.25,
    nfilter.privacy.old            = 5
  )

  builder <- DSI::newDSLoginBuilder()
  builder$append(
    server = "server_1",
    url    = "dslite.server",
    table  = "study_data",
    driver = "DSLiteDriver"
  )
  conns <- DSI::datashield.login(logins = builder$build(), assign = FALSE)
  DSI::datashield.assign.table(conns = conns, symbol = "study_data", table = "study_data")
  return(conns)
}
