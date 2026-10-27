library(DSLite)
# devtools::test() runs from tests/testthat/ so the relative path resolves.
# devtools::check() copies to a temp dir where it does not — fall back to the
# installed package (if available) or skip all tests.
server_src <- "../../../dsFedAte"
if (dir.exists(server_src)) {
  devtools::load_all(server_src, quiet = TRUE)
} else if (requireNamespace("dsFedAte", quietly = TRUE)) {
  library(dsFedAte)
} else {
  setupDSLite <- function() {
    testthat::skip("dsFedAte is not installed and source path is unavailable")
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
  datashield_env <- new.env(parent = emptyenv())
  options(datashield.env = datashield_env)
  assign("dslite.server", dslite.server, envir = datashield_env)

  dslite.server$config(defaultDSConfiguration(include = c("dsBase", "dsFedAte")))
  dslite.server$aggregateMethod("localAteDS",                        "localAteDS")
  dslite.server$aggregateMethod("listDisclosureSettingsDS",           "listDisclosureSettingsDS")
  dslite.server$aggregateMethod("checkPermissivePrivacyControlLevel", "checkPermissivePrivacyControlLevel")
  dslite.server$assignMethod("addStudyIdDS",                         "addStudyIdDS")

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

# Multi-server setup for study-effects tests ----------------------------------
#
# Three servers with different baseline shifts h_k and treatment probabilities
# p_k, chosen so the unadjusted GD estimator is substantially biased:
#
#   h  = c(3, 0, -3)      study effects (h_2 = 0 is baseline but still present)
#   p  = c(0.8, 0.5, 0.2) unequal treatment → H confounds W → Y
#
# Because beta(1) == beta(0), the true ATE is exactly c(1) - c(0) = 0.5
# regardless of X.  The expected bias in the unadjusted estimator is ~2.4,
# so the two estimators are easily distinguished in tests.

setupDSLiteStudyEffects <- function() {
  set.seed(7)
  n  <- 150
  h  <- c(3, 0, -3)
  pk <- c(0.8, 0.5, 0.2)
  c0 <- 1.0;  c1 <- 1.5
  b  <- c(0.3, -0.2)   # same slopes in both arms → ATE = c1 - c0 = 0.5

  make_site <- function(k) {
    X  <- cbind(rnorm(n), rnorm(n))
    W  <- rbinom(n, 1, pk[k])
    Y  <- ifelse(W == 1,
                 c1 + h[k] + X %*% b,
                 c0 + h[k] + X %*% b) + rnorm(n, sd = 0.5)
    data.frame(Y = as.numeric(Y), W = W, X1 = X[, 1], X2 = X[, 2])
  }

  sites <- lapply(seq_along(h), make_site)

  datashield_env <- new.env(parent = emptyenv())
  options(datashield.env = datashield_env)

  servers <- lapply(seq_along(sites), function(k) {
    srv_name <- paste0("dslite.study.server.", k)
    srv      <- DSLite::newDSLiteServer(tables = list(study_data = sites[[k]]))
    srv$config(defaultDSConfiguration(include = c("dsBase", "dsFedAte")))
    srv$aggregateMethod("localAteDS",                        "localAteDS")
    srv$aggregateMethod("listDisclosureSettingsDS",           "listDisclosureSettingsDS")
    srv$aggregateMethod("checkPermissivePrivacyControlLevel", "checkPermissivePrivacyControlLevel")
    srv$assignMethod("addStudyIdDS",                         "addStudyIdDS")
    assign(srv_name, srv, envir = datashield_env)
    srv_name
  })

  builder <- DSI::newDSLoginBuilder()
  for (k in seq_along(servers)) {
    builder$append(
      server = paste0("study", k),
      url    = servers[[k]],
      table  = "study_data",
      driver = "DSLiteDriver"
    )
  }
  conns <- DSI::datashield.login(logins = builder$build(), assign = FALSE)
  DSI::datashield.assign.table(conns = conns, symbol = "study_data",
                               table = "study_data")
  return(conns)
}

