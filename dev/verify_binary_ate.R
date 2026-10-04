library(DSLite)
library(DSI)
library(dsBase)
library(dsBaseClient)
library(devtools)

devtools::load_all(".")
devtools::load_all("../dsFedAte")

# ============================================================================
# Data-generating process
# True model: logit P(Y=1 | X, W=w) = alpha_w + beta' X
# with alpha0=1.0, alpha1=1.5, beta=(0.3, -0.2) identical in both arms.
# True ATE = E[plogis(1.5 + X'beta) - plogis(1.0 + X'beta)] > 0.
# ============================================================================

set.seed(42)

K           <- 3
ALPHA0      <- 1.0;  ALPHA1 <- 1.5
BETA        <- c(0.3, -0.2)
n_per_study <- c(200, 150, 250)
p_treat     <- c(0.5, 0.6, 0.4)

make_site <- function(k) {
  n  <- n_per_study[k]
  X  <- cbind(X1 = rnorm(n), X2 = rnorm(n))
  W  <- rbinom(n, 1, p_treat[k])
  lp <- ifelse(W == 1, ALPHA1, ALPHA0) + X %*% BETA
  Y  <- rbinom(n, 1, plogis(lp))
  data.frame(Y = Y, W = W, X1 = X[, 1], X2 = X[, 2])
}

data_list   <- lapply(seq_len(K), make_site)
pooled_data <- do.call(rbind, data_list)

cat("=== Data summary ===\n")
for (k in seq_len(K)) {
  cat(sprintf("  Site %d: n=%d, n_treated=%d, prevalence=%.3f\n",
              k, nrow(data_list[[k]]),
              sum(data_list[[k]]$W),
              mean(data_list[[k]]$Y)))
}
cat("\n")

# ============================================================================
# Pooled (centralised) G-computation
# ============================================================================

covariate_names <- c("X1", "X2")

pooled_treated <- pooled_data[pooled_data$W == 1, ]
pooled_control <- pooled_data[pooled_data$W == 0, ]

model1_pooled <- glm(Y ~ X1 + X2, family = binomial, data = pooled_treated)
model0_pooled <- glm(Y ~ X1 + X2, family = binomial, data = pooled_control)

alpha1_hat <- coef(model1_pooled)["(Intercept)"]
alpha0_hat <- coef(model0_pooled)["(Intercept)"]
beta1_hat  <- coef(model1_pooled)[covariate_names]
beta0_hat  <- coef(model0_pooled)[covariate_names]

X_pooled     <- as.matrix(pooled_data[, covariate_names])
mu1          <- plogis(alpha1_hat + X_pooled %*% beta1_hat)
mu0          <- plogis(alpha0_hat + X_pooled %*% beta0_hat)
ate_pooled   <- mean(mu1 - mu0)

cat(sprintf("Pooled G-computation ATE: %.6f\n\n", ate_pooled))

# ============================================================================
# Federated estimate via DSLite
# ============================================================================

options(
  datashield.env             = .GlobalEnv,
  datashield.privacyControlLevel = "banana",
  nfilter.tab                = 3,
  nfilter.subset             = 3,
  nfilter.glm                = 0.33,
  nfilter.string             = 80,
  nfilter.stringShort        = 20,
  nfilter.kNN                = 3,
  nfilter.levels.density     = 0.33,
  nfilter.levels.max         = 40,
  nfilter.noise              = 0.25,
  nfilter.privacy.old        = 5
)

builder <- DSI::newDSLoginBuilder()
for (k in seq_len(K)) {
  srv_name <- paste0("dslite.server.bin.", k)
  srv      <- DSLite::newDSLiteServer(tables = list(study_data = data_list[[k]]))
  srv$config(defaultDSConfiguration(include = c("dsBase", "dsFedAte")))
  srv$aggregateMethod("localAteDS",                        "localAteDS")
  srv$aggregateMethod("listDisclosureSettingsDS",           "listDisclosureSettingsDS")
  srv$aggregateMethod("checkPermissivePrivacyControlLevel", "checkPermissivePrivacyControlLevel")
  srv$assignMethod("addStudyIdDS",                         "addStudyIdDS")
  assign(srv_name, srv, envir = .GlobalEnv)
  builder$append(server = paste0("site", k), url = srv_name,
                 table = "study_data", driver = "DSLiteDriver")
}

conns <- DSI::datashield.login(logins = builder$build(), assign = FALSE)
DSI::datashield.assign.table(conns = conns, symbol = "study_data",
                              table = "study_data")

result_fed <- ds.EstimateAte(
  datasources      = conns,
  base_name        = "study_data",
  outcome_name     = "Y",
  treatment_name   = "W",
  covariate_names  = covariate_names,
  outcome_type     = "binary",
  compute_variance = FALSE
)

cat(sprintf("Federated ATE:            %.6f\n", result_fed$ate))
cat(sprintf("  per-site: %s\n\n",
            paste(sprintf("site%d=%.6f", seq_len(K), result_fed$ate_k),
                  collapse=", ")))

# ============================================================================
# Distance
# ============================================================================

distance <- abs(result_fed$ate - ate_pooled)
cat(sprintf("Absolute distance (federated vs pooled): %.2e\n", distance))
cat(sprintf("(Should be <= floating-point precision, ~1e-14)\n\n"))

DSI::datashield.logout(conns)
