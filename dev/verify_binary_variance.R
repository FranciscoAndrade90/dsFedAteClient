library(DSLite)
library(DSI)
library(dsBase)
library(dsBaseClient)
library(devtools)

devtools::load_all(".")
devtools::load_all("../dsRctFed")

# ============================================================================
# DGP: K=3 sites, heterogeneous n and treatment probability
# logit P(Y=1|X,W=w) = alpha_w + beta'X  (same beta in both arms)
# ============================================================================

set.seed(42)

K           <- 3
ALPHA0      <- 1.0;  ALPHA1 <- 1.5
BETA        <- c(0.3, -0.2)
n_per_study <- c(200, 150, 250)
p_treat     <- c(0.5, 0.6, 0.4)
covariate_names <- c("X1", "X2")

make_site <- function(k) {
  n  <- n_per_study[k]
  X  <- cbind(X1 = rnorm(n), X2 = rnorm(n))
  W  <- rbinom(n, 1, p_treat[k])
  lp <- ifelse(W == 1, ALPHA1, ALPHA0) + X %*% BETA
  Y  <- rbinom(n, 1, plogis(lp))
  data.frame(Y = Y, W = W, X1 = X[, 1], X2 = X[, 2])
}

data_list <- lapply(seq_len(K), make_site)
dat       <- do.call(rbind, data_list)
n         <- nrow(dat)

cat(sprintf("Total n = %d across %d sites\n\n", n, K))

# ============================================================================
# Set up DSLite servers
# ============================================================================

options(
  datashield.env                 = .GlobalEnv,
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
for (k in seq_len(K)) {
  srv_name <- paste0("dslite.server.bin.var.", k)
  srv      <- DSLite::newDSLiteServer(tables = list(study_data = data_list[[k]]))
  srv$config(defaultDSConfiguration(include = c("dsBase", "dsRctFed")))
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

# ============================================================================
# Step 1: run ds.EstimateAte() — this also leaves study_data_treated and
# study_data_control on every server, which we reuse below.
# ============================================================================

result_fed <- ds.EstimateAte(
  datasources      = conns,
  base_name        = "study_data",
  outcome_name     = "Y",
  treatment_name   = "W",
  covariate_names  = covariate_names,
  outcome_type     = "binary",
  compute_variance = TRUE
)

cat("=== Federated ds.EstimateAte() ===\n")
cat(sprintf("  ATE:       %.10f\n",   result_fed$ate))
cat(sprintf("  Variance:  %.10f\n\n", result_fed$variance))

# ============================================================================
# Step 2: extract the exact same federated coefficients by re-running
# ds.glm() on the arm subsets that ds.EstimateAte() already created.
# Both computations now use the same θ̂, so any remaining discrepancy
# is purely due to the federated decomposition — not IRLS convergence.
# ============================================================================

formula_str <- paste("Y ~", paste(covariate_names, collapse = " + "))

model1_fed <- suppressWarnings(dsBaseClient::ds.glm(
  formula = formula_str, data = "study_data_treated",
  family = "binomial", datasources = conns
))
model0_fed <- suppressWarnings(dsBaseClient::ds.glm(
  formula = formula_str, data = "study_data_control",
  family = "binomial", datasources = conns
))

alpha1_fed <- model1_fed$coefficients["(Intercept)", "Estimate"]
alpha0_fed <- model0_fed$coefficients["(Intercept)", "Estimate"]
beta1_fed  <- model1_fed$coefficients[covariate_names, "Estimate"]
beta0_fed  <- model0_fed$coefficients[covariate_names, "Estimate"]

cat("=== Federated GLM coefficients ===\n")
cat(sprintf("  alpha1=%.8f  beta1=%s\n", alpha1_fed,
            paste(round(beta1_fed, 8), collapse = ", ")))
cat(sprintf("  alpha0=%.8f  beta0=%s\n\n", alpha0_fed,
            paste(round(beta0_fed, 8), collapse = ", ")))

# ============================================================================
# Step 3: manual pooled variance using the SAME federated coefficients
# applied to the merged dataset — this is the ground truth for comparison.
# ============================================================================

X_aug <- cbind(1, as.matrix(dat[, covariate_names]))
mu1   <- as.numeric(plogis(X_aug %*% c(alpha1_fed, beta1_fed)))
mu0   <- as.numeric(plogis(X_aug %*% c(alpha0_fed, beta0_fed)))
delta <- mu1 - mu0
ate_manual <- mean(delta)

v1   <- mu1 * (1 - mu1)
v0   <- mu0 * (1 - mu0)
idx1 <- dat$W == 1L

X_aug1 <- X_aug[ idx1, ];  X_aug0 <- X_aug[!idx1, ]
H1_manual <- crossprod(X_aug1, v1[ idx1] * X_aug1)
H0_manual <- crossprod(X_aug0, v0[!idx1] * X_aug0)

g1_manual <-  colMeans(v1 * X_aug)
g0_manual <- -colMeans(v0 * X_aug)

var_param_manual <- as.numeric(
  crossprod(g1_manual, solve(H1_manual, g1_manual)) +
  crossprod(g0_manual, solve(H0_manual, g0_manual))
)
var_samp_manual <- var(delta) / n
var_manual      <- var_param_manual + var_samp_manual

cat("=== Manual pooled (same federated coefficients) ===\n")
cat(sprintf("  ATE:        %.10f\n",   ate_manual))
cat(sprintf("  Var(param): %.6e\n",   var_param_manual))
cat(sprintf("  Var(samp):  %.6e\n",   var_samp_manual))
cat(sprintf("  Total var:  %.10f\n\n", var_manual))

# ============================================================================
# Step 4: comparison
# ============================================================================

ate_diff <- abs(result_fed$ate      - ate_manual)
var_diff <- abs(result_fed$variance - var_manual)

cat("=== Federated vs manual pooled (same coefficients) ===\n")
cat(sprintf("  |ATE diff|:      %.2e  (expect <= 1e-10)\n", ate_diff))
cat(sprintf("  |Variance diff|: %.2e  (expect <= 1e-10)\n", var_diff))

stopifnot(ate_diff < 1e-10, var_diff < 1e-10)
cat("\nAll checks passed: federated == pooled at machine precision.\n")

# ============================================================================
# For reference: IRLS convergence gap between ds.glm() and R's glm()
# ============================================================================

model1_r <- glm(Y ~ X1 + X2, family = binomial,
                data = dat[dat$W == 1, ])
model0_r <- glm(Y ~ X1 + X2, family = binomial,
                data = dat[dat$W == 0, ])

coef_diff1 <- max(abs(c(alpha1_fed, beta1_fed) -
                      c(coef(model1_r)["(Intercept)"], coef(model1_r)[covariate_names])))
coef_diff0 <- max(abs(c(alpha0_fed, beta0_fed) -
                      c(coef(model0_r)["(Intercept)"], coef(model0_r)[covariate_names])))
cat(sprintf("\n[Info] IRLS convergence gap (ds.glm vs glm): arm1=%.2e, arm0=%.2e\n",
            coef_diff1, coef_diff0))
cat("       (Non-zero because federated IRLS converges to ~3e-9 tolerance)\n")

DSI::datashield.logout(conns)
