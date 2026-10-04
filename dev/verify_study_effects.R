# Verify that ds.EstimateAte(..., study_effects = TRUE) gives the same result
# as the equivalent pooled OLS on the augmented design matrix.
#
# Setting: K=3 servers with heterogeneous treatment probabilities and
# site-specific intercept shifts.  The unadjusted estimator is biased here;
# the adjusted estimator (study_effects = TRUE) should match the pooled OLS
# on the design matrix augmented with K-1 study dummy columns.

library(DSLite)
library(DSI)
devtools::load_all("../dsFedAte")
devtools::load_all(".")

# ── 1. Generate data ───────────────────────────────────────────────────────────

set.seed(42)

K        <- 3L
p_k      <- c(0.8, 0.5, 0.2)   # heterogeneous treatment probabilities
h_k      <- c(3, 0, -3)        # site-specific intercept shifts (same for both arms)
true_ate <- 0.5

alpha0 <- 1.0;  alpha1 <- alpha0 + true_ate
beta0  <- c(0.3, -0.2)
beta1  <- c(0.5,  0.1)

make_site <- function(n, p, h) {
  df   <- data.frame(X1 = rnorm(n), X2 = rnorm(n), W = rbinom(n, 1, p))
  df$Y <- h + ifelse(df$W == 1,
    alpha1 + df$X1 * beta1[1] + df$X2 * beta1[2],
    alpha0 + df$X1 * beta0[1] + df$X2 * beta0[2]
  ) + rnorm(n, sd = 0.8)
  df
}

site1 <- make_site(150, p_k[1], h_k[1])
site2 <- make_site(150, p_k[2], h_k[2])
site3 <- make_site(150, p_k[3], h_k[3])

sites     <- list(site1, site2, site3)
pooled    <- do.call(rbind, sites)

# ── 2. Federated computation via DSLite ───────────────────────────────────────

srv <- lapply(sites, function(df) {
  s <- DSLite::newDSLiteServer(tables = list(study_data = df))
  s$config(defaultDSConfiguration(include = c("dsBase", "dsFedAte")))
  s$aggregateMethod("localAteDS",                        "localAteDS")
  s$aggregateMethod("listDisclosureSettingsDS",           "listDisclosureSettingsDS")
  s$aggregateMethod("checkPermissivePrivacyControlLevel", "checkPermissivePrivacyControlLevel")
  s$assignMethod("addStudyIdDS", "addStudyIdDS")
  s
})

# DSLite requires server objects reachable by name from .GlobalEnv
for (k in seq_len(K)) assign(paste0("srv", k), srv[[k]], envir = .GlobalEnv)

options(datashield.env            = .GlobalEnv,
        datashield.privacyControlLevel = "banana",
        nfilter.subset             = 3,
        nfilter.glm                = 0.33)

builder <- DSI::newDSLoginBuilder()
for (k in seq_len(K))
  builder$append(server = paste0("site", k),
                 url    = paste0("srv", k),
                 table  = "study_data",
                 driver = "DSLiteDriver")

conns <- DSI::datashield.login(logins = builder$build(), assign = FALSE)
for (k in seq_len(K))
  DSI::datashield.assign.table(conns[paste0("site", k)],
                               symbol = "study_data",
                               table  = "study_data")

fed <- suppressWarnings(
  ds.EstimateAte(
    datasources     = conns,
    base_name       = "study_data",
    outcome_name    = "Y",
    treatment_name  = "W",
    covariate_names = c("X1", "X2"),
    study_effects   = TRUE,
    heterogeneous   = FALSE
  )
)

DSI::datashield.logout(conns)

# ── 3. Pooled (centralised) computation on augmented design matrix ─────────────
#
# Mirrors exactly what ds.EstimateAte does server-side:
#   - Append study_id2, study_id3 to each site's data (addStudyIdDS logic)
#   - Fit arm-specific OLS on augmented design matrix
#   - Compute ATE and variance with d = length(original covariates) + (K-1)

add_study_dummies <- function(df, k, K) {
  for (j in seq(2L, K))
    df[[paste0("study_id", j)]] <- as.integer(k == j)
  df
}

sites_aug <- lapply(seq_len(K), function(k) add_study_dummies(sites[[k]], k, K))
pooled_aug <- do.call(rbind, sites_aug)

cov_names_aug <- c("X1", "X2", paste0("study_id", seq(2L, K)))  # X1, X2, study_id2, study_id3
d_aug         <- length(cov_names_aug)                            # 4 for K=3
formula_aug   <- paste0("Y ~ ", paste(cov_names_aug, collapse = " + "))

treated_aug <- pooled_aug[pooled_aug$W == 1, ]
control_aug <- pooled_aug[pooled_aug$W == 0, ]

fit1 <- lm(formula_aug, data = treated_aug)
fit0 <- lm(formula_aug, data = control_aug)

cen_alpha1 <- coef(fit1)[["(Intercept)"]]
cen_alpha0 <- coef(fit0)[["(Intercept)"]]
cen_beta1  <- coef(fit1)[cov_names_aug]
cen_beta0  <- coef(fit0)[cov_names_aug]
cen_gamma  <- cen_beta1 - cen_beta0

x_bar_aug <- colMeans(pooled_aug[, cov_names_aug])
cen_ate   <- as.numeric((cen_alpha1 - cen_alpha0) + sum(x_bar_aug * cen_gamma))

# Variance — same formula as .compute_variance(), with d = d_aug
n     <- nrow(pooled_aug)
n_1   <- nrow(treated_aug)
p_hat <- n_1 / n

rss_1  <- sum(residuals(fit1)^2)
rss_0  <- sum(residuals(fit0)^2)
sigma2 <- (rss_1 + rss_0) / (n - 2L * (d_aug + 1L))

# second_term: γ̂ᵀ Σ̂ γ̂ / n  (Σ̂ = sample cov of augmented X, uses n-1 denom)
Sigma_hat   <- cov(pooled_aug[, cov_names_aug])
first_term  <- sigma2 / (n * p_hat * (1 - p_hat))
second_term <- as.numeric(t(cen_gamma) %*% Sigma_hat %*% cen_gamma) / n
cen_var     <- first_term + second_term

# ── 4. Comparison ─────────────────────────────────────────────────────────────

cat("\n── ATE ───────────────────────────────────────────────\n")
cat(sprintf("  Federated : %.8f\n", fed$ate))
cat(sprintf("  Pooled    : %.8f\n", cen_ate))
cat(sprintf("  Diff      : %.2e\n", abs(fed$ate - cen_ate)))

cat("\n── Variance ──────────────────────────────────────────\n")
cat(sprintf("  Federated : %.8f\n", fed$variance))
cat(sprintf("  Pooled    : %.8f\n", cen_var))
cat(sprintf("  Diff      : %.2e\n", abs(fed$variance - cen_var)))

cat("\n── 95%% CI ────────────────────────────────────────────\n")
z <- 1.96 * sqrt(cen_var)
cat(sprintf("  Federated : [%.6f, %.6f]\n", fed$ci_lower, fed$ci_upper))
cat(sprintf("  Pooled    : [%.6f, %.6f]\n", cen_ate - z, cen_ate + z))

cat("\n── Sanity: unadjusted ATE (should be biased) ─────────\n")
conns2 <- DSI::datashield.login(logins = builder$build(), assign = FALSE)
for (k in seq_len(K))
  DSI::datashield.assign.table(conns2[paste0("site", k)],
                               symbol = "study_data", table = "study_data")
fed_unadj <- suppressWarnings(
  ds.EstimateAte(
    datasources     = conns2,
    base_name       = "study_data",
    outcome_name    = "Y",
    treatment_name  = "W",
    covariate_names = c("X1", "X2"),
    study_effects   = FALSE,
    heterogeneous   = FALSE
  )
)
DSI::datashield.logout(conns2)
cat(sprintf("  Unadjusted ATE : %.8f  (true ATE = %.1f)\n", fed_unadj$ate, true_ate))
cat(sprintf("  Adjusted ATE   : %.8f  (true ATE = %.1f)\n", fed$ate,       true_ate))

cat("\n── Components ────────────────────────────────────────\n")
cat(sprintf("  sigma2   (pooled) : %.8f\n", sigma2))
cat(sprintf("  1st term (pooled) : %.8f\n", first_term))
cat(sprintf("  2nd term (pooled) : %.8f\n", second_term))
cat(sprintf("  d (aug)           : %d  (= %d covariates + %d study dummies)\n",
            d_aug, 2L, K - 1L))
