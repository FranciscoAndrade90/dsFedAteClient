# Verify that the federated variance equals the centralized pooled-OLS variance.
#
# Uses K=3 servers with unequal sizes.  Compares ds.EstimateAte() output to
# the equivalent computation on the merged dataset.

library(DSLite)
library(DSI)
devtools::load_all("../dsFedAte")
devtools::load_all(".")

# ── 1. Generate data across three servers ─────────────────────────────────────

set.seed(123)

make_site <- function(n, alpha0 = 1, alpha1 = 1.5,
                      beta0 = c(0.3, -0.2), beta1 = c(0.5, 0.1)) {
  df      <- data.frame(X1 = rnorm(n), X2 = rnorm(n), W = rbinom(n, 1, 0.5))
  df$Y    <- ifelse(df$W == 1,
    alpha1 + df$X1 * beta1[1] + df$X2 * beta1[2],
    alpha0 + df$X1 * beta0[1] + df$X2 * beta0[2]
  ) + rnorm(n, sd = 0.8)
  df
}

site1 <- make_site(120)
site2 <- make_site(80)
site3 <- make_site(200)

# ── 2. Federated computation via DSLite ───────────────────────────────────────

srv1 <- DSLite::newDSLiteServer(tables = list(study_data = site1))
srv2 <- DSLite::newDSLiteServer(tables = list(study_data = site2))
srv3 <- DSLite::newDSLiteServer(tables = list(study_data = site3))

for (srv in list(srv1, srv2, srv3)) {
  srv$config(defaultDSConfiguration(include = c("dsBase", "dsFedAte")))
  srv$aggregateMethod("localAteDS",                        "localAteDS")
  srv$aggregateMethod("listDisclosureSettingsDS",           "listDisclosureSettingsDS")
  srv$aggregateMethod("checkPermissivePrivacyControlLevel", "checkPermissivePrivacyControlLevel")
}

assign("srv1", srv1, envir = .GlobalEnv)
assign("srv2", srv2, envir = .GlobalEnv)
assign("srv3", srv3, envir = .GlobalEnv)
options(datashield.env        = .GlobalEnv,
        datashield.privacyControlLevel = "banana",
        nfilter.subset        = 3,
        nfilter.glm           = 0.33)

builder <- DSI::newDSLoginBuilder()
builder$append(server = "site1", url = "srv1", table = "study_data", driver = "DSLiteDriver")
builder$append(server = "site2", url = "srv2", table = "study_data", driver = "DSLiteDriver")
builder$append(server = "site3", url = "srv3", table = "study_data", driver = "DSLiteDriver")
conns <- DSI::datashield.login(logins = builder$build(), assign = FALSE)
for (s in c("site1", "site2", "site3")) {
  DSI::datashield.assign.table(conns[s], symbol = "study_data", table = "study_data")
}

fed <- suppressWarnings(
  ds.EstimateAte(datasources     = conns,
                 base_name       = "study_data",
                 outcome_name    = "Y",
                 treatment_name  = "W",
                 covariate_names = c("X1", "X2"),
                 heterogeneous   = FALSE)
)

DSI::datashield.logout(conns)

# ── 3. Centralized computation on pooled data ─────────────────────────────────

pooled <- rbind(site1, site2, site3)
n      <- nrow(pooled)
d      <- 2L

treated <- pooled[pooled$W == 1, ]
control <- pooled[pooled$W == 0, ]

fit1   <- lm(Y ~ X1 + X2, data = treated)
fit0   <- lm(Y ~ X1 + X2, data = control)

alpha1 <- coef(fit1)[["(Intercept)"]]
alpha0 <- coef(fit0)[["(Intercept)"]]
beta1  <- coef(fit1)[c("X1", "X2")]
beta0  <- coef(fit0)[c("X1", "X2")]
gamma  <- beta1 - beta0

x_bar <- colMeans(pooled[, c("X1", "X2")])
cen_ate <- (alpha1 - alpha0) + sum(x_bar * gamma)

# Variance
n_1    <- nrow(treated)
p_hat  <- n_1 / n
rss_1  <- sum(residuals(fit1)^2)
rss_0  <- sum(residuals(fit0)^2)
sigma2 <- (rss_1 + rss_0) / (n - 2L * (d + 1L))

Sigma_hat   <- cov(pooled[, c("X1", "X2")])   # uses n-1 denominator
first_term  <- sigma2 / (n * p_hat * (1 - p_hat))
second_term <- as.numeric(t(gamma) %*% Sigma_hat %*% gamma) / n
cen_var     <- first_term + second_term

# ── 4. Comparison ─────────────────────────────────────────────────────────────

cat("\n── ATE ───────────────────────────────────────────────\n")
cat(sprintf("  Federated : %.8f\n", fed$ate))
cat(sprintf("  Central   : %.8f\n", cen_ate))
cat(sprintf("  Diff      : %.2e\n",  abs(fed$ate - cen_ate)))

cat("\n── Variance ──────────────────────────────────────────\n")
cat(sprintf("  Federated : %.8f\n", fed$variance))
cat(sprintf("  Central   : %.8f\n", cen_var))
cat(sprintf("  Diff      : %.2e\n",  abs(fed$variance - cen_var)))

cat("\n── 95%% CI ────────────────────────────────────────────\n")
z <- 1.96 * sqrt(cen_var)
cat(sprintf("  Federated : [%.6f, %.6f]\n", fed$ci_lower, fed$ci_upper))
cat(sprintf("  Central   : [%.6f, %.6f]\n", cen_ate - z, cen_ate + z))

cat("\n── Components ────────────────────────────────────────\n")
cat(sprintf("  sigma2_hat (central)  : %.8f\n", sigma2))
cat(sprintf("  first_term (central)  : %.8f\n", first_term))
cat(sprintf("  second_term (central) : %.8f\n", second_term))
