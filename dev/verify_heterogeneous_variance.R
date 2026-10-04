# Verify that ds.EstimateAte(heterogeneous = TRUE) matches the centralised
# plug-in variance (Proposition 2 of the companion variance paper).
#
# Uses K = 3 servers with different covariate means (H -> X active) AND
# different treatment probabilities (H -> W active), so the leverage
# correction delta^(w) > 0 and the homogeneous closed form is wrong.
# The heterogeneous estimator should match the centralised formula at
# machine precision.

library(DSLite)
library(DSI)
devtools::load_all("../dsFedAte")
devtools::load_all(".")

# ── 1. Generate heterogeneous data ────────────────────────────────────────────
#
# Site 1: mu = (2, 0),  p = 0.8   (high treatment, X1 shifted right)
# Site 2: mu = (0, 0),  p = 0.5   (balanced)
# Site 3: mu = (-2, 1), p = 0.2   (low treatment, X1 shifted left, X2 shifted)
#
# Shared outcome model: Y = alpha_w + X1*beta_w[1] + X2*beta_w[2] + eps
# No study effects (H -> Y inactive), so the unadjusted estimator is unbiased.

set.seed(42)

make_site <- function(n, mu1, mu2, p,
                      alpha0 = 1.0, alpha1 = 1.5,
                      beta0  = c(0.3, -0.2), beta1 = c(0.5, 0.1)) {
  X1 <- rnorm(n, mean = mu1)
  X2 <- rnorm(n, mean = mu2)
  W  <- rbinom(n, 1, p)
  Y  <- ifelse(W == 1,
    alpha1 + X1 * beta1[1] + X2 * beta1[2],
    alpha0 + X1 * beta0[1] + X2 * beta0[2]
  ) + rnorm(n, sd = 0.8)
  data.frame(Y = Y, W = W, X1 = X1, X2 = X2)
}

site1 <- make_site(150, mu1 =  2, mu2 =  0, p = 0.8)
site2 <- make_site(120, mu1 =  0, mu2 =  0, p = 0.5)
site3 <- make_site(130, mu1 = -2, mu2 =  1, p = 0.2)

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
options(datashield.env                 = .GlobalEnv,
        datashield.privacyControlLevel = "banana",
        nfilter.subset                 = 3,
        nfilter.glm                    = 0.33)

builder <- DSI::newDSLoginBuilder()
builder$append(server = "site1", url = "srv1", table = "study_data", driver = "DSLiteDriver")
builder$append(server = "site2", url = "srv2", table = "study_data", driver = "DSLiteDriver")
builder$append(server = "site3", url = "srv3", table = "study_data", driver = "DSLiteDriver")
conns <- DSI::datashield.login(logins = builder$build(), assign = FALSE)
for (s in c("site1", "site2", "site3")) {
  DSI::datashield.assign.table(conns[s], symbol = "study_data", table = "study_data")
}

fed_het <- suppressWarnings(
  ds.EstimateAte(datasources     = conns,
                 base_name       = "study_data",
                 outcome_name    = "Y",
                 treatment_name  = "W",
                 covariate_names = c("X1", "X2"),
                 heterogeneous   = TRUE)
)

fed_hom <- suppressWarnings(
  ds.EstimateAte(datasources     = conns,
                 base_name       = "study_data",
                 outcome_name    = "Y",
                 treatment_name  = "W",
                 covariate_names = c("X1", "X2"),
                 heterogeneous   = FALSE)
)

DSI::datashield.logout(conns)

# ── 3. Centralised computation (Proposition 2) ────────────────────────────────
#
# Term (I) via eq. (15): sigma2_hat * X_bar' ((G^1)^{-1} + (G^0)^{-1}) X_bar'
# Term (II) as pooled sample variance of the projected index X_i' Delta-theta.

pooled  <- rbind(site1, site2, site3)
n       <- nrow(pooled)
d       <- 2L

treated <- pooled[pooled$W == 1, ]
control <- pooled[pooled$W == 0, ]

fit1    <- lm(Y ~ X1 + X2, data = treated)
fit0    <- lm(Y ~ X1 + X2, data = control)

alpha1  <- coef(fit1)[["(Intercept)"]]
alpha0  <- coef(fit0)[["(Intercept)"]]
beta1   <- coef(fit1)[c("X1", "X2")]
beta0   <- coef(fit0)[c("X1", "X2")]
gamma   <- beta1 - beta0

# ATE
X_bar    <- colMeans(pooled[, c("X1", "X2")])
cen_ate  <- (alpha1 - alpha0) + sum(X_bar * gamma)

# sigma2 with pooled df
rss_1   <- sum(residuals(fit1)^2)
rss_0   <- sum(residuals(fit0)^2)
sigma2  <- (rss_1 + rss_0) / (n - 2L * (d + 1L))

# Term (I): centralised plug-in via Gram inverses
X1_aug  <- model.matrix(~ X1 + X2, data = treated)
X0_aug  <- model.matrix(~ X1 + X2, data = control)
G1_inv  <- solve(t(X1_aug) %*% X1_aug)
G0_inv  <- solve(t(X0_aug) %*% X0_aug)
X_prime <- c(1, X_bar)
term_I  <- sigma2 * as.numeric(t(X_prime) %*% (G1_inv + G0_inv) %*% X_prime)

# Term (II): pooled sample variance of the projected index
proj    <- as.numeric(as.matrix(pooled[, c("X1", "X2")]) %*% gamma)
term_II <- var(proj) / n    # var() uses n-1 denominator = Delta-beta' Sigma_hat Delta-beta / n

cen_var <- term_I + term_II

# ── 4. Comparison ─────────────────────────────────────────────────────────────

cat("\n── ATE ───────────────────────────────────────────────────────────────\n")
cat(sprintf("  Federated (het) : %.10f\n", fed_het$ate))
cat(sprintf("  Centralized     : %.10f\n", cen_ate))
cat(sprintf("  Diff            : %.2e\n",  abs(fed_het$ate - cen_ate)))

cat("\n── Variance ──────────────────────────────────────────────────────────\n")
cat(sprintf("  Federated (het) : %.10f\n", fed_het$variance))
cat(sprintf("  Federated (hom) : %.10f\n", fed_hom$variance))
cat(sprintf("  Centralized     : %.10f\n", cen_var))
cat(sprintf("  Het - Central   : %.2e\n",  abs(fed_het$variance - cen_var)))
cat(sprintf("  Hom - Central   : %.2e  (expect >> 0 due to leverage)\n",
            abs(fed_hom$variance - cen_var)))

cat("\n── Term decomposition ────────────────────────────────────────────────\n")
n_1     <- nrow(treated)
p_hat   <- n_1 / n
hom_I   <- sigma2 / (n * p_hat * (1 - p_hat))
cat(sprintf("  Term I  (centralised)  : %.10f\n", term_I))
cat(sprintf("  Term I  (homog form)   : %.10f\n", hom_I))
cat(sprintf("  Term II                : %.10f\n", term_II))
cat(sprintf("  Leverage inflation     : %.2e\n",  term_I - hom_I))

cat("\n── Pass / fail ───────────────────────────────────────────────────────\n")
tol <- 1e-10
ok_ate <- abs(fed_het$ate - cen_ate)       < tol
ok_var <- abs(fed_het$variance - cen_var)  < tol
cat(sprintf("  ATE match at %.0e : %s\n", tol, if (ok_ate) "PASS" else "FAIL"))
cat(sprintf("  Var match at %.0e : %s\n", tol, if (ok_var) "PASS" else "FAIL"))
