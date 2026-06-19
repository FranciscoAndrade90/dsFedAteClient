# Compare three variance estimators for the federated ATE:
#
#   (A) Federated     — pooled sigma2 = (RSS_1 + RSS_0) / (n - 2(d+1))
#   (B) Central pooled — same pooled sigma2, computed on merged data
#   (C) Central natural — arm-specific sigma2 embedded in vcov(lm(...))
#
# (A) and (B) should agree at machine precision.
# (C) is asymptotically equivalent but differs in finite samples.
#
# Tested under two scenarios:
#   S1 — balanced arms, homogeneous sites  (small expected difference)
#   S2 — unbalanced arms, heterogeneous sites (larger expected difference)

library(DSLite)
library(DSI)
devtools::load_all("../dsRctFed")
devtools::load_all(".")

# ── Helper: build DSLite connections from a list of data frames ───────────────

make_conns <- function(sites) {
  srvs <- lapply(sites, function(df)
    DSLite::newDSLiteServer(tables = list(study_data = df)))
  for (srv in srvs) {
    srv$config(defaultDSConfiguration(include = c("dsBase", "dsRctFed")))
    srv$aggregateMethod("localAteDS",                        "localAteDS")
    srv$aggregateMethod("listDisclosureSettingsDS",           "listDisclosureSettingsDS")
    srv$aggregateMethod("checkPermissivePrivacyControlLevel", "checkPermissivePrivacyControlLevel")
  }
  for (j in seq_along(srvs))
    assign(paste0("srv_pn_", j), srvs[[j]], envir = .GlobalEnv)
  options(datashield.env                 = .GlobalEnv,
          datashield.privacyControlLevel = "banana",
          nfilter.subset                 = 3,
          nfilter.glm                    = 0.33)
  builder <- DSI::newDSLoginBuilder()
  for (j in seq_along(srvs))
    builder$append(server  = paste0("site", j),
                   url     = paste0("srv_pn_", j),
                   table   = "study_data",
                   driver  = "DSLiteDriver")
  conns <- DSI::datashield.login(logins = builder$build(), assign = FALSE)
  for (j in seq_along(srvs))
    DSI::datashield.assign.table(conns[paste0("site", j)],
                                 symbol = "study_data", table = "study_data")
  conns
}

# ── Helper: all three variance estimates on one dataset ───────────────────────

compare_estimators <- function(sites, heterogeneous, label) {

  cat(sprintf("\n%s\n%s\n", label, strrep("─", nchar(label))))

  # (A) Federated
  conns <- make_conns(sites)
  fed <- suppressWarnings(
    ds.EstimateAte(datasources     = conns,
                   base_name       = "study_data",
                   outcome_name    = "Y",
                   treatment_name  = "W",
                   covariate_names = c("X1", "X2"),
                   heterogeneous   = heterogeneous)
  )
  DSI::datashield.logout(conns)

  # Centralised ingredients
  pooled  <- do.call(rbind, sites)
  n       <- nrow(pooled)
  d       <- 2L
  treated <- pooled[pooled$W == 1, ]
  control <- pooled[pooled$W == 0, ]

  fit1    <- lm(Y ~ X1 + X2, data = treated)
  fit0    <- lm(Y ~ X1 + X2, data = control)

  alpha1  <- coef(fit1)[["(Intercept)"]]; beta1 <- coef(fit1)[c("X1","X2")]
  alpha0  <- coef(fit0)[["(Intercept)"]]; beta0 <- coef(fit0)[c("X1","X2")]
  gamma   <- beta1 - beta0
  x_bar   <- colMeans(pooled[, c("X1","X2")])
  x_prime <- c(1, x_bar)

  # Sampling term (common to B and C)
  proj      <- as.numeric(as.matrix(pooled[, c("X1","X2")]) %*% gamma)
  proj_mean <- mean(proj)
  if (heterogeneous) {
    # delta-method sampling term: tss / (n*(n-1))
    delta_i   <- (alpha1 - alpha0) + as.numeric(as.matrix(pooled[, c("X1","X2")]) %*% gamma)
    ate_hat   <- mean(delta_i)
    tss       <- sum((delta_i - ate_hat)^2)
    samp_term <- tss / (n * (n - 1L))
  } else {
    samp_term <- var(proj) / n     # var() uses n-1: = Delta-beta' Sigma_hat Delta-beta / n
  }

  # (B) Central pooled sigma2
  rss_1   <- sum(residuals(fit1)^2)
  rss_0   <- sum(residuals(fit0)^2)
  sigma2_pooled <- (rss_1 + rss_0) / (n - 2L * (d + 1L))

  if (heterogeneous) {
    X1_aug  <- model.matrix(~ X1 + X2, data = treated)
    X0_aug  <- model.matrix(~ X1 + X2, data = control)
    G1_inv  <- solve(t(X1_aug) %*% X1_aug)
    G0_inv  <- solve(t(X0_aug) %*% X0_aug)
    param_B <- sigma2_pooled *
               (as.numeric(t(x_prime) %*% G1_inv %*% x_prime) +
                as.numeric(t(x_prime) %*% G0_inv %*% x_prime))
  } else {
    n_1     <- nrow(treated)
    p_hat   <- n_1 / n
    param_B <- sigma2_pooled / (n * p_hat * (1 - p_hat))
  }
  var_B <- param_B + samp_term

  # (C) Central natural: arm-specific sigma2, same structural formula as (A/B).
  #     het=FALSE: simplified formula sigma2_w/(n*p*(1-p)) — no leverage term,
  #                matching the assumption made by the federated simplified path.
  #     het=TRUE:  full x_prime' vcov(fit_w) x_prime, which includes leverage
  #                automatically via the Gram inverse (Lemma 2 of the paper).
  rss_1        <- sum(residuals(fit1)^2)
  rss_0        <- sum(residuals(fit0)^2)
  sigma2_1     <- rss_1 / (nrow(treated) - (d + 1L))
  sigma2_0     <- rss_0 / (nrow(control) - (d + 1L))
  if (heterogeneous) {
    param_C <- as.numeric(t(x_prime) %*% vcov(fit1) %*% x_prime) +
               as.numeric(t(x_prime) %*% vcov(fit0) %*% x_prime)
  } else {
    n_1     <- nrow(treated)
    p_hat   <- n_1 / n
    # arm-specific analog of sigma2/(n*p*(1-p)) = sigma2*(1/n_1 + 1/n_0)
    param_C <- sigma2_1 / n_1 + sigma2_0 / nrow(control)
  }
  var_C <- param_C + samp_term

  cat(sprintf("  (A) Federated (pooled sigma2) : %.10f\n", fed$variance))
  cat(sprintf("  (B) Central   (pooled sigma2) : %.10f\n", var_B))
  cat(sprintf("  (C) Central   (natural vcov)  : %.10f\n", var_C))
  cat(sprintf("  |A - B|                       : %.2e\n", abs(fed$variance - var_B)))
  cat(sprintf("  |A - C|                       : %.2e\n", abs(fed$variance - var_C)))
  cat(sprintf("  |B - C|                       : %.2e\n", abs(var_B - var_C)))
  cat(sprintf("  sigma2_pooled                 : %.6f\n", sigma2_pooled))
  cat(sprintf("  sigma2 arm-1                  : %.6f\n", rss_1 / (nrow(treated) - (d+1L))))
  cat(sprintf("  sigma2 arm-0                  : %.6f\n", rss_0 / (nrow(control) - (d+1L))))
}

# ── Scenario 1: balanced arms, homogeneous sites ──────────────────────────────

set.seed(1)
make_balanced <- function(n) {
  df   <- data.frame(X1 = rnorm(n), X2 = rnorm(n), W = rbinom(n, 1, 0.5))
  df$Y <- ifelse(df$W == 1, 1.5 + 0.5*df$X1 + 0.1*df$X2,
                             1.0 + 0.3*df$X1 - 0.2*df$X2) + rnorm(n, sd = 0.8)
  df
}
sites_s1 <- list(make_balanced(120), make_balanced(80), make_balanced(200))

compare_estimators(sites_s1, heterogeneous = FALSE,
  "Scenario 1 — balanced arms, homogeneous sites, heterogeneous=FALSE")
compare_estimators(sites_s1, heterogeneous = TRUE,
  "Scenario 1 — balanced arms, homogeneous sites, heterogeneous=TRUE")

# ── Scenario 2: unbalanced arms, heterogeneous sites ─────────────────────────

set.seed(2)
make_site <- function(n, p, mu1, mu2) {
  X1 <- rnorm(n, mu1); X2 <- rnorm(n, mu2)
  W  <- rbinom(n, 1, p)
  Y  <- ifelse(W == 1, 1.5 + 0.5*X1 + 0.1*X2,
                        1.0 + 0.3*X1 - 0.2*X2) + rnorm(n, sd = 0.8)
  data.frame(Y=Y, W=W, X1=X1, X2=X2)
}
sites_s2 <- list(make_site(150, p=0.8, mu1= 2, mu2= 0),
                 make_site(120, p=0.5, mu1= 0, mu2= 0),
                 make_site(130, p=0.2, mu1=-2, mu2= 1))

compare_estimators(sites_s2, heterogeneous = FALSE,
  "Scenario 2 — unbalanced arms, heterogeneous sites, heterogeneous=FALSE")
compare_estimators(sites_s2, heterogeneous = TRUE,
  "Scenario 2 — unbalanced arms, heterogeneous sites, heterogeneous=TRUE")
