library(devtools)
devtools::load_all(".")
devtools::load_all("../dsFedAte")

# ============================================================================
# DGP (same as verify_binary_ate.R)
# logit P(Y=1 | X, W=w) = alpha_w + beta'X
# beta identical in both arms so the true ATE is entirely due to alpha1 - alpha0
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
dat         <- do.call(rbind, data_list)
n           <- nrow(dat)
cov_names   <- c("X1", "X2")

# ============================================================================
# Step 1: arm-specific logistic regressions (pooled)
# ============================================================================

dat_1 <- dat[dat$W == 1, ]
dat_0 <- dat[dat$W == 0, ]

model1 <- glm(Y ~ X1 + X2, family = binomial, data = dat_1)
model0 <- glm(Y ~ X1 + X2, family = binomial, data = dat_0)

cat("=== Arm-specific logistic regression coefficients ===\n")
cat("Treated arm:  alpha1 =", coef(model1)[1],
    " beta1 =", coef(model1)[-1], "\n")
cat("Control arm:  alpha0 =", coef(model0)[1],
    " beta0 =", coef(model0)[-1], "\n\n")

# ============================================================================
# Step 2: G-computation ATE
# ============================================================================

X_aug  <- cbind(1, as.matrix(dat[, cov_names]))   # n x (d+1), all observations

mu1    <- as.numeric(plogis(X_aug %*% coef(model1)))   # mu_1(X_i) for all i
mu0    <- as.numeric(plogis(X_aug %*% coef(model0)))   # mu_0(X_i) for all i

ate    <- mean(mu1 - mu0)
cat(sprintf("G-computation ATE: %.6f\n\n", ate))

# ============================================================================
# Step 3: Hessian matrices
#
# H_w = sum_{i: W_i=w} mu_w(X_i)(1-mu_w(X_i)) * x~_i x~_i'
#
# This is the observed Fisher information for arm w.
# vcov(model_w) returns H_w^{-1} directly for binomial GLM.
# We verify the manual computation against vcov().
# ============================================================================

# Augmented design matrices for each arm
X_aug1 <- cbind(1, as.matrix(dat_1[, cov_names]))   # n_1 x (d+1)
X_aug0 <- cbind(1, as.matrix(dat_0[, cov_names]))   # n_0 x (d+1)

# IRLS weights: mu_w(X_i)(1 - mu_w(X_i)) for arm-w observations
w1 <- as.numeric(fitted(model1) * (1 - fitted(model1)))
w0 <- as.numeric(fitted(model0) * (1 - fitted(model0)))

# Hessians (observed Fisher information matrices, positive definite)
H1 <- t(X_aug1) %*% (w1 * X_aug1)   # (d+1) x (d+1)
H0 <- t(X_aug0) %*% (w0 * X_aug0)

cat("=== Hessian H1 (treated arm) ===\n")
print(round(H1, 4))
cat("\n=== Hessian H0 (control arm) ===\n")
print(round(H0, 4))

# Verify: H_w^{-1} should equal vcov(model_w)
cat("\nMax discrepancy H1^{-1} vs vcov(model1):",
    max(abs(solve(H1) - vcov(model1))), "\n")
cat("Max discrepancy H0^{-1} vs vcov(model0):",
    max(abs(solve(H0) - vcov(model0))), "\n\n")

# ============================================================================
# Step 4: Gradient of ATE w.r.t. model parameters
#
# g_w = (1/n) sum_{i=1}^{n} mu_w(X_i)(1-mu_w(X_i)) * x~_i
#
# Note: sum is over ALL observations (both arms), not just arm w.
# ============================================================================

v1 <- mu1 * (1 - mu1)   # mu_1(X_i)(1 - mu_1(X_i)) for all i
v0 <- mu0 * (1 - mu0)   # mu_0(X_i)(1 - mu_0(X_i)) for all i

g1 <-  colMeans(v1 * X_aug)   # (d+1) vector,  sign: d(tau)/d(theta1)
g0 <- -colMeans(v0 * X_aug)   # (d+1) vector,  sign: d(tau)/d(theta0)

cat("=== Gradient g1 (wrt treated-arm params) ===\n")
print(round(g1, 6))
cat("\n=== Gradient g0 (wrt control-arm params) ===\n")
print(round(g0, 6))
cat("\n")

# ============================================================================
# Step 5: Delta-method variance (parameter estimation uncertainty)
#
# Var_param = g1' H1^{-1} g1 + g0' H0^{-1} g0
# ============================================================================

H1_inv <- solve(H1)
H0_inv <- solve(H0)

var_param <- as.numeric(t(g1) %*% H1_inv %*% g1 +
                        t(g0) %*% H0_inv %*% g0)

# ============================================================================
# Step 6: Sampling variance (finite covariate sample)
#
# Var_samp = Var(mu_1(X_i) - mu_0(X_i)) / n
# ============================================================================

var_samp <- var(mu1 - mu0) / n

# ============================================================================
# Step 7: Total variance and 95% CI
# ============================================================================

variance <- var_param + var_samp
se       <- sqrt(variance)
ci_lower <- ate - 1.96 * se
ci_upper <- ate + 1.96 * se

cat("=== Variance decomposition ===\n")
cat(sprintf("  Parameter estimation (delta method): %.6e\n", var_param))
cat(sprintf("  Sampling (finite covariate sample):  %.6e\n", var_samp))
cat(sprintf("  Total:                               %.6e\n", variance))
cat(sprintf("  SE:                                  %.6f\n\n", se))

cat("=== 95%% confidence interval ===\n")
cat(sprintf("  ATE:      %.6f\n", ate))
cat(sprintf("  CI lower: %.6f\n", ci_lower))
cat(sprintf("  CI upper: %.6f\n", ci_upper))
