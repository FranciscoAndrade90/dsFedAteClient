devtools::load_all("/home/francisco/DATASHIELD/dsRctFedClient")
devtools::load_all("/home/francisco/DATASHIELD/dsRctFed")

library(DSLite)
library(MASS)

# Parameters ------------------------------------------------------------------

set.seed(123)

K <- 5  # number of studies
d <- 10 # dimension
n_per_study <- 20 * d
p_treatment <- c(0.9, 0.9, 0.9, 0.1, 0.1)

sigma_noise <- 1.0
c_0 <- -2
c_1 <- -1.85

beta_1 <- seq(-1.75, 0.5, by = 0.25)
beta_0 <- seq(-1.8, 0, by = 0.2)

mu_X    <- c(rep(1, d %/% 2), rep(-1, d %/% 2))
I_d     <- diag(d)
J       <- matrix(1, nrow = d, ncol = d)
Sigma_X <- 0.5 * I_d + 0.5 * J

true_ATE <- (c_1 - c_0) + sum(mu_X * (beta_1 - beta_0))
print(true_ATE)

# Data generation -------------------------------------------------------------

data_list <- list()

for (k in 1:K) {
  n_k <- n_per_study
  X_k <- mvrnorm(n = n_k, mu = mu_X, Sigma = Sigma_X)
  W_k <- rbinom(n = n_k, size = 1, prob = p_treatment[k])

  epsilon_0 <- rnorm(n_k, mean = 0, sd = sigma_noise)
  epsilon_1 <- rnorm(n_k, mean = 0, sd = sigma_noise)

  Y_0 <- c_0 + X_k %*% beta_0 + epsilon_0
  Y_1 <- c_1 + X_k %*% beta_1 + epsilon_1
  Y_obs <- W_k * Y_1 + (1 - W_k) * Y_0

  data_list[[k]] <- data.frame(
    subject_id = 1:n_k,
    W = W_k,
    Y = as.vector(Y_obs),
    X_k
  )
  colnames(data_list[[k]]) <- c("subject_id", "W", "Y", paste0("X", 1:d))
}

pooled_data <- do.call(rbind, data_list)

# DSLite servers (one per study) ----------------------------------------------

for (k in 1:K) {
  server_name <- paste0("dslite.server", k)
  server_obj  <- newDSLiteServer(tables = list(study_data = data_list[[k]]))
  server_obj$config(defaultDSConfiguration(include = c("dsBase", "dsRctFed")))
  server_obj$aggregateMethod("localAteDS", "localAteDS")
  assign(server_name, server_obj, envir = .GlobalEnv)
}

builder <- DSI::newDSLoginBuilder()
for (k in 1:K) {
  builder$append(
    server = paste0("study", k),
    url    = paste0("dslite.server", k),
    table  = "study_data",
    driver = "DSLiteDriver"
  )
}
logindata <- builder$build()

conns <- DSI::datashield.login(logins = logindata, assign = FALSE)
datashield.assign.table(conns = conns, table = "study_data", symbol = "base")

# Federated ATE estimation ----------------------------------------------------

base_name       <- "base"
outcome_name    <- "Y"
treatment_name  <- "W"
covariate_names <- paste0("X", 1:d)

federated_ate <- ds.EstimateAte(
  datasources     = conns,
  base_name       = base_name,
  outcome_name    = outcome_name,
  treatment_name  = treatment_name,
  covariate_names = covariate_names
)

# Pooled (centralised) estimator for comparison --------------------------------

pooled_data_treated <- pooled_data[pooled_data$W == 1, ]
pooled_data_control <- pooled_data[pooled_data$W == 0, ]

formula <- paste0(outcome_name, " ~ ", paste(covariate_names, collapse = " + "))

model_treated <- glm(formula = formula, family = gaussian, data = pooled_data_treated)
model_control <- glm(formula = formula, family = gaussian, data = pooled_data_control)

coef_treated <- model_treated$coefficients
coef_control <- model_control$coefficients

alpha1 <- coef_treated["(Intercept)"]
alpha0 <- coef_control["(Intercept)"]
beta1  <- coef_treated[covariate_names]
beta0  <- coef_control[covariate_names]

tau_diff   <- alpha1 - alpha0
gamma_diff <- beta1 - beta0

X_mat    <- as.matrix(pooled_data[, covariate_names, drop = FALSE])
x_bar    <- colMeans(X_mat, na.rm = TRUE)
pooled_ate <- as.numeric(tau_diff + sum(x_bar * gamma_diff))

cat("Federated ATE:  ", federated_ate$ate, "\n")
cat("95% CI:         [", federated_ate$ci_lower, ",", federated_ate$ci_upper, "]\n")
cat("Variance:       ", federated_ate$variance, "\n")
print(abs(pooled_ate - federated_ate$ate))

DSI::datashield.logout(conns)
