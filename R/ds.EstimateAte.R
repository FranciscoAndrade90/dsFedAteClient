ds.EstimateAte <- function(connections, base_name, outcome_name, treatment_name, covariate_names) {
  # Regression ---------------------------------------------------------------------

  formula <- paste0(outcome_name, " ~ ", paste(covariate_names, collapse = " + "))

  ## Regression over the treated group ---------------------------------------------

  treatment_indicator <- paste0(base_name, "$", treatment_name)

  model_treated <- dsBaseClient::ds.glm(
    formula     = formula,
    data        = base_name,
    family      = "gaussian",
    weights     = treatment_indicator,
    datasources = connections
  )

  ## Regression ver the control group ---------------------------------------------

  control_indicator <- paste0("1 - ", base_name, "$", treatment_name)

  model_control <- dsBaseClient::ds.glm(
    formula     = formula,
    data        = base_name,
    family      = "gaussian",
    weights     = control_indicator,
    datasources = connections
  )

  ## Extract Coefficients ----------------------------------------------------------

  coef_treated <- model_treated$coefficients
  coef_control <- model_control$coefficients

  alpha1 <- coef_treated["(Intercept)", "Estimate"]
  alpha0 <- coef_control["(Intercept)", "Estimate"]
  beta1 <- coef_treated[covariate_names, "Estimate"]
  beta0 <- coef_control[covariate_names, "Estimate"]

  tau <- alpha1 - alpha0
  gamma <- beta1 - beta0

  # Send coefficients and ask for local ATE ----------------------------------------

  local_results <- ds.localAte(
    connections = connections,
    base_name = base_name,
    covariate_names = covariate_names,
    treatment_name = treatment_name,
    tau = tau,
    gamma = gamma
  )

  # Final aggregation --------------------------------------------------------------

  n <- sum(local_results$n_k)
  sw_weights <- n_k / n
  tau <- sum(sw_weights * local_results$tau_k)
}


ds.localAte <- function(connections, base_name, covariate_names, treatment_name, tau, gamma) {
  K <- length(connections)
  n_k <- integer(K)
  tau_k <- numeric(K)

  covariate_names_str <- paste(covariate_names, collapse = ",")
  gamma_str <- paste(as.character(gamma), collapse = ",")
  tau_str <- as.character(tau)

  cally <- call("localAteDS", base_name, covariate_names_str, treatment_name, tau_str, gamma_str)

  results <- DSI::datashield.aggregate(conns = connections, expr = cally)

  for (k in seq_len(K)) {
    res <- results[[k]]
    n_k[k] <- res$n_k
    tau_k[k] <- res$tau_k
  }

  list(
    n_k = n_k,
    tau_k = tau_k
  )
}
