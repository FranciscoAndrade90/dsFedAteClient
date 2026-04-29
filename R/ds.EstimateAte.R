#' Estimate average treatment effect in federated RCT data
#'
#' Client-side function that orchestrates the federated computation of an average
#' treatment effect (ATE) for data collected in a randomized controlled trial (RCT)
#' setting using DataSHIELD. The function uses \code{\link[dsBaseClient]{ds.glm}} to perform
#' separate linear regressions of the outcome on covariates for treatment and
#' control groups across distributed data sources. Local ATEs are computed at each
#' server using the regression coefficients, then aggregated proportionally to each
#' server's sample size to produce the final pooled estimate.
#' @param connections A list of \code{\link[DSI]{DSConnection-class}} objects obtained after login.
#' @param base_name A character string specifying the name of the data frame on
#'   the server side. Must be identical across all servers.
#' @param outcome_name A character string specifying the name of the outcome
#'   variable in the server-side data frame.
#' @param treatment_name A character string specifying the name of the treatment
#'   variable in the server-side data frame.
#' @param covariate_names A character vector specifying the names of covariate
#'   variables in the server-side data frame.
#'
#' @returns A numeric value representing the pooled average treatment effect across all federated data sources.
#' @export
ds.EstimateAte <- function(connections, base_name, outcome_name, treatment_name, covariate_names) {
  # Regression over the treated/control groups--------------------------------------

  treatment_indicator <- paste0(base_name, "$", treatment_name)
  treated_name <- paste0(base_name, "_treated")
  control_name <- paste0(base_name, "_control")
  formula <- paste0(outcome_name, " ~ ", paste(covariate_names, collapse = " + "))

  dsBaseClient::ds.dataFrameSubset(
    df.name          = base_name,
    V1.name          = treatment_indicator,
    V2.name          = "1",
    Boolean.operator = "==",
    newobj           = treated_name,
    datasources      = connections
  )

  dsBaseClient::ds.dataFrameSubset(
    df.name          = base_name,
    V1.name          = treatment_indicator,
    V2.name          = "0",
    Boolean.operator = "==",
    newobj           = control_name,
    datasources      = connections
  )

  model_treated <- suppressWarnings(dsBaseClient::ds.glm(
    formula     = formula,
    data        = treated_name,
    family      = "gaussian",
    datasources = connections
  ))

  model_control <- suppressWarnings(dsBaseClient::ds.glm(
    formula     = formula,
    data        = control_name,
    family      = "gaussian",
    datasources = connections
  ))
  ## Extract Coefficients ----------------------------------------------------------

  coef_treated <- model_treated$coefficients
  coef_control <- model_control$coefficients

  alpha1 <- coef_treated["(Intercept)", "Estimate"]
  alpha0 <- coef_control["(Intercept)", "Estimate"]
  beta1 <- coef_treated[covariate_names, "Estimate"]
  beta0 <- coef_control[covariate_names, "Estimate"]

  tau_diff <- alpha1 - alpha0
  gamma_diff <- beta1 - beta0

  # Send coefficients and ask for local ATE ----------------------------------------

  local_results <- ds.localAte(
    connections = connections,
    base_name = base_name,
    covariate_names = covariate_names,
    treatment_name = treatment_name,
    tau = tau_diff,
    gamma = gamma_diff
  )

  # Final aggregation --------------------------------------------------------------

  n_k <- local_results$n_k
  n <- sum(n_k)
  sw_weights <- n_k / n
  ate_hat <- sum(sw_weights * local_results$tau_k)
  return(ate_hat)
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

  return(list(
    n_k = n_k,
    tau_k = tau_k
  ))
}
