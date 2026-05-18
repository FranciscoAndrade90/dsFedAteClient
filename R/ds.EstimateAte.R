#' Estimate average treatment effect in federated RCT data
#'
#' Client-side function that orchestrates the federated computation of an average
#' treatment effect (ATE) for data collected in a randomized controlled trial (RCT)
#' setting using DataSHIELD. The function uses \code{\link[dsBaseClient]{ds.glm}} to perform
#' separate linear regressions of the outcome on covariates for treatment and
#' control groups across distributed data sources. Local ATEs are computed at each
#' server using the regression coefficients, then aggregated proportionally to each
#' server's sample size to produce the final pooled estimate.
#'
#' @param datasources A list of \code{\link[DSI]{DSConnection-class}} objects obtained after
#'   login. If \code{NULL}, uses \code{\link[DSI]{datashield.connections_find}}.
#' @param base_name A character string specifying the name of the data frame on
#'   the server side. Must be identical across all servers.
#' @param outcome_name A character string specifying the name of the outcome
#'   variable in the server-side data frame.
#' @param treatment_name A character string specifying the name of the treatment
#'   variable in the server-side data frame.
#' @param covariate_names A character vector specifying the names of covariate
#'   variables in the server-side data frame.
#'
#' @return A numeric value representing the pooled average treatment effect
#'   across all federated data sources.
#' @importFrom DSI datashield.connections_find
#' @examples
#' \dontrun{
#' conns <- DSI::datashield.login(logins = logindata, assign = TRUE)
#' ate <- ds.EstimateAte(
#'   datasources     = conns,
#'   base_name       = "study_data",
#'   outcome_name    = "Y",
#'   treatment_name  = "W",
#'   covariate_names = c("X1", "X2")
#' )
#' DSI::datashield.logout(conns)
#' }
#' @export
ds.EstimateAte <- function(datasources = NULL, base_name, outcome_name, treatment_name, covariate_names) {
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }

  .check_datasources(datasources)
  .check_args(base_name, outcome_name, treatment_name, covariate_names)

  # Regression over the treated/control groups ----------------------------------

  treatment_indicator <- paste0(base_name, "$", treatment_name)
  treated_name        <- paste0(base_name, "_treated")
  control_name        <- paste0(base_name, "_control")
  formula             <- paste0(outcome_name, " ~ ", paste(covariate_names, collapse = " + "))

  dsBaseClient::ds.dataFrameSubset(
    df.name          = base_name,
    V1.name          = treatment_indicator,
    V2.name          = "1",
    Boolean.operator = "==",
    newobj           = treated_name,
    datasources      = datasources
  )

  dsBaseClient::ds.dataFrameSubset(
    df.name          = base_name,
    V1.name          = treatment_indicator,
    V2.name          = "0",
    Boolean.operator = "==",
    newobj           = control_name,
    datasources      = datasources
  )

  model_treated <- dsBaseClient::ds.glm(
    formula     = formula,
    data        = treated_name,
    family      = "gaussian",
    datasources = datasources
  )

  model_control <- dsBaseClient::ds.glm(
    formula     = formula,
    data        = control_name,
    family      = "gaussian",
    datasources = datasources
  )

  ## Extract coefficients -------------------------------------------------------

  coef_treated <- model_treated$coefficients
  coef_control <- model_control$coefficients

  alpha1 <- coef_treated["(Intercept)", "Estimate"]
  alpha0 <- coef_control["(Intercept)", "Estimate"]
  beta1  <- coef_treated[covariate_names, "Estimate"]
  beta0  <- coef_control[covariate_names, "Estimate"]

  tau_diff   <- alpha1 - alpha0
  gamma_diff <- beta1 - beta0

  # Send coefficients and ask for local ATE ------------------------------------

  local_results <- .localAte(
    datasources     = datasources,
    base_name       = base_name,
    covariate_names = covariate_names,
    treatment_name  = treatment_name,
    tau             = tau_diff,
    gamma           = gamma_diff
  )

  # Final aggregation ----------------------------------------------------------

  n_k        <- local_results$n_k
  n          <- sum(n_k)
  sw_weights <- n_k / n
  ate_hat    <- sum(sw_weights * local_results$tau_k)

  return(ate_hat)
}


#' @keywords internal
#' @noRd
#' @importFrom DSI datashield.aggregate
.localAte <- function(datasources, base_name, covariate_names, treatment_name, tau, gamma) {
  K     <- length(datasources)
  n_k   <- integer(K)
  tau_k <- numeric(K)

  covariate_names_str <- paste(covariate_names, collapse = ",")
  gamma_str           <- paste(as.character(gamma), collapse = ",")
  tau_str             <- as.character(tau)

  cally   <- call("localAteDS", base_name, covariate_names_str, treatment_name, tau_str, gamma_str)
  results <- DSI::datashield.aggregate(conns = datasources, expr = cally)

  for (k in seq_len(K)) {
    res      <- results[[k]]
    n_k[k]   <- res$n_k
    tau_k[k] <- res$tau_k
  }

  return(list(
    n_k   = n_k,
    tau_k = tau_k
  ))
}


#' @keywords internal
#' @noRd
#' @importFrom methods is
#' @importFrom cli cli_abort
.check_datasources <- function(datasources) {
  if (!(is.list(datasources) &&
        all(unlist(lapply(datasources, function(d) methods::is(d, "DSConnection")))))) {
    cli::cli_abort("'datasources' must be a list of DSConnection objects.")
  }
}


#' @keywords internal
#' @noRd
#' @importFrom cli cli_abort
.check_args <- function(base_name, outcome_name, treatment_name, covariate_names) {
  .check_string(base_name, "base_name")
  .check_string(outcome_name, "outcome_name")
  .check_string(treatment_name, "treatment_name")
  if (is.null(covariate_names)) {
    cli::cli_abort("`covariate_names` must not be NULL")
  }
  if (!is.character(covariate_names) || length(covariate_names) == 0L) {
    cli::cli_abort(
      c(
        "x" = "`covariate_names` must be a non-empty character vector",
        "i" = "You provided an object of class {class(covariate_names)}"
      )
    )
  }
}


#' @keywords internal
#' @noRd
#' @importFrom cli cli_abort
.check_string <- function(x, arg_name) {
  if (is.null(x)) {
    cli::cli_abort("`{arg_name}` must not be NULL")
  }
  if (!is.character(x) || length(x) != 1L) {
    cli::cli_abort(
      c(
        "x" = "`{arg_name}` must be a single character string",
        "i" = "You provided an object of class {class(x)}"
      )
    )
  }
}
