#' Estimate average treatment effect in federated RCT data
#'
#' Client-side function that orchestrates the federated computation of an
#' average treatment effect (ATE) and its variance for data collected in a
#' randomised controlled trial (RCT) setting using DataSHIELD.
#'
#' The function fits separate arm-specific linear models via
#' \code{\link[dsBaseClient]{ds.glm}}, broadcasts the pooled coefficients back
#' to the servers, and collects the quantities needed for the delta-method
#' variance in a single additional aggregate call (\code{localAteDS}).
#'
#' The variance estimator has two components (see \code{federated_variance.pdf}
#' for derivation):
#' \enumerate{
#'   \item OLS noise: \eqn{\hat\sigma^2 / (n\,\hat p(1-\hat p))}, where
#'     \eqn{\hat\sigma^2} uses the pooled df \eqn{n - 2(d+1)}.
#'   \item Covariate-mean uncertainty:
#'     \eqn{\Delta\hat\beta^\top\hat\Sigma\Delta\hat\beta / n}, pooled from
#'     per-server scatter scalars.
#' }
#'
#' @param datasources A list of \code{\link[DSI]{DSConnection-class}} objects.
#'   If \code{NULL}, uses \code{\link[DSI]{datashield.connections_find}}.
#' @param base_name A character string naming the server-side data frame
#'   (must be identical across all servers).
#' @param outcome_name A character string naming the outcome variable.
#' @param treatment_name A character string naming the binary treatment
#'   indicator.
#' @param covariate_names A character vector of covariate names.
#'
#' @return A named list with seven elements:
#'   \describe{
#'     \item{ate}{Pooled GD ATE estimate \eqn{\hat\tau}.}
#'     \item{variance}{Delta-method variance estimate \eqn{\hat V(\hat\tau)}.}
#'     \item{ci_lower}{Lower bound of the approximate 95-percent confidence interval.}
#'     \item{ci_upper}{Upper bound of the approximate 95-percent confidence interval.}
#'     \item{n_total}{Total sample size across all servers.}
#'     \item{n_k}{Named integer vector of per-server sample sizes.}
#'     \item{ate_k}{Named numeric vector of per-server local ATE estimates.}
#'   }
#' @importFrom DSI datashield.connections_find
#' @examples
#' \dontrun{
#' conns <- DSI::datashield.login(logins = logindata, assign = TRUE)
#' result <- ds.EstimateAte(
#'   datasources     = conns,
#'   base_name       = "study_data",
#'   outcome_name    = "Y",
#'   treatment_name  = "W",
#'   covariate_names = c("X1", "X2")
#' )
#' DSI::datashield.logout(conns)
#' }
#' @export
ds.EstimateAte <- function(datasources = NULL, base_name, outcome_name,
                           treatment_name, covariate_names) {
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }

  .check_datasources(datasources)
  .check_args(base_name, outcome_name, treatment_name, covariate_names)

  # Subset into arm-specific server-side objects --------------------------------

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

  # Federated arm-specific OLS --------------------------------------------------

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

  .check_glm_validity(model_treated, model_control, covariate_names)

  # Extract individual arm coefficients -----------------------------------------

  coef_treated <- model_treated$coefficients
  coef_control <- model_control$coefficients

  alpha1 <- coef_treated["(Intercept)", "Estimate"]
  alpha0 <- coef_control["(Intercept)", "Estimate"]
  beta1  <- coef_treated[covariate_names, "Estimate"]
  beta0  <- coef_control[covariate_names, "Estimate"]

  # Collect local statistics from all servers -----------------------------------

  local_results <- .localAte(
    datasources     = datasources,
    base_name       = base_name,
    outcome_name    = outcome_name,
    covariate_names = covariate_names,
    treatment_name  = treatment_name,
    alpha0          = alpha0,
    beta0           = beta0,
    alpha1          = alpha1,
    beta1           = beta1
  )

  n_k        <- local_results$n_k
  n_k_1      <- local_results$n_k_1
  tau_k      <- local_results$tau_k
  s_k_scalar <- local_results$s_k_scalar
  rss_k      <- local_results$rss_k
  x_bar_k    <- local_results$x_bar_k   # K × d matrix

  names(n_k)   <- names(datasources)
  names(tau_k) <- names(datasources)

  # Aggregation -----------------------------------------------------------------

  n      <- sum(n_k)
  n_1    <- sum(n_k_1)
  gamma  <- beta1 - beta0

  sw_weights <- n_k / n
  ate_hat    <- sum(sw_weights * tau_k)

  .check_sample_sizes(n_k - n_k_1, n_k_1, length(covariate_names))

  variance <- .compute_variance(
    n          = n,
    n_1        = n_1,
    rss_k      = rss_k,
    s_k_scalar = s_k_scalar,
    n_k        = n_k,
    x_bar_k    = x_bar_k,
    gamma       = gamma
  )

  se       <- sqrt(variance)
  ci_lower <- ate_hat - 1.96 * se
  ci_upper <- ate_hat + 1.96 * se

  return(list(
    ate      = ate_hat,
    variance = variance,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    n_total  = n,
    n_k      = n_k,
    ate_k    = tau_k
  ))
}


# Internal helpers -------------------------------------------------------------

#' @keywords internal
#' @noRd
#' @importFrom DSI datashield.aggregate
.localAte <- function(datasources, base_name, outcome_name, covariate_names,
                      treatment_name, alpha0, beta0, alpha1, beta1) {
  K <- length(datasources)
  d <- length(beta0)

  n_k        <- integer(K)
  n_k_1      <- integer(K)
  tau_k      <- numeric(K)
  s_k_scalar <- numeric(K)
  rss_k      <- numeric(K)
  x_bar_k    <- matrix(0.0, nrow = K, ncol = d)

  covariate_names_str <- paste(covariate_names, collapse = ",")
  alpha0_str          <- as.character(alpha0)
  alpha1_str          <- as.character(alpha1)
  beta0_str           <- paste(as.character(beta0), collapse = ",")
  beta1_str           <- paste(as.character(beta1), collapse = ",")

  cally <- call(
    "localAteDS",
    base_name, outcome_name, covariate_names_str, treatment_name,
    alpha0_str, beta0_str, alpha1_str, beta1_str
  )
  results <- DSI::datashield.aggregate(conns = datasources, expr = cally)

  for (k in seq_len(K)) {
    res            <- results[[k]]
    n_k[k]        <- res$n_k
    n_k_1[k]      <- res$n_k_1
    tau_k[k]      <- res$tau_k
    s_k_scalar[k] <- res$s_k_scalar
    rss_k[k]      <- res$rss_k
    x_bar_k[k, ]  <- res$x_bar_k
  }

  list(
    n_k        = n_k,
    n_k_1      = n_k_1,
    tau_k      = tau_k,
    s_k_scalar = s_k_scalar,
    rss_k      = rss_k,
    x_bar_k    = x_bar_k
  )
}


#' @keywords internal
#' @noRd
#' @importFrom cli cli_warn
.check_sample_sizes <- function(n_k_0, n_k_1, d) {
  threshold <- 5L * (d + 1L)
  n_0 <- sum(n_k_0)
  n_1 <- sum(n_k_1)
  if (n_0 < threshold || n_1 < threshold) {
    cli::cli_warn(c(
      "!" = "Asymptotic approximation may be unreliable.",
      "i" = paste0(
        "Federated arm sizes n^(0)=", n_0, ", n^(1)=", n_1,
        " should both exceed 5*(d+1)=", threshold, "."
      )
    ))
  }
}


#' @keywords internal
#' @noRd
.compute_variance <- function(n, n_1, rss_k, s_k_scalar, n_k, x_bar_k, gamma) {
  d <- length(gamma)

  # First term: OLS noise with pooled df = n - 2(d+1)
  p_hat      <- n_1 / n
  sigma2_hat <- sum(rss_k) / (n - 2L * (d + 1L))
  first_term <- sigma2_hat / (n * p_hat * (1 - p_hat))

  # Second term: covariate-mean uncertainty Δβᵀ Σ̂ Δβ / n
  # proj_local_k = X̄_k ᵀ Δβ̂ (K-vector)
  proj_local  <- as.numeric(x_bar_k %*% gamma)
  proj_global <- sum(n_k * proj_local) / n
  tss         <- sum(s_k_scalar + n_k * (proj_local - proj_global)^2)
  second_term <- tss / (n * (n - 1L))

  first_term + second_term
}


#' @keywords internal
#' @noRd
#' @importFrom cli cli_abort
.check_args <- function(base_name, outcome_name, treatment_name, covariate_names) {
  .check_string(base_name,    "base_name")
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
.check_glm_validity <- function(model_treated, model_control, covariate_names) {
  .failed <- function(model) {
    errors <- model$errorMessage
    if (is.null(errors)) return(character(0))
    rownames(errors)[errors != "No errors"]
  }
  failed <- union(.failed(model_treated), .failed(model_control))
  if (length(failed) == 0L) return(invisible(NULL))
  n_params <- length(covariate_names) + 1L
  cli::cli_abort(c(
    "x" = "DataSHIELD blocked the GLM in {length(failed)} study/studies: {paste(failed, collapse = ', ')}.",
    "i" = "At least one treatment subset in those studies has too few observations relative to the {n_params} model parameters.",
    "i" = "DataSHIELD's nfilter.glm threshold requires each arm to have at least n_params / nfilter.glm observations (default: {n_params} / 0.33 = {ceiling(n_params / 0.33)}).",
    ">" = "Consider reducing the number of covariates, or exclude the affected studies via the `datasources` argument."
  ))
}

