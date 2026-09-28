#' Estimate ATE from multiple RCTs
#'
#' Client-side function that orchestrates the federated estimation of an
#' average treatment effect (ATE) — and, optionally, its variance and
#' confidence interval — from multiple independent RCTs, each held on a
#' separate data server. Randomisation holds within each study; the combined
#' dataset is not itself a single RCT.
#'
#' The function fits separate arm-specific linear models via
#' \code{\link[dsBaseClient]{ds.glm}}, broadcasts the pooled coefficients back
#' to the servers, and collects local statistics via a single aggregate call
#' to \code{localAteDS}.  When \code{compute_variance = FALSE}, only the three
#' scalars strictly needed for the point estimate are returned per server
#' (\eqn{n_k}, \eqn{n_k^{(1)}}, \eqn{\hat\tau_k}), omitting the residual
#' sum of squares and projected scatter that constitute additional disclosure.
#'
#' \strong{Server-side objects created.}
#' Two symbols are always created per server: \code{base_name_treated}
#' (treated-arm subset) and \code{base_name_control} (control-arm subset),
#' both used by \code{ds.glm}.  When \code{study_effects = TRUE} is set,
#' a third symbol is also created: \code{base_name_study_adj}, holding the
#' data frame augmented with the K-1 study dummy columns.  All three
#' persist in the server workspace after the function returns, as is typical
#' for DataSHIELD assign operations.  The caller's original
#' \code{base_name} object is never modified.
#'
#' When \code{compute_variance = TRUE} (the default), a delta-method variance
#' estimate and 95% confidence interval are also returned.  The formula used
#' depends on the \code{heterogeneous} argument; see below.
#'
#' @param datasources A list of \code{\link[DSI]{DSConnection-class}} objects.
#'   If \code{NULL}, the default connections are used.
#' @param base_name A character string naming the server-side data frame
#'   (must be identical across all servers).
#' @param outcome_name A character string naming the outcome variable.
#' @param treatment_name A character string naming the binary treatment
#'   indicator.
#' @param covariate_names A character vector of covariate names.
#' @param compute_variance Logical. If \code{TRUE} (default), the function
#'   returns a delta-method variance estimate and a 95% confidence interval in
#'   addition to the ATE point estimate.  To compute these, each server returns
#'   additional statistics (residual sum of squares and projected scatter
#'   scalar), and — when \code{heterogeneous = TRUE} — the per-arm Gram
#'   inverses are obtained by calling \code{ds.glm} with
#'   \code{viewVarCov = TRUE}, and each server additionally returns covariate
#'   column sums.  All of this constitutes an extra disclosure step.  Set to
#'   \code{FALSE} to return only the ATE point estimate and minimise the
#'   information shared by each server.
#' @param study_effects Logical. If \code{FALSE} (default), the no-study-effects
#'   model is assumed (potential outcome model defined in equation (2) of Khellaf et al., 2025):
#'   there is no arrow H -> Y in the DAG, meaning study membership does not
#'   directly shift the outcome.  If \code{TRUE}, the study-effects model is
#'   used (potential outcome model (16) of Khellaf et al., 2025): an arrow
#'   H -> Y is present, so site-specific intercept shifts confound the estimate
#'   when treatment probabilities also vary across sites.
#'   K-1 integer dummy columns (\code{study_id2}, \ldots, \code{study_idK})
#'   are appended to each server's data frame by \code{addStudyIdDS()}, and
#'   the federated GLMs are fitted on the augmented design matrix, absorbing
#'   the shifts.  The asymptotic variance is identical to the no-study-effects
#'   case (Appendix B.4.2 of Khellaf et al., 2025); the degrees of freedom
#'   adjust from n - 2(d+1) to n - 2(d+K) automatically.  Has no effect and
#'   emits a warning when only one datasource is present.
#' @param heterogeneous Logical. Controls which variance formula is used.
#'   \code{TRUE} (default) makes no assumption about the participating RCTs:
#'   each site may have a different covariate distribution and a different
#'   treatment probability.  The leverage-corrected OLS noise component is
#'   computed from the per-arm Gram inverses returned by
#'   \code{ds.glm(viewVarCov = TRUE)}, using the global covariate mean
#'   assembled from per-server column sums returned by \code{localAteDS()}.
#'   \code{FALSE} assumes that all RCTs share the same covariate distribution
#'   and uses a simpler closed-form variance that is exact only under this
#'   homogeneity assumption (constant covariate means and treatment
#'   probabilities across sites).
#'   Has no effect when \code{compute_variance = FALSE}.
#' @param heteroscedastic_arms Logical. If \code{FALSE} (default), the noise
#'   variance \eqn{\sigma^2} is assumed identical in both arms and estimated
#'   from the pooled residual sum of squares across treated and control
#'   observations.  If \code{TRUE}, arm-specific noise variances
#'   \eqn{\sigma_1^2} and \eqn{\sigma_0^2} are estimated separately from each
#'   arm's residual sum of squares, and the variance formula is adjusted
#'   accordingly.  Has no effect when \code{compute_variance = FALSE}.
#'
#' @return A list. When \code{compute_variance = TRUE} (the default),
#'   seven elements:
#'   \describe{
#'     \item{ate}{Federated ATE estimate \eqn{\hat\tau}.}
#'     \item{variance}{Delta-method variance estimate \eqn{\hat V(\hat\tau)}.}
#'     \item{ci_lower}{Lower 95-percent confidence interval bound.}
#'     \item{ci_upper}{Upper 95-percent confidence interval bound.}
#'     \item{n_total}{Total sample size across all servers.}
#'     \item{n_k}{Integer vector of per-server sample sizes.}
#'     \item{ate_k}{Numeric vector of per-server local ATE estimates.}
#'   }
#'   When \code{compute_variance = FALSE}, four elements:
#'   \code{ate}, \code{n_total}, \code{n_k}, \code{ate_k}.
#'   The structure is the same regardless of \code{study_effects}.
#' @references
#' Khellaf, R., Bellet, A., and Josse, J. (2025).
#' Federated Causal Inference: Multi-Study ATE Estimation beyond Meta-Analysis.
#' \emph{Proceedings of the 28th International Conference on Artificial
#' Intelligence and Statistics (AISTATS 2025)}.
#' @importFrom DSI datashield.connections_find
#' @importFrom dsBaseClient ds.glm ds.dataFrameSubset
#' @examples
#' \dontrun{
#' conns <- DSI::datashield.login(logins = logindata, assign = TRUE)
#' # Full result with variance and CI
#' result <- ds.EstimateAte(
#'   datasources = conns,
#'   base_name = "study_data",
#'   outcome_name = "Y",
#'   treatment_name = "W",
#'   covariate_names = c("X1", "X2")
#' )
#' # Point estimate only (less disclosure)
#' result_ate <- ds.EstimateAte(
#'   datasources = conns,
#'   base_name = "study_data",
#'   outcome_name = "Y",
#'   treatment_name = "W",
#'   covariate_names = c("X1", "X2"),
#'   compute_variance = FALSE
#' )
#' DSI::datashield.logout(conns)
#' }
#' @export
ds.EstimateAte <- function(datasources = NULL, base_name, outcome_name,
                           treatment_name, covariate_names,
                           compute_variance = TRUE,
                           study_effects = FALSE,
                           heterogeneous = TRUE,
                           heteroscedastic_arms = FALSE) {
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }

  .check_datasources(datasources)
  .check_args(base_name, outcome_name, treatment_name, covariate_names)

  # Study-effects setup ---------------------------------------------------------

  K <- length(datasources)

  if (study_effects && K == 1L) {
    cli::cli_warn("`study_effects = TRUE` has no effect with a single datasource.")
    study_effects <- FALSE
  }

  if (study_effects) {
    # Integer dummy columns (not a factor) are used so that ds.glm() does not
    # trigger R's factor-contrast machinery: model.frame() drops empty factor
    # levels to a 1-level factor on each server, causing "contrasts can be
    # applied only to factors with >= 2 levels".
    base_name_adj   <- paste0(base_name, "_study_adj")
    .add_study_id(datasources, base_name, K, base_name_adj)
    base_name       <- base_name_adj
    covariate_names <- c(covariate_names, paste0("study_id", seq(2L, K)))
  }

  # Subset into arm-specific server-side objects --------------------------------

  treatment_indicator <- paste0(base_name, "$", treatment_name)
  treated_name        <- paste0(base_name, "_treated")
  control_name        <- paste0(base_name, "_control")
  glm_formula <- paste0(outcome_name, " ~ ", paste(covariate_names, collapse = " + "))

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

  need_vcov <- compute_variance && heterogeneous

  model_treated <- suppressWarnings(dsBaseClient::ds.glm(
    formula     = glm_formula,
    data        = treated_name,
    family      = "gaussian",
    viewVarCov  = need_vcov,
    datasources = datasources
  ))
  model_control <- suppressWarnings(dsBaseClient::ds.glm(
    formula     = glm_formula,
    data        = control_name,
    family      = "gaussian",
    viewVarCov  = need_vcov,
    datasources = datasources
  ))

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
    datasources      = datasources,
    base_name        = base_name,
    outcome_name     = outcome_name,
    covariate_names  = covariate_names,
    treatment_name   = treatment_name,
    alpha0           = alpha0,
    beta0            = beta0,
    alpha1           = alpha1,
    beta1            = beta1,
    compute_variance = compute_variance,
    heterogeneous    = heterogeneous
  )

  n_k   <- local_results$n_k
  n_k_1 <- local_results$n_k_1
  tau_k <- local_results$tau_k

  names(n_k)   <- names(datasources)
  names(tau_k) <- names(datasources)

  # Aggregation -----------------------------------------------------------------

  n       <- sum(n_k)
  n_1     <- sum(n_k_1)
  ate_hat <- sum(n_k * tau_k) / n

  d <- length(covariate_names)  # with dummies when study_effects = TRUE
  .check_sample_sizes(n_k - n_k_1, n_k_1, d)

  if (!compute_variance) {
    return(list(
      ate     = ate_hat,
      n_total = n,
      n_k     = n_k,
      ate_k   = tau_k
    ))
  }

  # When study_effects = TRUE, d already includes the K-1 dummies, so the
  # denominator automatically becomes n - 2(d_orig + K) as required.
  # The same logic applies to arm-specific df: n_w - (d + 1).
  n_0 <- n - n_1
  if (heteroscedastic_arms) {
    sigma2_hat_1 <- sum(local_results$rss_k_1) / (n_1 - (d + 1L))
    sigma2_hat_0 <- sum(local_results$rss_k_0) / (n_0 - (d + 1L))
  } else {
    sigma2_hat <- sum(local_results$rss_k_1 + local_results$rss_k_0) /
                  (n - 2L * (d + 1L))
  }

  if (heterogeneous) {
    variance <- .compute_variance_delta(
      n            = n,
      n_k          = n_k,
      tau_k        = tau_k,
      ate_hat      = ate_hat,
      g_treated    = local_results$g1,
      g_control    = local_results$g0,
      vcov_treated = if (heteroscedastic_arms) sigma2_hat_1 * model_treated$VarCovMatrix
                     else                    sigma2_hat   * model_treated$VarCovMatrix,
      vcov_control = if (heteroscedastic_arms) sigma2_hat_0 * model_control$VarCovMatrix
                     else                    sigma2_hat   * model_control$VarCovMatrix,
      s_k          = local_results$s_k
    )
  } else {
    p_hat <- n_1 / n
    tss   <- sum(local_results$s_k + n_k * (tau_k - ate_hat)^2)
    if (heteroscedastic_arms) {
      variance <- sigma2_hat_1 / (n * p_hat) + sigma2_hat_0 / (n * (1 - p_hat)) +
                  tss / (n * (n - 1L))
    } else {
      variance <- sigma2_hat / (n * p_hat * (1 - p_hat)) + tss / (n * (n - 1L))
    }
  }

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
#' @importFrom DSI datashield.assign
.add_study_id <- function(datasources, base_name, K, newobj) {
  for (k in seq_along(datasources)) {
    cally <- call("addStudyIdDS", base_name, as.character(k), as.character(K))
    DSI::datashield.assign(datasources[k], newobj, cally)
  }
  invisible(NULL)
}

#' @keywords internal
#' @noRd
#' @importFrom DSI datashield.aggregate
.localAte <- function(datasources, base_name, outcome_name, covariate_names,
                      treatment_name, alpha0, beta0, alpha1, beta1,
                      compute_variance = TRUE, heterogeneous = TRUE) {
  cally <- call(
    "localAteDS",
    base_name, outcome_name,
    paste(covariate_names, collapse = ","),
    treatment_name,
    as.character(alpha0),
    paste(as.character(beta0), collapse = ","),
    as.character(alpha1),
    paste(as.character(beta1), collapse = ","),
    compute_variance,
    heterogeneous
  )
  results <- DSI::datashield.aggregate(conns = datasources, expr = cally)

  out <- list(
    n_k   = vapply(results, `[[`, 1L,  "n_k"),
    n_k_1 = vapply(results, `[[`, 1L,  "n_k_1"),
    tau_k = vapply(results, `[[`, 1.0, "tau_k")
  )
  if (compute_variance) {
    out$s_k     <- vapply(results, `[[`, 1.0, "s_k")
    out$rss_k_1 <- vapply(results, `[[`, 1.0, "rss_k_1")
    out$rss_k_0 <- vapply(results, `[[`, 1.0, "rss_k_0")
  }
  if (compute_variance && heterogeneous) {
    n      <- sum(out$n_k)
    # Both gradients equal (1, x_bar): d_tau/d_theta1 = +(1, x_bar) and
    # d_tau/d_theta0 = -(1, x_bar), but the sign cancels in the quadratic form
    # g^T V g, so we use the same positive vector for both arms.
    g_bar  <- Reduce("+", lapply(results, `[[`, "g1_raw_k")) / n
    out$g1 <- g_bar
    out$g0 <- g_bar
  }
  out
}

#' @keywords internal
#' @noRd
.compute_variance_delta <- function(n, n_k, tau_k, ate_hat,
                                    g_treated, g_control,
                                    vcov_treated, vcov_control,
                                    s_k) {
  param_term <- as.numeric(g_treated %*% vcov_treated %*% g_treated) +
                as.numeric(g_control %*% vcov_control %*% g_control)
  tss        <- sum(s_k + n_k * (tau_k - ate_hat)^2)
  samp_term  <- tss / (n * (n - 1L))
  param_term + samp_term
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
