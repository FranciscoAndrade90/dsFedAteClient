#' Estimate average treatment effect in federated RCT data
#'
#' Client-side function that orchestrates the federated computation of an
#' average treatment effect (ATE) and, optionally, its variance
#' and confidence interval, for data collected in a randomised controlled trial
#' (RCT) setting using DataSHIELD.
#'
#' The function fits separate arm-specific linear models via
#' \code{\link[dsBaseClient]{ds.glm}}, broadcasts the pooled coefficients back
#' to the servers, and collects local statistics in a single additional
#' aggregate call (\code{localAteDS}).  Setting \code{compute_variance = FALSE}
#' limits the information returned by each server to the three scalars strictly
#' needed for the point estimate (\eqn{n_k}, \eqn{n_k^{(1)}},
#' \eqn{\hat\tau_k}), omitting the residual sum of squares and projected scatter
#' that are required for the variance but constitute additional disclosure.
#'
#' \strong{Server-side objects created.}  The function always creates two
#' symbols on every server: \code{<base_name>_treated} and
#' \code{<base_name>_control} (the arm-specific subsets used by
#' \code{ds.glm}).  When \code{study_effects = TRUE}, a third symbol
#' \code{<base_name>_study_adj} is also created, holding the original data
#' frame augmented with the K-1 integer study dummy columns.  All three
#' symbols persist in the server workspace after the function returns,
#' consistent with standard DataSHIELD practice.  The caller's original
#' \code{base_name} object is never modified.
#'
#' When \code{compute_variance = TRUE} (the default), a delta-method variance
#' estimate and 95\% confidence interval are also returned.  The formula used
#' depends on the \code{heterogeneous} argument; see below.
#'
#' @param datasources A list of \code{\link[DSI]{DSConnection-class}} objects.
#'   If \code{NULL}, uses \code{\link[DSI]{datashield.connections_find}}.
#' @param base_name A character string naming the server-side data frame
#'   (must be identical across all servers).
#' @param outcome_name A character string naming the outcome variable.
#' @param treatment_name A character string naming the binary treatment
#'   indicator.
#' @param covariate_names A character vector of covariate names.
#' @param compute_variance Logical. If \code{TRUE} (default), each server
#'   returns additional statistics beyond the point estimate (residual sum of
#'   squares and projected scatter scalar), and — when \code{heterogeneous = TRUE}
#'   — the per-arm Gram inverses are obtained by calling \code{ds.glm} with
#'   \code{viewVarCov = TRUE}, and each server additionally returns covariate
#'   column sums.  All of this constitutes an extra disclosure step.  Set to
#'   \code{FALSE} to request only the point estimate and minimise the
#'   information shared by each server.
#' @param study_effects Logical. If \code{FALSE} (default), the no-study-effects
#'   model is assumed (potential outcome model (2) of Khellaf et al., 2025):
#'   there is no arrow H -> Y in the DAG, meaning study membership does not
#'   directly shift the outcome.  If \code{TRUE}, the study-effects model is
#'   used (potential outcome model (16) of Khellaf et al., 2025): an arrow
#'   H -> Y is present, so site-specific intercept shifts confound the estimate
#'   when treatment probabilities also vary across sites.  \code{addStudyIdDS()} appends K-1 integer dummy columns
#'   (\code{study_id2}, \ldots, \code{study_idK}) to each server's data frame,
#'   and the federated GLMs are fitted on the augmented design matrix, absorbing
#'   the shifts.  The asymptotic variance is identical to the no-study-effects
#'   case (Appendix B.4.2 of Khellaf et al., 2025); the degrees of freedom
#'   adjust from n - 2(d+1) to n - 2(d+K) automatically.  Has no effect and
#'   emits a warning when only one datasource is present, or when
#'   \code{outcome_type = "binary"} (not currently supported for binary outcomes).
#' @param heterogeneous Logical. If \code{TRUE} (default), the leverage-corrected
#'   OLS noise component of the variance is computed from the per-arm Gram inverses returned by
#'   \code{ds.glm(viewVarCov = TRUE)}, using the global covariate mean assembled
#'   from per-server column sums returned by \code{localAteDS()}.  This is valid
#'   under full heterogeneity (different covariate distributions and treatment
#'   probabilities across sites) and is algebraically identical to the
#'   centralised plug-in variance.  Set to \code{FALSE} to use the simpler homogeneous closed form,
#'   which is exact only when treatment probabilities or covariate means are
#'   constant across sites.  Has no effect when \code{compute_variance = FALSE}.
#' @param outcome_type A character string, either \code{"continuous"} (default)
#'   or \code{"binary"}. Controls the family used in \code{ds.glm}
#'   (\code{"gaussian"} vs \code{"binomial"}) and the ATE formula in
#'   \code{localAteDS} (linear projection vs G-computation). When
#'   \code{"binary"} and \code{compute_variance = TRUE}, a delta-method
#'   variance is returned using the local Fisher information matrices and
#'   gradient vectors from each site; the \code{heterogeneous} argument has
#'   no effect for binary outcomes (the Hessian-based formula handles
#'   heterogeneity automatically).  Setting \code{outcome_type = "binary"}
#'   forces \code{study_effects = FALSE} with a warning.
#'
#' @return When \code{compute_variance = TRUE}, a named list with seven
#'   elements:
#'   \describe{
#'     \item{ate}{Federated ATE estimate \eqn{\hat\tau}.}
#'     \item{variance}{Delta-method variance estimate \eqn{\hat V(\hat\tau)}.}
#'     \item{ci_lower}{Lower bound of the approximate 95-percent confidence interval.}
#'     \item{ci_upper}{Upper bound of the approximate 95-percent confidence interval.}
#'     \item{n_total}{Total sample size across all servers.}
#'     \item{n_k}{Named integer vector of per-server sample sizes.}
#'     \item{ate_k}{Named numeric vector of per-server local ATE estimates.}
#'   }
#'   When \code{compute_variance = FALSE}, a named list with four elements:
#'   \code{ate}, \code{n_total}, \code{n_k}, \code{ate_k}.
#'   The structure is the same regardless of \code{study_effects}.
#' @references
#' Khellaf, R., Bellet, A., and Josse, J. (2025).
#' Federated Causal Inference: Multi-Study ATE Estimation beyond Meta-Analysis.
#' \emph{Proceedings of the 28th International Conference on Artificial
#' Intelligence and Statistics (AISTATS 2025)}.
#' @importFrom DSI datashield.connections_find
#' @examples
#' \dontrun{
#' conns <- DSI::datashield.login(logins = logindata, assign = TRUE)
#' # Full result with variance and CI
#' result <- ds.EstimateAte(
#'   datasources     = conns,
#'   base_name       = "study_data",
#'   outcome_name    = "Y",
#'   treatment_name  = "W",
#'   covariate_names = c("X1", "X2")
#' )
#' # Point estimate only (less disclosure)
#' result_ate <- ds.EstimateAte(
#'   datasources      = conns,
#'   base_name        = "study_data",
#'   outcome_name     = "Y",
#'   treatment_name   = "W",
#'   covariate_names  = c("X1", "X2"),
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
                           outcome_type = "continuous") {
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }

  .check_datasources(datasources)
  .check_args(base_name, outcome_name, treatment_name, covariate_names)

  outcome_type <- match.arg(outcome_type, c("continuous", "binary"))
  family <- if (outcome_type == "continuous") "gaussian" else "binomial"

  # Study-effects setup ---------------------------------------------------------

  K <- length(datasources)

  if (study_effects && outcome_type == "binary") {
    cli::cli_warn("`study_effects = TRUE` is not supported for binary outcomes and will be ignored.")
    study_effects <- FALSE
  }

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
  glm_formula         <- paste0(outcome_name, " ~ ", paste(covariate_names, collapse = " + "))

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

  need_vcov <- compute_variance && (outcome_type == "binary" || heterogeneous)

  model_treated <- .fit_arm_glm(glm_formula, treated_name, family, need_vcov, datasources)
  model_control <- .fit_arm_glm(glm_formula, control_name, family, need_vcov, datasources)

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
    heterogeneous    = heterogeneous,
    outcome_type     = outcome_type
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

  sigma2_hat <- if (outcome_type == "binary") 1.0
               else sum(local_results$rss_k) / (n - 2L * (d + 1L))

  if (heterogeneous || outcome_type == "binary") {
    variance <- .compute_variance_delta(
      n            = n,
      n_k          = n_k,
      tau_k        = tau_k,
      ate_hat      = ate_hat,
      g_treated    = local_results$g1,
      g_control    = local_results$g0,
      vcov_treated = sigma2_hat * model_treated$VarCovMatrix,
      vcov_control = sigma2_hat * model_control$VarCovMatrix,
      s_k          = local_results$s_k
    )
  } else {
    p_hat    <- n_1 / n
    tss      <- sum(local_results$s_k + n_k * (tau_k - ate_hat)^2)
    variance <- sigma2_hat / (n * p_hat * (1 - p_hat)) + tss / (n * (n - 1L))
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
#' @importFrom dsBaseClient ds.glm
.fit_arm_glm <- function(glm_formula, data_name, family, need_vcov, datasources) {
  suppressWarnings(dsBaseClient::ds.glm(
    formula     = glm_formula,
    data        = data_name,
    family      = family,
    viewVarCov  = need_vcov,
    datasources = datasources
  ))
}

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
                      compute_variance = TRUE, heterogeneous = TRUE,
                      outcome_type = "continuous") {
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
    heterogeneous,
    outcome_type
  )
  results <- DSI::datashield.aggregate(conns = datasources, expr = cally)

  out <- list(
    n_k   = vapply(results, `[[`, 1L,  "n_k"),
    n_k_1 = vapply(results, `[[`, 1L,  "n_k_1"),
    tau_k = vapply(results, `[[`, 1.0, "tau_k")
  )
  if (compute_variance) {
    out$s_k <- vapply(results, `[[`, 1.0, "s_k")
  }
  if (compute_variance && outcome_type == "continuous") {
    out$rss_k <- vapply(results, `[[`, 1.0, "rss_k")
  }
  if (compute_variance && (heterogeneous || outcome_type == "binary")) {
    n      <- sum(out$n_k)
    out$g1 <-  Reduce("+", lapply(results, `[[`, "g1_raw_k")) / n
    out$g0 <- -Reduce("+", lapply(results, `[[`, "g0_raw_k")) / n
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
