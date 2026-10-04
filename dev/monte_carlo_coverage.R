# Monte Carlo coverage check for ds.EstimateAte variance estimator.
#
# Runs B simulations per scenario and checks that the nominal 95% CI covers
# the true ATE at approximately 95%.  With B = 500 the SE of estimated coverage
# is ~1%, so a 3-SE tolerance accepts coverage in roughly [92%, 98%].
#
# Scenarios
#   1. Continuous / homogeneous  / K = 1   (heterogeneous = FALSE)
#   2. Continuous / heterogeneous / K = 3  (heterogeneous = TRUE,  different site
#                                           means and treatment probabilities)
#   3. Binary     / K = 1                  (outcome_type = "binary")

library(DSLite)
library(DSI)
devtools::load_all("../dsFedAte")
devtools::load_all(".")

options(
  datashield.env                 = .GlobalEnv,
  datashield.privacyControlLevel = "banana",
  nfilter.tab                    = 3,
  nfilter.subset                 = 3,
  nfilter.glm                    = 0.33,
  nfilter.string                 = 80,
  nfilter.stringShort            = 20,
  nfilter.kNN                    = 3,
  nfilter.levels.density         = 0.33,
  nfilter.levels.max             = 40,
  nfilter.noise                  = 0.25,
  nfilter.privacy.old            = 5
)

# ── Shared helpers ─────────────────────────────────────────────────────────────

make_conns <- function(sites) {
  srv_names <- paste0("mc_srv_", seq_along(sites))
  builder   <- DSI::newDSLoginBuilder()
  for (k in seq_along(sites)) {
    srv <- DSLite::newDSLiteServer(tables = list(study_data = sites[[k]]))
    srv$config(defaultDSConfiguration(include = c("dsBase", "dsFedAte")))
    srv$aggregateMethod("localAteDS",                        "localAteDS")
    srv$aggregateMethod("listDisclosureSettingsDS",           "listDisclosureSettingsDS")
    srv$aggregateMethod("checkPermissivePrivacyControlLevel", "checkPermissivePrivacyControlLevel")
    srv$assignMethod("addStudyIdDS",                         "addStudyIdDS")
    assign(srv_names[k], srv, envir = .GlobalEnv)
    builder$append(server = srv_names[k], url = srv_names[k],
                   table = "study_data", driver = "DSLiteDriver")
  }
  conns <- DSI::datashield.login(logins = builder$build(), assign = FALSE)
  DSI::datashield.assign.table(conns = conns, symbol = "study_data",
                               table = "study_data")
  conns
}

run_coverage <- function(B, make_sites, true_ate, label, ...) {
  cat(sprintf("\n── %s (%d sims) ──\n", label, B))
  t0     <- proc.time()[["elapsed"]]
  covers <- logical(B)
  for (b in seq_len(B)) {
    if (b %% 100 == 0) cat(sprintf("  %d / %d\n", b, B))
    conns  <- make_conns(make_sites(seed = b))
    result <- tryCatch(
      suppressWarnings(ds.EstimateAte(
        datasources    = conns,
        base_name      = "study_data",
        outcome_name   = "Y",
        treatment_name = "W",
        ...
      )),
      error = function(e) NULL
    )
    DSI::datashield.logout(conns)
    if (!is.null(result) && !is.null(result$ci_lower)) {
      covers[b] <- result$ci_lower <= true_ate && true_ate <= result$ci_upper
    }
  }
  elapsed  <- proc.time()[["elapsed"]] - t0
  coverage <- mean(covers)
  se       <- sqrt(coverage * (1 - coverage) / B)
  pass     <- abs(coverage - 0.95) <= 3 * se
  cat(sprintf("  Coverage : %.1f%%  (SE ±%.1f%%, target 95%%)  →  %s  [%.0fs]\n",
              100 * coverage, 100 * se, if (pass) "PASS" else "FAIL", elapsed))
  invisible(list(coverage = coverage, se = se, pass = pass))
}

B <- 500

# ── Scenario 1: Continuous / homogeneous / K=1 ────────────────────────────────
#
# Y = 1 + 0.5*W + 0.3*X1 - 0.2*X2 + eps,  eps ~ N(0, 0.8^2)
# ATE = 0.5 exactly (intercept shift).

TRUE_ATE_CONT <- 0.5

r1 <- run_coverage(
  B          = B,
  make_sites = function(seed, n = 300) {
    set.seed(seed)
    df   <- data.frame(X1 = rnorm(n), X2 = rnorm(n), W = rbinom(n, 1, 0.5))
    df$Y <- 1.0 + TRUE_ATE_CONT * df$W + 0.3 * df$X1 - 0.2 * df$X2 +
            rnorm(n, sd = 0.8)
    list(df)
  },
  true_ate        = TRUE_ATE_CONT,
  label           = "Continuous / homogeneous / K=1",
  covariate_names = c("X1", "X2"),
  heterogeneous   = FALSE
)

# ── Scenario 2: Continuous / heterogeneous / K=3 ──────────────────────────────
#
# Same outcome model (beta1 == beta0 so ATE = alpha1 - alpha0 = 0.5 exactly),
# but with different site covariate means and treatment probabilities so that
# the leverage correction in the heterogeneous formula is needed.
#
# Site 1: X_bar = (1,  0),  p = 0.7
# Site 2: X_bar = (0,  0),  p = 0.5
# Site 3: X_bar = (-1, 1),  p = 0.3

r2 <- run_coverage(
  B          = B,
  make_sites = function(seed) {
    set.seed(seed)
    ns  <- c(120, 100, 130)
    mus <- list(c(1, 0), c(0, 0), c(-1, 1))
    ps  <- c(0.7, 0.5, 0.3)
    b   <- c(0.3, -0.2)
    lapply(seq_along(ns), function(k) {
      X1 <- rnorm(ns[k], mean = mus[[k]][1])
      X2 <- rnorm(ns[k], mean = mus[[k]][2])
      W  <- rbinom(ns[k], 1, ps[k])
      Y  <- 1.0 + TRUE_ATE_CONT * W + b[1] * X1 + b[2] * X2 +
            rnorm(ns[k], sd = 0.8)
      data.frame(Y = Y, W = W, X1 = X1, X2 = X2)
    })
  },
  true_ate        = TRUE_ATE_CONT,
  label           = "Continuous / heterogeneous / K=3",
  covariate_names = c("X1", "X2"),
  heterogeneous   = TRUE
)

# ── Scenario 3: Binary / K=1 ──────────────────────────────────────────────────
#
# logit P(Y=1|X,W=w) = alpha_w + X'beta,  X ~ N(0,I)
# True ATE = E_X[sigma(alpha1 + X'beta) - sigma(alpha0 + X'beta)]
# Computed by large-sample Monte Carlo (1e6 draws).

A0_BIN <- -0.5
A1_BIN <-  0.2
B_BIN  <- c(0.4, -0.3)

{
  set.seed(0)
  X_pop        <- cbind(rnorm(1e6), rnorm(1e6))
  TRUE_ATE_BIN <- mean(plogis(A1_BIN + X_pop %*% B_BIN) -
                       plogis(A0_BIN + X_pop %*% B_BIN))
  cat(sprintf("\nBinary true ATE (MC, 1e6 draws): %.6f\n", TRUE_ATE_BIN))
}

r3 <- run_coverage(
  B          = B,
  make_sites = function(seed, n = 400) {
    set.seed(seed)
    df   <- data.frame(X1 = rnorm(n), X2 = rnorm(n), W = rbinom(n, 1, 0.5))
    df$Y <- ifelse(
      df$W == 1,
      rbinom(n, 1, plogis(A1_BIN + B_BIN[1] * df$X1 + B_BIN[2] * df$X2)),
      rbinom(n, 1, plogis(A0_BIN + B_BIN[1] * df$X1 + B_BIN[2] * df$X2))
    )
    list(df)
  },
  true_ate        = TRUE_ATE_BIN,
  label           = "Binary / K=1",
  covariate_names = c("X1", "X2"),
  outcome_type    = "binary"
)

# ── Summary ────────────────────────────────────────────────────────────────────

cat("\n══ Monte Carlo coverage summary ══════════════════════════════\n")
results <- list(
  "Continuous / homogeneous  / K=1" = r1,
  "Continuous / heterogeneous / K=3" = r2,
  "Binary / K=1"                    = r3
)
for (nm in names(results)) {
  r <- results[[nm]]
  cat(sprintf("  %-38s  %.1f%%  %s\n",
              nm, 100 * r$coverage, if (r$pass) "PASS" else "FAIL"))
}
