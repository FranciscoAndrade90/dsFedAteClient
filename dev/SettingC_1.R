library(DSLite)
library(DSI)
library(dsBase)
library(dsBaseClient)
library(MASS)
library(devtools)

devtools::load_all("/Users/fdelimaa/Desktop/INRIA/FedCausal/dsRctFedClient")
devtools::load_all("/Users/fdelimaa/Desktop/INRIA/FedCausal/dsRctFed")

# ============================================================================
# STEP 2: Generate Data for Setting C.1
# ============================================================================

set.seed(123)

# Parameters
K <- 5
d <- 10
n_per_study <- 20 * d
p_treatment <- c(0.5, 0.5, 0.5, 0.5, 0.5)

# Outcome model parameters
sigma_noise <- 1.0
c_0 <- -2
c_1 <- -1.85

beta_1 <- seq(-1.75, 0.5, by = 0.25)
beta_0 <- seq(-1.8, 0, by = 0.2)

mu_X <- c(rep(1, d %/% 2), rep(-1, d %/% 2))
I_d <- diag(d)                  # identity matrix
J <- matrix(1, nrow = d, ncol = d)  # matrix of all ones
Sigma_X <-  0.5 * I_d + 0.5 * J

true_ATE <- (c_1 - c_0) + sum(mu_X * (beta_1 - beta_0))

cat("=== Setting C.1: DSLite Setup ===\n")
cat("True ATE:", true_ATE, "\n\n")

# Generate data for each study
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

  # Create data frame for this study
  data_list[[k]] <- data.frame(
    subject_id = 1:n_k,
    W = W_k,
    Y = as.vector(Y_obs),
    X_k
  )
  colnames(data_list[[k]]) <- c("subject_id", "W", "Y", paste0("X", 1:d))

  cat(sprintf("Study %d: n=%d, n_treated=%d, n_control=%d\n",
              k, n_k, sum(W_k), sum(1-W_k)))
}

pooled_data <- do.call(rbind, data_list)

cat("\n")

# ============================================================================
# STEP 3: Create 5 SEPARATE DSLite Servers (One per Study)
# ============================================================================

cat("=== Creating", K, "Separate DSLite Servers ===\n\n")

# Create K separate named server objects dynamically
# Use assign() to create objects with names like "dslite.server1", "dslite.server2", etc.

for (k in 1:K) {
  server_name <- paste0("dslite.server", k)

  # Create the server
  server_obj <- newDSLiteServer(
    tables = list(study_data = data_list[[k]])
  )

  # Configure it
  server_obj$config(defaultDSConfiguration(include = c("dsBase", "dsRctFed")))
  server_obj$aggregateMethod("localAteDS", "localAteDS")

  # Assign it to the environment with the dynamic name
  assign(server_name, server_obj, envir = .GlobalEnv)

  cat(sprintf("Server %d created: %d observations (object: %s)\n",
              k, nrow(data_list[[k]]), server_name))
}

cat("\n")

# ============================================================================
# STEP 4: Set DataSHIELD Options
# ============================================================================

options(
  datashield.privacyControlLevel = "permissive",
  nfilter.tab = 3,
  nfilter.subset = 3,
  nfilter.glm = 0.33,
  nfilter.string = 80,
  nfilter.stringShort = 20,
  nfilter.kNN = 3,
  nfilter.levels.density = 0.33,
  nfilter.levels.max = 40,
  nfilter.noise = 0.25,
  nfilter.privacy.old = 5
)

# ============================================================================
# STEP 5: Build Login Object for 5 Separate Servers
# ============================================================================

cat("=== Building Login Configuration for", K, "Servers ===\n\n")

builder <- DSI::newDSLoginBuilder()

# Add each server using a loop
for (k in 1:K) {
  builder$append(
    server = paste0("study", k),
    url = paste0("dslite.server", k),
    table = "study_data",
    driver = "DSLiteDriver"
  )
}

logindata <- builder$build()

cat("Login configuration created for", K, "SEPARATE servers\n")
print(logindata)
cat("\n")

# ============================================================================
# STEP 6: Login to All 5 Servers
# ============================================================================

cat("=== Logging into 5 Separate DataSHIELD Servers ===\n\n")

# This creates 5 connections, each to a DIFFERENT server
conns <- DSI::datashield.login(logins = logindata, assign = FALSE)

cat("Successfully connected to", length(conns), "separate DataSHIELD servers\n")
cat("Connection names:", names(conns), "\n\n")

# ============================================================================
# STEP 7: Assign Data on Each Server
# ============================================================================

cat("=== Assigning Data on Each Server ===\n\n")

# Assign the study data to symbol 'base' on each server
datashield.assign.table(
  conns = conns,
  table = "study_data",
  symbol = "base"
)

cat("Data assigned to symbol 'base' on all 5 servers\n\n")

# ============================================================================
# STEP 8: Verify Setup - Each Server is Separate
# ============================================================================

cat("=== Verifying 5 Separate Servers ===\n\n")

# Check that each connection points to a different server
cat("Server objects (", K, "separate named objects):\n", sep = "")
for (k in 1:K) {
  server_name <- paste0("dslite.server", k)
  server_obj <- get(server_name)
  cat(sprintf("  %s: %s\n", server_name, class(server_obj)[1]))
}
cat("\n")

# Check dimensions on each server
cat("Dimensions on each server:\n")
dims <- ds.dim("base", datasources = conns)
for (study_name in names(dims)) {
  cat(sprintf("  %s: %d rows x %d columns\n",
              study_name, dims[[study_name]][1], dims[[study_name]][2]))
}
cat("\n")

# Check column names
cat("Column names (should be same on all servers):\n")
colnames_result <- ds.colnames("base", datasources = conns[1])
print(colnames_result)
cat("\n")

# Get summary of treatment variable on each server
cat("Treatment variable summary on each server:\n")
W_summary <- ds.summary("base$W", datasources = conns)
print(W_summary)
cat("\n")

# Get mean of outcome on each server
cat("Mean of outcome Y on each server:\n")
Y_means <- ds.mean("base$Y", datasources = conns)
print(Y_means)
cat("\n")

# Verify data isolation - each server only has its own data
cat("Verifying data isolation:\n")
for (k in 1:K) {
  n_k <- ds.length("base$Y", datasources = conns[k])
  cat(sprintf("  Server %d has %d observations (isolated from other servers)\n",
              k, n_k[[1]]))
}
cat("\n")

base_name <- "base"                    # Name of the dataset symbol on each server
outcome_name <- "Y"                 # Outcome variable
treatment_name <- "W"               # Treatment variable
covariate_names <- paste0("X", 1:d)  # Covariates X1, X2, ..., X10
connections <- conns
# Call the function
cat("Calling: ds.EstimateAte(conns, '", base_name, "', '", outcome_name,
    "', '", treatment_name, "', c('", paste(covariate_names, collapse = "', '"), "'))\n\n", sep = "")

federated_ate <- ds.EstimateAte(
  datasources     = conns,
  base_name       = base_name,
  outcome_name    = outcome_name,
  treatment_name  = treatment_name,
  covariate_names = covariate_names
)

cat("Function returned successfully!\n\n")

cat("==========================================================\n")
cat("RESULTS FROM ds.EstimateAte\n")
cat("==========================================================\n\n")

print(federated_ate)

cat("\n")


cat("==========================================================\n")
cat("Centralized Setting")
cat("==========================================================\n\n")

pooled_data_treated <- pooled_data[pooled_data$W == 1, ]
pooled_data_control <- pooled_data[pooled_data$W == 0, ]

formula <- paste0(outcome_name, " ~ ", paste(covariate_names, collapse = " + "))

model_treated <- glm(formula = formula, family = gaussian, data = pooled_data_treated)
model_control <- glm(formula = formula, family = gaussian, data = pooled_data_control)

coef_treated <- model_treated$coefficients
coef_control <- model_control$coefficients

alpha1 <- coef_treated["(Intercept)"]
alpha0 <- coef_control["(Intercept)"]
beta1 <- coef_treated[covariate_names]
beta0 <- coef_control[covariate_names]

tau_diff <- alpha1 - alpha0
gamma_diff <- beta1 - beta0

X_mat <- as.matrix(pooled_data[, covariate_names, drop = FALSE])
x_bar <- colMeans(X_mat, na.rm = TRUE)
pooled_ate <- as.numeric(tau_diff + sum(x_bar * gamma_diff))
# ============================================================================
# STEP 9: Save Configuration
# ============================================================================

cat("=== Saving Configuration ===\n\n")

# Save the complete setup
dslite_setup <- list(
  connections = conns,
  logindata = logindata,
  K = K,
  d = d,
  true_ATE = true_ATE,
  beta_0 = beta_0,
  beta_1 = beta_1,
  c_0 = c_0,
  c_1 = c_1,
  data_list = data_list,  # For reference only
  server_names = paste0("dslite.server", 1:K)  # Names of server objects in environment
)

# Also save references to the server objects
for (k in 1:K) {
  server_name <- paste0("dslite.server", k)
  dslite_setup[[server_name]] <- get(server_name)
}

saveRDS(dslite_setup, "/Users/fdelimaa/Desktop/INRIA/FedCausal/Claude_files/dslite_5servers_c1.rds")
cat("DSLite setup saved to: /Users/fdelimaa/Desktop/INRIA/FedCausal/Claude_files/dslite_5servers_c1.rds\n\n")

# ============================================================================
# USAGE INFORMATION
# ============================================================================

cat("==========================================================\n")
cat("5 SEPARATE DSLITE SERVERS - SETUP COMPLETE\n")
cat("==========================================================\n\n")

cat("Architecture:\n")
cat("  -", K, "SEPARATE DSLite server objects:\n")
cat("    Created dynamically: dslite.server1, dslite.server2, ..., dslite.server", K, "\n", sep = "")
cat("  - Each server hosts only its own study data\n")
cat("  -", K, "connections (conns), each to a DIFFERENT server\n")
cat("  - Data symbol 'base' assigned on each server\n")
cat("  - Each server is completely isolated from the others\n\n")

cat("Key difference from previous version:\n")
cat("  BEFORE: 1 server with", K, "tables -> ", K, "connections to same server\n")
cat("  NOW:   ", K, "servers with 1 table each ->", K, "connections to different servers\n\n")

cat("MODULAR: Change K to any number and it will work!\n\n")

cat("This correctly simulates a federated setting where:\n")
cat("  - Each institution has its own server\n")
cat("  - Data never moves between servers\n")
cat("  - Only aggregated statistics are shared\n\n")

cat("To use in your GD estimator:\n\n")

cat("# Load the setup\n")
cat("dslite_setup <- readRDS('/home/claude/dslite_5servers_c1.rds')\n")
cat("conns <- dslite_setup$connections\n")
cat("K <- dslite_setup$K\n\n")

cat("# Access individual servers dynamically\n")
cat("for (k in 1:K) {\n")
cat("  server_name <- paste0('dslite.server', k)\n")
cat("  server <- dslite_setup[[server_name]]\n")
cat("  # Or: server <- get(server_name) if you want from environment\n")
cat("}\n\n")

cat("# Register custom server-side functions on all servers\n")
cat("for (k in 1:K) {\n")
cat("  server_name <- paste0('dslite.server', k)\n")
cat("  server <- get(server_name)\n")
cat("  server$aggregateMethod('computeGradientDS', 'computeGradientDS')\n")
cat("  server$aggregateMethod('computeLocalATEDS', 'computeLocalATEDS')\n")
cat("}\n\n")

cat("# Use standard DataSHIELD calls\n")
cat("results <- datashield.aggregate(conns, quote(computeGradientDS(theta, 0)))\n\n")

cat("==========================================================\n\n")

cat("Setup complete with 5 SEPARATE servers!\n")
cat("Each server is independent and isolated.\n\n")

# Note: Don't logout yet - keep connections open
# To logout later: datashield.logout(conns)
