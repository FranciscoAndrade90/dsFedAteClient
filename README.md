
# dsFedAteClient

<!-- badges: start -->
<!-- badges: end -->

Client-side 'DataSHIELD' package for federated estimation of the average
treatment effect (ATE) across multiple independent studies, each held on
a separate data server. `dsFedAteClient` runs on the analyst's machine;
the companion server-side package `dsFedAte` must be installed on each
data server.

## Installation

```r
install.packages("dsFedAteClient")
```

## Example

```r
library(dsFedAteClient)

conns <- DSI::datashield.login(logins = logindata, assign = TRUE,
                               symbol = "base")

result <- ds.EstimateAte(
  datasources     = conns,
  base_name       = "base",
  outcome_name    = "Y",
  treatment_name  = "W",
  covariate_names = c("X1", "X2", "X3"),
  heterogeneous   = TRUE
)

result$ate
result$ci_lower; result$ci_upper
```

A full worked example — including arm-specific variance and study-effects
adjustment, run locally with `DSLite` — is provided in the vignette:

```r
vignette("federated_rct", package = "dsFedAteClient")
```
