# Systemic risk e-backtesting

This repository contains the script `systemic_risk_ebacktesting.r` for simulation, forecasting, and e-backtesting of systemic risk measures.

## What the script does

- simulates correlated AR(1)-GARCH(1,1) return series with Student-t innovations,
- backtests on losses,
- supports three forecast modes:
  - `oracle_true`
  - `oracle_empirical`
  - `ml_t`
- supports three systemic tests:
  - `MES`
  - `CoVaR`
  - `CoES`
- runs componentwise GREM,
- combines components via Bonferroni and mixture procedures.

The script only defines functions. It does not automatically run examples when sourced.

## Package requirement

Install `rugarch` before running the code:

```r
install.packages("rugarch")
```

Then load the script:

```r
source("systemic_risk_ebacktesting.r")
```

## Example runs


### Example 1: single replication

Choose one test and one forecast mode, then you get the detection result.

```r
source("systemic_risk_ebacktesting.r")

ex1 <- run_one_replication(
  test = "CoES",
  mode = "ml_t",
  f = 0.7,
  apply_to = "r_only",
  delta = 0.2,
  N_mc = 1500L,
  N_mc_cond = 1500L,
  refit_every = 10L,
  seed = 2026
)

cat("\n[Example 1] Single run summary:\n")
print(ex1$meta)
cat("Bonferroni detected:", ex1$results$bonferroni$detected, " at t=", ex1$results$bonferroni$t, "\n")
cat("Triggered component(s):", paste(ex1$results$bonferroni$who, collapse = ","), "\n")
cat("Mixture detected:", ex1$results$mixture$detected, " at t=", ex1$results$mixture$t, "\n")
```

### Example 2: small R with three distortions

Keep `R` small for a quick check.

```r
source("systemic_risk_ebacktesting.r")

ex2 <- run_many_replications(
  R = 30,
  seed0 = 1000,
  test = "MES",
  mode = "oracle_true",
  f_set = c(0.9, 1.0, 1.1),
  apply_to = "r_only",
  delta = 0.2,
  N_mc = 1500L,
  N_mc_cond = 1500L,
  refit_every = 10L,
  sim_rho = 0.6,
  sim_df  = 4
)
```

