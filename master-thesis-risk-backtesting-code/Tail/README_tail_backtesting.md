# Tail risk e-backtesting

This repository contains the script `tail_risk_ebacktesting.r` for simulation, forecasting, and sequential e-backtesting of the tail-risk measures.

## What the script does

- simulates a univariate AR(1)-GARCH(1,1) loss process with standardized Student-t innovations,
- constructs forecasts for tail-risk functionals,
- backtests the four tail-risk targets:
  - `ES`
  - `TailSD`
  - `MedianShortfall`
  - `RVaR`
- supports the three forecast modes:
  - `oracle_true`
  - `oracle_empirical`
  - `ml_t`
- applies scale distortions with
  - `apply_to = "none"`
  - `apply_to = "r_only"`
  - `apply_to = "both"`
- runs componentwise GREM,
- combines component e-processes via Bonferroni and mixture rules,
- returns both delta-based summaries and c-threshold summaries for `c = 2, 5, 10`.

## Package requirement

Install `rugarch` before running the code:

```r
install.packages("rugarch")
```

Then source the script:

```r
source("tail_risk_ebacktesting.r")
```

## Main functions

### `run_one_rep_tail()`

Runs one replication for a single tail-risk test and one forecast mode.


### `run_many_tail()`

Runs repeated replications for one test / mode combination over a set of distortion factors `f_set`.

Useful when you want a small local check before scaling up.

### `run_tail_thesis_grid()`

Runs the full grid on demand. 
By default it iterates over:

- tests: `ES`, `TailSD`, `RVaR`, `MedianShortfall`
- modes: `oracle_true`, `ml_t`, `oracle_empirical`
- scenarios:
  - size check: `f = 1.0`, `apply_to = "none"`
  - alternatives: `f in c(0.9, 1.1)`, `apply_to = "r_only"`
  - alternatives: `f in c(0.9, 1.1)`, `apply_to = "both"`

You can also pass `save_path = "...rds"` if you want the returned object saved automatically.

## Example runs

The examples below are small so that they are practica checks.

### Example 1: single replication for TailSD under the student-t forecast mode

```r
source("tail_risk_ebacktesting.r")

ex1 <- run_one_rep_tail(
  test = "TailSD",
  mode = "ml_t",
  f = 0.7,
  apply_to = "both",
  delta = 0.05,
  N_mc = 1000L,
  N_cond = 1000L,
  refit_every = 10L,
  seed = 2026
)

cat("\n[Example 1] Meta settings:\n")
print(ex1$meta)

cat("\nBonferroni detection time:", ex1$results$bonferroni$t, "\n")
cat("Bonferroni triggering component(s):",
    paste(ex1$results$bonferroni$who, collapse = ", "), "\n")
cat("Mixture detection time:", ex1$results$mixture$t, "\n")
cat("Mixture detection time at c = 10:", ex1$results$mix_c$c10$t, "\n")
```


### Example 2: small power-oriented RVaR run

```r
source("tail_risk_ebacktesting.r")

pow_rvar <- run_many_tail(
  R = 10L,
  seed0 = 2000L,
  test = "RVaR",
  mode = "ml_t",
  f_set = c(0.9, 1.1),
  apply_to = "both",
  delta = 0.05,
  N_mc = 1000L,
  N_cond = 1000L,
  refit_every = 10L,
  verbose = FALSE
)

pow_rvar[["f_0.9"]]$mix_c$c10$summary
pow_rvar[["f_1.1"]]$mix_c$c10$summary
```

### Example 3: mini grid

This reproduces the grid logic on a much smaller scale.

```r
source("tail_risk_ebacktesting.r")

grid_small <- run_tail_thesis_grid(
  R = 5L,
  tests_grid = c("ES", "RVaR"),
  modes_grid = c("oracle_true", "ml_t"),
  N_mc = 800L,
  N_cond = 800L,
  refit_every = 10L,
  verbose = FALSE
)

names(grid_small)
grid_small[["RVaR__ml_t__ALT_f09_f11_both"]][["f_0.9"]]$mix_c$c2$summary
```
