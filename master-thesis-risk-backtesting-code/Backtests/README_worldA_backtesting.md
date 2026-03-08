# World A backtesting code

This repository contains the World A simulation / forecasting / e-value / GREM(S0) code.

The main script is:

```r
source("worldA_backtesting.r")
```

The script only defines functions. It does not automatically run examples when sourced.
## Required packages

Install the packages once if needed:

```r
install.packages(c("rugarch", "fGarch", "dplyr"))
```

Then load the script:

```r
source("worldA_backtesting.r")
```

## Notes before running

- The examples below are intentionally small and are meant as quick checks before upload.
- For `Expectile`, the e-value construction uses a finite `a_lower`. This has to be chosen conservatively for the problem at hand.
- The script sets `set.seed(123)` at the top, and `run_many_A()` also resets seeds internally via `seed0 + r`.

## Example 1: one replication — OCE_Huber (oracle_true)

Fast first check for a single path.

```r
rep_oce <- run_one_rep_A(
  risk_type = "OCE_Huber",
  forecast_mode = "oracle_true",
  f = 0.9,
  dist_apply_to = "r_only",   # try also "all" for OCE risks
  N_mc = 2000
)

detect_time(rep_oce$GREM$M, c = 5)
```

## Example 2: one replication — Expectile (ml_t)

Important: `Expectile` needs a finite `a_lower` so that the denominator in `e_expectile()` stays positive.

```r
rep_ex <- run_one_rep_A(
  risk_type = "Expectile",
  forecast_mode = "ml_t",
  f = 1.0,
  dist_apply_to = "r_only",
  N_mc = 3000,
  tau_ex = 0.975,
  a_lower = -10
)


detect_time(rep_ex$GREM$M, c = 5)
```

## Example 3: small multi-rep run — OCE_Huber

Small repetition run for a single distortion level.

```r
res_oce_small <- run_many_A(
  R = 50,
  seed0 = 123,
  risk_type = "OCE_Huber",
  forecast_mode = "oracle_true",
  f_set = c(0.9),
  dist_apply_to = "r_only",
  N_mc = 2000,
  huber_c = 1.5
)

res_oce_small[["f_0.9"]]$det_summary
```

## Example 4: small multi-rep run — Expectile

Again, `a_lower` needs to be set explicitly.

```r
res_ex_small <- run_many_A(
  R = 50,
  seed0 = 123,
  risk_type = "Expectile",
  forecast_mode = "ml_t",
  f_set = c(0.9, 1.0, 1.1),
  dist_apply_to = "r_only",
  N_mc = 3000,
  tau_ex = 0.975,
  a_lower = -10
)
# Example: results for f = 1
res_ex_small[[2]]$det_summary
```

## Example 5: restricted grid run — OCE_Huber and Expectile

This is a compact test of the grid wrapper on two risk measures and two forecast modes.

```r
grid_two_risks <- run_grid_worldA_allrisks(
  R = 50,
  seed0 = 123,
  risk_types = c("OCE_Huber", "Expectile"),
  forecast_modes = c("oracle_true", "ml_t"),
  f_set = c(0.9),
  N_mc = 2000,
  huber_c = 1.5,
  tau_ex = 0.975,
  a_lower = -10,
  verbose_print = FALSE
)

head(grid_two_risks$summary)
```

