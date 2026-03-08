# Interval-wise VaR/ES backtesting code

This repository contains the simulation and backtesting code for the comparison of GREE, GREL, GREM, WEALTH, and COV.

The main script is:

```r
source("intervalwise_backtesting_methods.r")
```

The script only defines functions. It does not automatically run examples when sourced.

## Required packages

Install the packages once if needed:

```r
install.packages(c("rugarch", "fGarch", "ggplot2", "dplyr", "tidyr"))
```

Then load the script:

```r
source("intervalwise_backtesting_methods.r")
```

## What the script covers

The script is long because it combines several simulation setups and several betting strategies in one place.

### Worlds

- `world = "A"`: stationary AR(1)-GARCH(1,1) losses with skewed-t innovations
- `world = "B"`: structural break in GARCH persistence during the backtest
- `world = "A7"`, `"B7"`, `"C7"`: direct-forecast Example-7 worlds from Wang, Wang, Ziegel (2025)
- `world = "D"`: regime-switching stress test built from A7 / B7 / C7 segments

### Risks

- `risk = "VaR"`
- `risk = "ES"`

### Forecast modes

- `forecast_mode = "oracle_true"`
- `forecast_mode = "oracle_empirical"`
- `forecast_mode = "ml_norm"`
- `forecast_mode = "ml_t"`
- `forecast_mode = "direct_world"` for worlds that already contain direct forecast paths (`A7`, `B7`, `C7`, `D`)

Notes:
- For **World B**, the code uses the Gaussian-QML and empirical-constants construction internally.
- For **World A**, `oracle_true` computes the ES constant by Monte Carlo once and caches it. For quick checks, use a smaller `oracle_mc_N`.
- For **World D**, `direct_world` uses the forecast formula implied by the segment design.

### Betting strategies returned by `run_one_rep()` / `run_many()`

- `GREE`
- `GREL`
- `GREM`
- `WEALTH`
- `COV`

### Interval designs

- `preset = "S0"`: one rolling interval
- `preset = "S2"`: recent short blocks plus a few older blocks
- `preset = "S3"`: recent / middle / early split
- `interval_mode = "multi_window"`: rolling windows of multiple lengths

### Distortions

For VaR:
- `distort_mode_var = "none"`
- `distort_mode_var = "scale"`
- `distort_mode_var = "smooth_step"`

For ES:
- `distort_mode_es = "none"`
- `distort_mode_es = "both"`
- `distort_mode_es = "es_only"`
- `distort_mode_es = "smooth_step"`

The smooth-step option is the banking-game transition. You control it via:
- `step_f_low`
- `step_f_high`
- `step_center`
- `step_width`
- `step_apply_to` (for ES: `"both"` or `"es_only"`)

### Guard presets

Two guard presets are defined in the script:

```r
var_guard_strict
var_guard_smooth
```

These are mainly relevant for interval-wise modes such as `S2`, `S3`, or `multi_window`.

## Notes before running

- The examples below are intentionally small and are meant as quick checks before upload.
- `run_many()` can still take time even with small `R`, especially when you use `oracle_true` for World A.
- The script sets `set.seed(123)` at the top, and `run_many()` resets seeds internally via `seed0 + r`.
- `run_many()` returns one matrix of log e-process paths per method in `res$M_store` and one matrix of detection times per method in `res$det_store`.
- A compact summary table is obtained by:
  ```r
  summarize_detection(res$det_store, T = T_backtest)
  ```


## Example 1: one World A run — ES, oracle_true, S0

A fast first check for a single ES path under the stationary setup.

```r
ex_a_es <- run_one_rep(
  world = "A",
  risk = "ES",
  f = 0.9,
  distort_mode_es = "both",
  forecast_mode = "oracle_true",
  preset = "S0",
  weighting = "size_recency",
  tau = 150,
  oracle_mc_N = 20000L
)

detect_time(ex_a_es$GREM$M, c = 5)
plot_one_path_loss_forecast(ex_a_es, risk = "ES")
```

## Example 2: World A banking game — ES, S2, strict guard

This mirrors the smooth-step distortion setup in a small test run.  
To get a decreasing multiplier from 1.10 to 0.90, set `step_f_low = 1.10` and `step_f_high = 0.90`.

```r
res_bank <- run_many(
  R = 10,
  seed0 = 123,
  world = "A",
  risk = "ES",
  f = 1,
  distort_mode_es = "smooth_step",
  step_f_low = 1.10,
  step_f_high = 0.90,
  step_center = 200,
  step_width = 30,
  step_apply_to = "both",
  forecast_mode = "oracle_true",
  preset = "S2",
  weighting = "size_recency",
  tau = 150,
  var_guard = var_guard_strict,
  oracle_mc_N = 20000L
)

summarize_detection(res_bank$det_store, T = T_backtest)
```


## Example 3: World B regime break — VaR, S2

A small structural-break check.  
Here `b_star` fixes the break point inside the backtest window.

```r
res_break <- run_many(
  R = 10,
  seed0 = 200,
  world = "B",
  risk = "VaR",
  b_star = 150,
  f = 1,
  distort_mode_var = "none",
  preset = "S2",
  weighting = "size_recency",
  tau = 150,
  var_guard = var_guard_strict
)

summarize_detection(res_break$det_store, T = T_backtest)
sum_br <- summarize_mean_paths(res_break$M_store)
plot_avg_log_eprocess(sum_br, title = "World B", ribbon = FALSE)
```

## Example 4: World B7 — direct forecasts, ES, S2

This is a compact direct-world check for one of the Example-7, Wang, Wang, Ziegel (2025) setups. 
Note that even if f=1, distortion mode=none, the forecast is underreporting since one uses a "direct" misspecified forecast mode.

```r
res_b7 <- run_many(
  R = 10,
  seed0 = 321,
  world = "B7",
  risk = "ES",
  f = 1,
  distort_mode_es = "none",
  forecast_mode = "direct_world",
  preset = "S2",
  weighting = "size_recency",
  tau = 150,
  var_guard = var_guard_strict,
  ex7_theta = ex7_theta
)

summarize_detection(res_b7$det_store, T = T_backtest)
```

## Example 5: World D regime switch — A7 to C7, direct forecasts, ES, S2

This gives a small regime-switching test using a fixed forecast formula from the first segment.


```r
segments_demo <- list(
  list(world = "A7", t_end = 700),
  list(world = "C7", t_end = 1000)
)

res_d <- run_many(
  R = 10,
  seed0 = 500,
  world = "D",
  risk = "ES",
  f = 1,
  distort_mode_es = "none",
  forecast_mode = "direct_world", 
  preset = "S2",
  weighting = "size_recency",
  tau = 150,
  var_guard = var_guard_strict,
  segments = segments_demo
)

summarize_detection(res_d$det_store, T = T_backtest)
```


## Useful objects returned by the main functions

### `run_one_rep(...)`

Returns a list containing, among other entries:

- `GREE`, `GREL`, `GREM`, `WEALTH`, `COV`
- `series` with the realised losses and reported forecast path(s)
- `meta` / `forecast_mode` information

Examples:

```r
ex_a_es$GREE$M
ex_a_es$GREM$lambda
ex_a_es$series$ES
```

### `run_many(...)`

Returns:

- `M_store`: one `(R x (T_backtest + 1))` matrix per method
- `det_store`: one `(R x 3)` detection-time matrix per method for thresholds `c = 2, 5, 10`
- `forecast_avg`: average point-forecast paths across successful runs
- `forecast_scalar`: scalar averages of the average forecast paths
- `settings`: the run settings


