###############################################################################
# VaR/ES backtesting betting methods
#
# This script collects the simulation and backtesting setups used to compare
# GREE, GREL, GREM, WEALTH, and COV under different interval constructions.
#
# Included worlds
# - World A: stationary AR(1)-GARCH(1,1) with skewed-t innovations (see Wang, Wang, Ziegel (2025), Section 7.1)
# - World B: structural break in GARCH persistence (see Wang, Wang, Ziegel (2025), Section 7.2)
# - Worlds A7 / B7 / C7: direct-forecast example worlds (see Wang, Wang, Ziegel (2025), Example 7)
# - World D: regime-switching combination of A7 / B7 / C7 segments
#
# Included features
# - VaR and ES backtesting
# - oracle / empirical / ML forecast construction where applicable
# - fixed-scale and smooth-step distortions
# - interval presets S0 / S2 / S3 and a multi-window design
# - repeated replications, detection summaries, and plotting helpers
###############################################################################

suppressPackageStartupMessages({
  library(rugarch)  # GARCH fitting and filtering
  library(fGarch)   # skewed-t simulation and quantiles
  library(ggplot2)  # plotting helpers
  library(dplyr)    # summaries and data support
  library(tidyr)    # plotting helpers
})

set.seed(123)  # seed for reproducibility

# ============================================================
# 0) Global settings
# ============================================================

# Tail levels for the VaR / ES backtests
alpha_VaR <- 0.975
alpha_ES  <- 0.975

# Innovation parameters for World A
df_true   <- 5
skew_true <- 1.5

# Horizon split
T_train    <- 500
T_backtest <- 500
T_total    <- T_train + T_backtest

# Example-7 defaults
# ex7_theta = 4*pi/T_backtest gives two complete cycles over the 500-day backtest.
ex7_theta <- 4 * pi / T_backtest
ex7_noise_set <- c(0, as.vector(sapply(1:5, function(i) c(-i, i))) / 10)

# ============================================================
# 1) Simulation worlds
# ============================================================

# ------------------------------------------------------------
# World A
# Stationary AR(1) mean + GARCH(1,1) volatility with skewed-t innovations.
# Returns losses L_t, the true conditional mean mu_t, conditional sd sig_t,
# and the underlying innovation series Z_t.
# ------------------------------------------------------------

simulate_world_A <- function(T_total, df = df_true, skew = skew_true, burnin = 1000L) {
  T_all <- T_total + burnin

  L   <- numeric(T_all)
  mu  <- numeric(T_all)
  sig <- numeric(T_all)
  Z   <- numeric(T_all)

  # GARCH parameters
  g_omega <- 0.01
  g_alpha <- 0.10
  g_beta  <- 0.85

  # AR(1) mean parameters
  mu0 <- -0.05
  phi <-  0.30

  # Start at the unconditional volatility level
  sig[1] <- sqrt(g_omega / (1 - g_alpha - g_beta))
  mu[1]  <- mu0
  Z[1]   <- rsstd(1, mean = 0, sd = 1, nu = df, xi = skew)
  L[1]   <- mu[1] + sig[1] * Z[1]

  for (t in 2:T_all) {
    mu[t]  <- mu0 + phi * L[t - 1]
    sig[t] <- sqrt(g_omega + g_alpha * sig[t - 1]^2 * Z[t - 1]^2 + g_beta * sig[t - 1]^2)
    Z[t]   <- rsstd(1, mean = 0, sd = 1, nu = df, xi = skew)
    L[t]   <- mu[t] + sig[t] * Z[t]
  }

  idx <- (burnin + 1L):T_all

  list(
    L   = L[idx],
    mu  = mu[idx],
    sig = sig[idx],
    Z   = Z[idx]
  )
}

# ------------------------------------------------------------
# World B
# Structural break in GARCH persistence during the backtest window.
# The break location b_star is measured in backtest days:
#   b_star = 0   means break right at the start of backtest,
#   b_star = 200 means break after 200 backtest observations, etc.
# ------------------------------------------------------------

simulate_world_B <- function(T_total,
                             T_train_in    = T_train,
                             T_backtest_in = T_backtest,
                             b_star = 200,
                             burnin = 1000L,
                             omega = 0.00001,
                             alpha = 0.04,
                             beta0 = 0.70,
                             beta_jump = 0.251,
                             nu = 5,
                             skew = 0.95) {

  T_total    <- as.integer(T_total)
  burnin     <- as.integer(burnin)
  T_train    <- as.integer(T_train_in)
  T_backtest <- as.integer(T_backtest_in)

  if (is.null(b_star)) {
    b_star <- sample.int(T_backtest, 1) - 1L
  } else {
    b_star <- as.integer(b_star)
    if (b_star < 0L || b_star > (T_backtest - 1L)) {
      stop("b_star must be in {0, ..., T_backtest - 1}. Got: ", b_star)
    }
  }

  # Backtest starts at t = T_train + 1 in the full series.
  break_t_full <- as.integer(T_train + b_star)

  T_all <- T_total + burnin

  L      <- numeric(T_all)
  sig    <- numeric(T_all)
  Z      <- numeric(T_all)
  beta_t <- numeric(T_all)

  beta_t[1] <- beta0
  if (alpha + beta0 >= 1) stop("Need alpha + beta0 < 1 for unconditional initialisation.")

  sig[1] <- sqrt(omega / (1 - alpha - beta0))
  Z[1]   <- rsstd(1, mean = 0, sd = 1, nu = nu, xi = skew)
  L[1]   <- -sig[1] * Z[1]

  for (t in 2:T_all) {
    in_post_break <- (t > (burnin + break_t_full))
    beta_t[t] <- beta0 + beta_jump * as.numeric(in_post_break)

    sig2 <- omega + alpha * (L[t - 1]^2) + beta_t[t] * (sig[t - 1]^2)
    sig[t] <- sqrt(max(sig2, 0))

    Z[t] <- rsstd(1, mean = 0, sd = 1, nu = nu, xi = skew)
    L[t] <- -sig[t] * Z[t]
  }

  idx <- (burnin + 1L):T_all

  list(
    L = L[idx],
    mu = rep(0, T_total),
    sig = sig[idx],
    Z = Z[idx],
    b_star = b_star,
    break_t_full = break_t_full,
    beta_t = beta_t[idx]
  )
}

# ------------------------------------------------------------
# Example-7 worlds A7 / B7 / C7
# These worlds come with direct forecast paths VaR_full / ES_full.
# ------------------------------------------------------------

simulate_world_A7 <- function(T_total,
                              alpha = 0.95,
                              z_const = 1.48,
                              r_const = 1.86) {
  Z <- rnorm(T_total)
  t <- 1:T_total
  s <- 1 + t / T_total

  L   <- s * Z
  mu  <- rep(0, T_total)
  sig <- s

  VaR_full <- z_const * s
  ES_full  <- r_const * s

  list(
    L = L,
    mu = mu,
    sig = sig,
    Z = Z,
    VaR_full = VaR_full,
    ES_full  = ES_full,
    meta = list(type = "A7_WWZ_linear", alpha = alpha, z_const = z_const, r_const = r_const)
  )
}

simulate_world_B7 <- function(T_total,
                              alpha = 0.95,
                              theta = ex7_theta,
                              z_const = 1.48,
                              r_const = 1.86) {
  Z <- rnorm(T_total)
  t <- 1:T_total
  s <- 1 + sin(theta * t)

  L   <- s * Z
  mu  <- rep(0, T_total)
  sig <- s

  VaR_full <- z_const * s
  ES_full  <- r_const * s

  list(
    L = L,
    mu = mu,
    sig = sig,
    Z = Z,
    VaR_full = VaR_full,
    ES_full  = ES_full,
    meta = list(type = "B7_WWZ_sine", alpha = alpha, theta = theta, z_const = z_const, r_const = r_const)
  )
}

simulate_world_C7 <- function(T_total,
                              alpha = 0.95,
                              z0 = 1.64,
                              r0 = 2.06,
                              eps_set = c(0, as.vector(sapply(1:5, function(i) c(-i, i))) / 10)) {
  Z <- rnorm(T_total)

  L   <- Z
  mu  <- rep(0, T_total)
  sig <- rep(1, T_total)

  eps <- sample(eps_set, size = T_total, replace = TRUE)
  VaR_full <- z0 + eps
  ES_full  <- r0 + eps

  list(
    L = L,
    mu = mu,
    sig = sig,
    Z = Z,
    VaR_full = VaR_full,
    ES_full  = ES_full,
    eps = eps,
    meta = list(type = "C7_WWZ_iid_noisy", alpha = alpha, z0 = z0, r0 = r0, eps_set = eps_set)
  )
}

# ------------------------------------------------------------
# World D
# Piecewise combination of A7 / B7 / C7 over the full horizon.
# If stick_forecast_to_first = TRUE, the reported forecasts keep following the
# formula of the first segment even after the data regime has changed.
# ------------------------------------------------------------

simulate_world_D <- function(T_total,
                             segments,
                             alpha = 0.95,
                             theta = ex7_theta,
                             z_const = 1.48,
                             r_const = 1.86,
                             z0 = 1.64,
                             r0 = 2.06,
                             eps_set = c(0, as.vector(sapply(1:5, function(i) c(-i, i))) / 10),
                             stick_forecast_to_first = TRUE) {

  if (!is.list(segments) || length(segments) < 1L) {
    stop("segments must be a non-empty list.")
  }

  t_ends <- vapply(segments, function(s) as.integer(s$t_end), integer(1))
  if (max(t_ends) != T_total) stop("Last segment t_end must equal T_total.")
  if (any(diff(t_ends) <= 0)) stop("t_end values must be strictly increasing.")

  worlds <- vapply(segments, function(s) as.character(s$world), character(1))

  Z <- rnorm(T_total)
  t <- 1:T_total

  # For C7-type forecasts, keep one noise sequence over the full horizon.
  eps_full <- sample(eps_set, size = T_total, replace = TRUE)

  L <- numeric(T_total)
  mu <- rep(0, T_total)
  sig <- numeric(T_total)
  VaR_full <- numeric(T_total)
  ES_full  <- numeric(T_total)
  regime <- character(T_total)

  forecast_world <- if (stick_forecast_to_first) worlds[1] else NA_character_

  forecast_formula <- function(W, idx) {
    if (W == "A7") {
      s <- 1 + t[idx] / T_total
      return(list(VaR = z_const * s, ES = r_const * s))
    }
    if (W == "B7") {
      s <- 1 + sin(theta * t[idx])
      return(list(VaR = z_const * s, ES = r_const * s))
    }
    if (W == "C7") {
      return(list(VaR = z0 + eps_full[idx], ES = r0 + eps_full[idx]))
    }
    stop("Unknown world in forecast_formula: ", W)
  }

  t_start <- 1L
  for (j in seq_along(segments)) {
    t_end <- t_ends[j]
    idx <- t_start:t_end
    w_data <- worlds[j]

    if (w_data == "A7") {
      s <- 1 + t[idx] / T_total
      L[idx] <- s * Z[idx]
      sig[idx] <- s
      regime[idx] <- "A7"
    } else if (w_data == "B7") {
      s <- 1 + sin(theta * t[idx])
      L[idx] <- s * Z[idx]
      sig[idx] <- s
      regime[idx] <- "B7"
    } else if (w_data == "C7") {
      L[idx] <- Z[idx]
      sig[idx] <- 1
      regime[idx] <- "C7"
    } else {
      stop("Unknown data world: ", w_data)
    }

    w_fc <- if (stick_forecast_to_first) forecast_world else w_data
    fc <- forecast_formula(w_fc, idx)
    VaR_full[idx] <- fc$VaR
    ES_full[idx]  <- fc$ES

    t_start <- t_end + 1L
  }

  list(
    L = L,
    mu = mu,
    sig = sig,
    Z = Z,
    VaR_full = VaR_full,
    ES_full  = ES_full,
    eps = eps_full,
    regime = regime,
    meta = list(
      type = "D_WWZ_combo",
      alpha = alpha,
      theta = theta,
      segments = segments,
      stick_forecast_to_first = stick_forecast_to_first,
      forecast_world = forecast_world
    )
  )
}

# ------------------------------------------------------------
# Extract the backtest slice t = T_train + 1, ..., T_train + T_backtest.
# If a world provides direct forecast paths, include them as well.
# ------------------------------------------------------------

get_backtest_slice <- function(sim, T_train, T_backtest) {
  idx_bt <- (T_train + 1L):(T_train + T_backtest)

  out <- list(
    L_back   = sim$L[idx_bt],
    mu_back  = sim$mu[idx_bt],
    sig_back = sim$sig[idx_bt],
    idx_bt   = idx_bt
  )

  if (!is.null(sim$VaR_full)) out$VaR_back <- sim$VaR_full[idx_bt]
  if (!is.null(sim$ES_full))  out$ES_back  <- sim$ES_full[idx_bt]

  out
}

# ============================================================
# 2) Forecast helpers
# ============================================================

# ------------------------------------------------------------
# Empirical innovation constants from a presample of standardized residuals.
# These are used for the empirical forecast mode.
# ------------------------------------------------------------

estimate_innov_constants_presample <- function(L_full, mu_true, sig_true,
                                               alpha_VaR, alpha_ES,
                                               presample_length,
                                               q_type = 8) {
  stopifnot(length(L_full)  >= presample_length)
  stopifnot(length(mu_true) >= presample_length)
  stopifnot(length(sig_true) >= presample_length)

  idx <- 1:presample_length
  Z_hat <- (L_full[idx] - mu_true[idx]) / sig_true[idx]

  q_hat <- as.numeric(quantile(Z_hat, probs = alpha_VaR, type = q_type))

  q_es  <- as.numeric(quantile(Z_hat, probs = alpha_ES, type = q_type))
  tailZ <- Z_hat[Z_hat > q_es]
  es_hat <- if (length(tailZ) == 0L) q_es else mean(tailZ)

  list(q_hat = q_hat, es_hat = es_hat)
}

build_baseline_forecasts <- function(mu_true, sig_true,
                                     T_train, T_backtest,
                                     q_hat, es_hat) {
  idx_bt <- (T_train + 1L):(T_train + T_backtest)

  VaR_base <- mu_true[idx_bt] + sig_true[idx_bt] * q_hat
  ES_base  <- mu_true[idx_bt] + sig_true[idx_bt] * es_hat

  list(VaR = VaR_base, ES = ES_base)
}

make_empirical_baseline_forecasts <- function(sim,
                                              alpha_VaR, alpha_ES,
                                              T_train, T_backtest,
                                              q_type = 8) {
  presample_length <- T_train

  consts <- estimate_innov_constants_presample(
    L_full = sim$L,
    mu_true = sim$mu,
    sig_true = sim$sig,
    alpha_VaR = alpha_VaR,
    alpha_ES  = alpha_ES,
    presample_length = presample_length,
    q_type = q_type
  )

  fc <- build_baseline_forecasts(
    mu_true = sim$mu,
    sig_true = sim$sig,
    T_train = T_train,
    T_backtest = T_backtest,
    q_hat = consts$q_hat,
    es_hat = consts$es_hat
  )

  list(
    VaR = fc$VaR,
    ES  = fc$ES,
    q_hat = consts$q_hat,
    es_hat = consts$es_hat,
    presample_length = presample_length
  )
}

# ------------------------------------------------------------
# World A forecast modes
# - oracle_true: uses true conditional mean / scale and skewed-t constants
# - oracle_empirical: uses presample empirical constants
# - ml_norm / ml_t: rolling AR(1)-GARCH(1,1) fits
# ------------------------------------------------------------

.oracle_cache <- new.env(parent = emptyenv())

oracle_true_constants_sstd <- function(alpha, df, skew, N_mc = 200000L) {
  key <- paste0("a=", alpha, "|df=", df, "|sk=", skew, "|N=", N_mc)
  if (exists(key, envir = .oracle_cache, inherits = FALSE)) {
    return(get(key, envir = .oracle_cache, inherits = FALSE))
  }

  q <- as.numeric(fGarch::qsstd(alpha, mean = 0, sd = 1, nu = df, xi = skew))

  Z <- fGarch::rsstd(N_mc, mean = 0, sd = 1, nu = df, xi = skew)
  tailZ <- Z[Z > q]
  es <- if (length(tailZ) == 0L) q else mean(tailZ)

  out <- list(q = q, es = es, N_mc = N_mc)
  assign(key, out, envir = .oracle_cache)
  out
}

# Analytic constants for N(0, 1)
const_norm <- function(alpha) {
  q <- qnorm(alpha)
  es <- dnorm(q) / (1 - alpha)
  list(q = q, es = es)
}

# Analytic constants for standardized t innovations with variance 1.
const_std <- function(alpha, nu) {
  nu <- as.numeric(nu)

  if (!is.finite(nu) || nu <= 2) {
    nu2 <- max(nu, 2.1)
    q <- qt(alpha, df = nu2)
    esT <- (dt(q, df = nu2) / (1 - alpha)) * ((nu2 + q^2) / (nu2 - 1))
    return(list(q = q, es = esT, nu = nu, note = "nu<=2 fallback"))
  }

  s <- sqrt((nu - 2) / nu)
  qT <- qt(alpha, df = nu)
  esT <- (dt(qT, df = nu) / (1 - alpha)) * ((nu + qT^2) / (nu - 1))

  list(q = s * qT, es = s * esT, nu = nu)
}

fit_ar1_garch11 <- function(x, dist = c("norm", "std")) {
  dist <- match.arg(dist)

  spec <- rugarch::ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model     = list(armaOrder = c(1, 0), include.mean = TRUE),
    distribution.model = dist
  )

  rugarch::ugarchfit(spec = spec, data = x, solver = "hybrid")
}

# Rolling ML forecasts with blockwise refits.
# Within a refit block, parameters stay fixed and the recursion is updated
# using the realised losses.
build_ml_forecasts_worldA <- function(L_full, T_train, T_backtest,
                                      dist = c("norm", "std"),
                                      refit_every = 50L) {
  dist <- match.arg(dist)
  refit_every <- as.integer(refit_every)
  if (refit_every <= 0L) stop("refit_every must be positive.")

  VaR <- numeric(T_backtest)
  ES  <- numeric(T_backtest)

  get_consts_pair <- function(fit) {
    cf <- coef(fit)

    if (dist == "norm") {
      cV <- const_norm(alpha_VaR)
      cE <- const_norm(alpha_ES)
      return(list(qV = cV$q, esE = cE$es, dist = "norm", nu = NA_real_))
    }

    nu <- if ("shape" %in% names(cf)) as.numeric(cf[["shape"]]) else NA_real_
    cV <- const_std(alpha_VaR, nu)
    cE <- const_std(alpha_ES,  nu)

    list(qV = cV$q, esE = cE$es, dist = "std", nu = nu)
  }

  k <- 1L
  while (k <= T_backtest) {
    t_prev_full <- T_train + (k - 1L)

    win_start <- t_prev_full - T_train + 1L
    win_end   <- t_prev_full
    x_win <- L_full[win_start:win_end]

    fit <- tryCatch(
      fit_ar1_garch11(x_win, dist = dist),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      VaR[k:T_backtest] <- NA_real_
      ES[k:T_backtest]  <- NA_real_
      break
    }

    cf <- coef(fit)
    mu0 <- if ("mu" %in% names(cf)) as.numeric(cf[["mu"]]) else 0

    ar1_name <- intersect(names(cf), c("ar1", "ar[1]", "ar"))
    ar1 <- if (length(ar1_name)) as.numeric(cf[[ar1_name[1]]]) else 0

    omg <- if ("omega"  %in% names(cf)) as.numeric(cf[["omega"]])  else NA_real_
    a1  <- if ("alpha1" %in% names(cf)) as.numeric(cf[["alpha1"]]) else NA_real_
    b1  <- if ("beta1"  %in% names(cf)) as.numeric(cf[["beta1"]])  else NA_real_

    sig_series <- as.numeric(rugarch::sigma(fit))
    mu_series  <- as.numeric(rugarch::fitted(fit))

    L_prev   <- L_full[t_prev_full]
    sig_prev <- sig_series[length(sig_series)]
    mu_prev  <- mu_series[length(mu_series)]
    res_prev <- L_prev - mu_prev

    cc <- get_consts_pair(fit)
    qV  <- cc$qV
    esE <- cc$esE

    block_len <- min(refit_every, T_backtest - k + 1L)

    for (j in 1:block_len) {
      t_full <- t_prev_full + 1L

      mu_f  <- mu0 + ar1 * L_prev
      sig_f <- sqrt(pmax(omg + a1 * (res_prev^2) + b1 * (sig_prev^2), 0))

      VaR[k] <- mu_f + sig_f * qV
      ES[k]  <- mu_f + sig_f * esE

      L_prev   <- L_full[t_full]
      mu_prev  <- mu_f
      sig_prev <- sig_f
      res_prev <- L_prev - mu_prev

      t_prev_full <- t_full
      k <- k + 1L
      if (k > T_backtest) break
    }
  }

  list(VaR = VaR, ES = ES)
}

build_forecasts_worldA <- function(sim,
                                   forecast_mode = c("oracle_empirical", "oracle_true", "ml_norm", "ml_t"),
                                   T_train, T_backtest,
                                   alpha_VaR, alpha_ES,
                                   q_type = 8,
                                   refit_every = 50L,
                                   oracle_mc_N = 200000L) {
  forecast_mode <- match.arg(forecast_mode)

  if (forecast_mode == "oracle_empirical") {
    fc <- make_empirical_baseline_forecasts(
      sim = sim,
      alpha_VaR = alpha_VaR,
      alpha_ES = alpha_ES,
      T_train = T_train,
      T_backtest = T_backtest,
      q_type = q_type
    )
    return(list(VaR = fc$VaR, ES = fc$ES, mode = "oracle_empirical"))
  }

  if (forecast_mode == "oracle_true") {
    const <- oracle_true_constants_sstd(
      alpha = alpha_ES,
      df = df_true,
      skew = skew_true,
      N_mc = oracle_mc_N
    )

    qV <- as.numeric(fGarch::qsstd(alpha_VaR, mean = 0, sd = 1, nu = df_true, xi = skew_true))
    esE <- const$es

    idx_bt <- (T_train + 1L):(T_train + T_backtest)
    VaR <- sim$mu[idx_bt] + sim$sig[idx_bt] * qV
    ES  <- sim$mu[idx_bt] + sim$sig[idx_bt] * esE

    return(list(VaR = VaR, ES = ES, mode = "oracle_true", oracle_mc_N = oracle_mc_N))
  }

  if (forecast_mode == "ml_norm") {
    fc <- build_ml_forecasts_worldA(sim$L, T_train, T_backtest, dist = "norm", refit_every = refit_every)
    return(list(VaR = fc$VaR, ES = fc$ES, mode = "ml_norm", refit_every = refit_every))
  }

  if (forecast_mode == "ml_t") {
    fc <- build_ml_forecasts_worldA(sim$L, T_train, T_backtest, dist = "std", refit_every = refit_every)
    return(list(VaR = fc$VaR, ES = fc$ES, mode = "ml_t", refit_every = refit_every))
  }

  stop("Unknown forecast_mode: ", forecast_mode)
}

# ------------------------------------------------------------
# World B forecast construction
# Gaussian QML on the presample + empirical innovation constants,
# then fixed-parameter filtering over the full series.
# ------------------------------------------------------------

make_ec44_forecasts_worldB <- function(L_full,
                                       alpha_VaR, alpha_ES,
                                       presample_length,
                                       backtest_length,
                                       q_type = 8) {
  stopifnot(length(L_full) >= presample_length + backtest_length)

  x_pre <- L_full[1:presample_length]

  spec <- rugarch::ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model     = list(armaOrder = c(0, 0), include.mean = FALSE),
    distribution.model = "norm"
  )

  fit <- rugarch::ugarchfit(
    spec = spec,
    data = x_pre,
    solver = "hybrid"
  )

  cf <- coef(fit)
  fixed_list <- as.list(cf)

  spec_fix <- rugarch::ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model     = list(armaOrder = c(0, 0), include.mean = FALSE),
    distribution.model = "norm",
    fixed.pars = fixed_list
  )

  filt <- rugarch::ugarchfilter(spec = spec_fix, data = L_full)
  sig_hat_full <- as.numeric(rugarch::sigma(filt))

  z_hat_pre <- x_pre / sig_hat_full[1:presample_length]

  q_hat <- as.numeric(stats::quantile(z_hat_pre, probs = alpha_VaR, type = q_type))

  q_es  <- as.numeric(stats::quantile(z_hat_pre, probs = alpha_ES, type = q_type))
  tailZ <- z_hat_pre[z_hat_pre > q_es]
  es_hat <- if (length(tailZ) == 0L) q_es else mean(tailZ)

  idx_bt <- (presample_length + 1L):(presample_length + backtest_length)
  sig_bt <- sig_hat_full[idx_bt]

  VaR_bt <- sig_bt * q_hat
  ES_bt  <- sig_bt * es_hat

  list(
    VaR = VaR_bt,
    ES  = ES_bt,
    q_hat = q_hat,
    es_hat = es_hat,
    theta_hat = cf,
    sig_full = sig_hat_full
  )
}

# ============================================================
# 3) Distortions
# ============================================================

# Banking-game smooth transition from one multiplier to another.
smooth_step <- function(t, center, width) {
  1 / (1 + exp(-(t - center) / width))
}

# ------------------------------------------------------------
# VaR distortions
# - "none": keep the base forecast
# - "scale": multiply by a fixed factor f
# - "smooth_step": move gradually from f_low to f_high over time
# ------------------------------------------------------------

distort_var_forecast <- function(VaR_base,
                                 f = 1.0,
                                 mode = c("none", "scale", "smooth_step"),
                                 f_low = 0.90,
                                 f_high = 1.10,
                                 center = NULL,
                                 width = NULL) {
  mode <- match.arg(mode)
  T <- length(VaR_base)
  t <- 1:T

  if (mode == "none") {
    return(list(VaR = VaR_base, label = "none", mult = rep(1, T)))
  }

  if (mode == "scale") {
    return(list(
      VaR = f * VaR_base,
      label = sprintf("scale_%+g%%", (f - 1) * 100),
      mult = rep(f, T)
    ))
  }

  if (is.null(center)) center <- (T + 1) / 2
  if (is.null(width))  width  <- T / 10

  w <- smooth_step(t, center = center, width = width)
  mult <- f_low + (f_high - f_low) * w

  list(
    VaR = mult * VaR_base,
    label = sprintf("smooth_step_%g_to_%g_c%g_w%g", f_low, f_high, center, width),
    mult = mult
  )
}

# ------------------------------------------------------------
# ES distortions
# - "none": no distortion
# - "both": scale ES and VaR together
# - "es_only": scale only ES
# - "smooth_step": smooth transition over time, applied either to both
#                  components or to ES only
# ------------------------------------------------------------

distort_es_forecasts <- function(ES_base, VaR_base,
                                 f = 1.0,
                                 mode = c("none", "both", "es_only", "smooth_step"),
                                 f_low = 0.90,
                                 f_high = 1.10,
                                 center = NULL,
                                 width = NULL,
                                 apply_to = c("both", "es_only")) {
  mode <- match.arg(mode)
  apply_to <- match.arg(apply_to)

  T <- length(ES_base)
  t <- 1:T

  if (mode == "none") {
    return(list(
      ES = ES_base,
      VaR = VaR_base,
      label = "none",
      mult = rep(1, T)
    ))
  }

  if (mode == "both") {
    return(list(
      ES = f * ES_base,
      VaR = f * VaR_base,
      label = sprintf("both_%+g%%", (f - 1) * 100),
      mult = rep(f, T)
    ))
  }

  if (mode == "es_only") {
    return(list(
      ES = f * ES_base,
      VaR = VaR_base,
      label = sprintf("esonly_%+g%%", (f - 1) * 100),
      mult = rep(f, T)
    ))
  }

  if (is.null(center)) center <- (T + 1) / 2
  if (is.null(width))  width  <- T / 10

  w <- smooth_step(t, center = center, width = width)
  mult <- f_low + (f_high - f_low) * w

  ES  <- mult * ES_base
  VaR <- VaR_base
  if (apply_to == "both") VaR <- mult * VaR_base

  list(
    ES = ES,
    VaR = VaR,
    label = sprintf("smooth_step_%s_%g_to_%g_c%g_w%g", apply_to, f_low, f_high, center, width),
    mult = mult
  )
}

# ============================================================
# 4) E-values and e-process updates
# ============================================================

# VaR e-value:
# e_t = 1 / (1 - alpha) * 1{L_t > VaR_t}
e_value_var <- function(L_t, VaR_t, alpha_VaR) {
  (1 / (1 - alpha_VaR)) * as.numeric(L_t > VaR_t)
}

e_values_var <- function(L_back, VaR_back, alpha_VaR) {
  stopifnot(length(L_back) == length(VaR_back))
  (1 / (1 - alpha_VaR)) * as.numeric(L_back > VaR_back)
}

# ES e-value:
# e_t = (L_t - VaR_t)_+ / ((1 - alpha) * (ES_t - VaR_t))
e_value_es <- function(L_t, ES_t, VaR_t, alpha_ES) {
  den <- (1 - alpha_ES) * (ES_t - VaR_t)
  if (!is.finite(den) || den <= 0) return(Inf)
  pmax(L_t - VaR_t, 0) / den
}

e_values_es <- function(L_back, ES_back, VaR_back, alpha_ES) {
  stopifnot(length(L_back) == length(ES_back), length(L_back) == length(VaR_back))

  den <- (1 - alpha_ES) * (ES_back - VaR_back)
  num <- pmax(L_back - VaR_back, 0)

  out <- rep(Inf, length(L_back))
  ok <- is.finite(den) & (den > 0)
  out[ok] <- num[ok] / den[ok]

  out
}

# Sequential update of log-wealth.
update_e_process <- function(M_prev, e_t, lambda_t) {
  if (is.infinite(M_prev) && M_prev < 0) return(-Inf)
  if (!is.finite(M_prev)) return(Inf)

  if (!is.finite(lambda_t) || lambda_t < 0) lambda_t <- 0
  if (lambda_t == 0) return(M_prev)

  if (!is.finite(e_t)) return(Inf)

  inc <- 1 - lambda_t + lambda_t * e_t
  if (!is.finite(inc)) return(Inf)
  if (inc <= 0) return(-Inf)

  M_prev + log(inc)
}

# Stable weight for mixing two log-wealth processes.
weight_from_logs <- function(logA, logB) {
  if (is.infinite(logA) && logA < 0 && is.infinite(logB) && logB < 0) return(0.5)
  if (is.infinite(logA) && logA < 0) return(0)
  if (is.infinite(logB) && logB < 0) return(1)
  1 / (1 + exp(logB - logA))
}

# ============================================================
# 5) Interval construction and lambda rules
# ============================================================

# ------------------------------------------------------------
# Select historical intervals used to estimate candidate lambdas at time t.
# - S0: one rolling window
# - S2: many short recent blocks plus a few older blocks
# - S3: recent / middle / early split
# - multi_window: multiple rolling windows of different lengths
# ------------------------------------------------------------

select_intervals <- function(t,
                             mode = c("multi_window", "preset"),
                             J_vec = c(15L, 30L, 60L, 120L, 240L, 480L),
                             preset = c("S0", "S2", "S3"),
                             J_total = 500L,
                             J_recent = 100L,
                             J_mid = 150L) {
  mode <- match.arg(mode)
  if (t <= 1L) return(structure(list(integer(0)), past_end = 0L))

  past_end <- t - 1L

  safe_seq <- function(a, b) {
    a <- max(1L, as.integer(a))
    b <- min(past_end, as.integer(b))
    if (a > b) integer(0) else a:b
  }

  if (mode == "multi_window") {
    out <- list()
    for (J in J_vec) {
      out[[length(out) + 1L]] <- safe_seq(past_end - J + 1L, past_end)
    }
    return(structure(out, past_end = past_end))
  }

  preset <- match.arg(preset)

  if (preset == "S0") {
    out <- list(safe_seq(past_end - J_total + 1L, past_end))
    return(structure(out, past_end = past_end))
  }

  I_recent <- safe_seq(past_end - J_recent + 1L, past_end)

  if (preset == "S2") {
    J <- 10L
    recent_days <- 150L
    K_recent <- as.integer(ceiling(recent_days / J))
    old_ages <- c(200L, 300L, 450L)

    if (past_end <= J * K_recent) {
      ends <- seq(from = past_end, to = 1L, by = -J)
      out <- lapply(ends, function(e) safe_seq(e - J + 1L, e))
      out <- out[vapply(out, length, integer(1)) > 0]
      return(structure(out, past_end = past_end))
    }

    recent_ends <- past_end - J * (0:(K_recent - 1L))
    older_ends <- past_end - old_ages
    older_ends <- older_ends[older_ends >= 1L]

    ends_all <- unique(c(recent_ends, older_ends))
    ends_all <- sort(ends_all, decreasing = TRUE)

    out <- lapply(ends_all, function(e) safe_seq(e - J + 1L, e))
    out <- out[vapply(out, length, integer(1)) > 0]

    return(structure(out, past_end = past_end))
  }

  I_mid   <- safe_seq(past_end - (J_recent + J_mid) + 1L, past_end - J_recent)
  I_early <- safe_seq(past_end - J_total + 1L, past_end - (J_recent + J_mid))

  structure(list(I_recent, I_mid, I_early), past_end = past_end)
}

# Combine interval-wise lambda candidates into one lambda_t.
# The optional guard can fall back to a short recent prefix if the interval-wise
# lambda estimates disagree too much or show strong drift.
combine_interval_lambdas <- function(lambda_vec,
                                     interval_list,
                                     weighting = c("size", "equal", "size_recency"),
                                     tau = 150) {
  weighting <- match.arg(weighting)

  sizes <- vapply(interval_list, length, integer(1))
  ends  <- vapply(
    interval_list,
    function(idx) if (length(idx) == 0) NA_integer_ else max(idx),
    integer(1)
  )

  past_end <- attr(interval_list, "past_end")
  if (is.null(past_end) || !is.finite(past_end)) {
    past_end <- max(ends, na.rm = TRUE)
    if (!is.finite(past_end)) past_end <- 0L
  }

  ok <- (sizes > 0) & is.finite(lambda_vec)
  if (!any(ok)) return(0)

  vg <- attr(interval_list, "var_guard")
  if (!is.null(vg) && is.list(vg) && length(lambda_vec) >= 2) {
    diff_thresh <- if (!is.null(vg$diff_thresh)) as.numeric(vg$diff_thresh) else 0.20
    sd_thresh   <- if (!is.null(vg$sd_thresh))   as.numeric(vg$sd_thresh)   else 0.12

    compare_q       <- if (!is.null(vg$compare_q))       as.integer(vg$compare_q)       else 2L
    mean_gap_thresh <- if (!is.null(vg$mean_gap_thresh)) as.numeric(vg$mean_gap_thresh) else 0.10
    slope_thresh    <- if (!is.null(vg$slope_thresh))    as.numeric(vg$slope_thresh)    else 0.03

    lam <- lambda_vec
    lam[!ok] <- NA_real_

    n_fin <- sum(is.finite(lam))
    lam_recent <- lam[1]
    trigger <- FALSE

    if (is.finite(lam_recent) && n_fin >= 2) {
      dmax <- max(abs(lam - lam_recent), na.rm = TRUE)
      sdev <- stats::sd(lam, na.rm = TRUE)

      n <- length(lam)
      q <- max(1L, min(compare_q, floor(n / 2L)))
      recent_idx <- 1:q
      older_idx  <- (n - q + 1L):n

      mean_recent <- mean(lam[recent_idx], na.rm = TRUE)
      mean_older  <- mean(lam[older_idx],  na.rm = TRUE)
      mean_gap <- abs(mean_recent - mean_older)

      x <- seq_len(n)
      ok_lm <- is.finite(lam)
      slope <- 0
      if (sum(ok_lm) >= 3) {
        slope <- as.numeric(stats::coef(stats::lm(lam[ok_lm] ~ x[ok_lm]))[2])
        if (!is.finite(slope)) slope <- 0
      }

      if (is.finite(dmax) && dmax > diff_thresh) trigger <- TRUE
      if (is.finite(sdev) && sdev > sd_thresh) trigger <- TRUE
      if (is.finite(mean_gap) && mean_gap > mean_gap_thresh) trigger <- TRUE
      if (is.finite(slope) && abs(slope) > slope_thresh) trigger <- TRUE
    }

    if (trigger) {
      min_keep <- if (!is.null(vg$min_keep)) as.integer(vg$min_keep) else 2L
      max_keep <- if (!is.null(vg$max_keep)) as.integer(vg$max_keep) else n
      max_keep <- max(1L, min(max_keep, n))

      prefix_rule <- if (!is.null(vg$prefix_rule)) as.character(vg$prefix_rule) else "sd"
      pref_sd_th  <- if (!is.null(vg$prefix_sd_thresh)) as.numeric(vg$prefix_sd_thresh) else sd_thresh
      pref_df_th  <- if (!is.null(vg$prefix_diff_thresh)) as.numeric(vg$prefix_diff_thresh) else diff_thresh

      lam_ok <- lam
      lam_ok[!is.finite(lam_ok)] <- NA_real_

      Kkeep <- 1L
      if (prefix_rule == "diff") {
        ref <- lam_ok[1]
        if (is.finite(ref)) {
          for (k in 2:max_keep) {
            if (!is.finite(lam_ok[k])) break
            if (abs(lam_ok[k] - ref) <= pref_df_th) Kkeep <- k else break
          }
        }
      } else {
        for (k in 2:max_keep) {
          if (!all(is.finite(lam_ok[1:k]))) break
          if (stats::sd(lam_ok[1:k]) <= pref_sd_th) Kkeep <- k else break
        }
      }

      Kkeep <- max(min_keep, Kkeep)
      ok <- ok & (seq_along(lambda_vec) <= Kkeep)
      if (!any(ok)) return(0)
    }
  }

  lambda_use <- lambda_vec[ok]
  sizes_use  <- sizes[ok]
  ends_use   <- ends[ok]

  if (weighting == "equal") {
    return(mean(lambda_use))
  }

  if (weighting == "size") {
    w <- sizes_use / sum(sizes_use)
    return(sum(w * lambda_use))
  }

  tau <- as.numeric(tau)
  if (!is.finite(tau) || tau <= 0) tau <- 100

  age <- pmax(0, past_end - ends_use)
  w_raw <- sizes_use * exp(-age / tau)

  if (!all(is.finite(w_raw)) || sum(w_raw) <= 0) {
    w <- sizes_use / sum(sizes_use)
    return(sum(w * lambda_use))
  }

  w <- w_raw / sum(w_raw)
  sum(w * lambda_use)
}

# Taylor approximation for the growth-optimal lambda based on e-values.
lambda_from_e_taylor <- function(e_vec, lambda_max = 1) {
  if (length(e_vec) == 0L) return(0)
  if (any(is.infinite(e_vec))) return(lambda_max)

  e_vec <- e_vec[is.finite(e_vec)]
  if (length(e_vec) == 0L) return(0)

  z <- e_vec - 1
  m1 <- mean(z)
  m2 <- mean(z^2)

  if (!is.finite(m1) || !is.finite(m2) || m2 <= 0) return(0)

  lam <- m1 / m2
  if (!is.finite(lam) || lam <= 0) return(0)

  min(lambda_max, lam)
}

# ------------------------------------------------------------
# Interval-wise lambda construction for VaR and ES
# GREE uses historical realised e-values.
# GREL re-evaluates past observations with the current forecast at time t.
# ------------------------------------------------------------

lambda_var_at_t <- function(t,
                            L_back, VaR_back, alpha_VaR,
                            e_var_back,
                            interval_mode = "preset",
                            preset = "S0",
                            weighting = "size_recency",
                            tau = 150,
                            lambda_max = 1,
                            var_guard = NULL) {

  intervals <- select_intervals(t = t, mode = interval_mode, preset = preset)
  if (!is.null(var_guard)) attr(intervals, "var_guard") <- var_guard

  if (attr(intervals, "past_end") <= 0L) {
    return(list(lambda_gree = 0, lambda_grel = 0, intervals = intervals))
  }

  lambda_vec_gree <- vapply(intervals, function(idx) {
    lambda_from_e_taylor(e_var_back[idx], lambda_max = lambda_max)
  }, numeric(1))

  lambda_gree <- combine_interval_lambdas(
    lambda_vec = lambda_vec_gree,
    interval_list = intervals,
    weighting = weighting,
    tau = tau
  )

  VaR_current <- VaR_back[t]

  lambda_vec_grel <- vapply(intervals, function(idx) {
    e_re <- (1 / (1 - alpha_VaR)) * as.numeric(L_back[idx] > VaR_current)
    lambda_from_e_taylor(e_re, lambda_max = lambda_max)
  }, numeric(1))

  lambda_grel <- combine_interval_lambdas(
    lambda_vec = lambda_vec_grel,
    interval_list = intervals,
    weighting = weighting,
    tau = tau
  )

  list(
    lambda_gree = lambda_gree,
    lambda_grel = lambda_grel,
    lambda_vec_gree = lambda_vec_gree,
    lambda_vec_grel = lambda_vec_grel,
    intervals = intervals
  )
}

lambda_es_at_t <- function(t,
                           L_back, ES_back, VaR_back, alpha_ES,
                           e_es_back,
                           interval_mode = "preset",
                           preset = "S0",
                           weighting = "size_recency",
                           tau = 150,
                           lambda_max = 1,
                           var_guard = NULL) {

  intervals <- select_intervals(t = t, mode = interval_mode, preset = preset)
  if (!is.null(var_guard)) attr(intervals, "var_guard") <- var_guard

  if (attr(intervals, "past_end") <= 0L) {
    return(list(lambda_gree = 0, lambda_grel = 0, intervals = intervals))
  }

  lambda_vec_gree <- vapply(intervals, function(idx) {
    lambda_from_e_taylor(e_es_back[idx], lambda_max = lambda_max)
  }, numeric(1))

  lambda_gree <- combine_interval_lambdas(
    lambda_vec = lambda_vec_gree,
    interval_list = intervals,
    weighting = weighting,
    tau = tau
  )

  ES_current  <- ES_back[t]
  VaR_current <- VaR_back[t]
  den <- (1 - alpha_ES) * (ES_current - VaR_current)

  if (!is.finite(den) || den <= 0) {
    lambda_grel <- lambda_max
    lambda_vec_grel <- rep(lambda_max, length(intervals))

    return(list(
      lambda_gree = lambda_gree,
      lambda_grel = lambda_grel,
      lambda_vec_gree = lambda_vec_gree,
      lambda_vec_grel = lambda_vec_grel,
      intervals = intervals
    ))
  }

  lambda_vec_grel <- vapply(intervals, function(idx) {
    e_re <- pmax(L_back[idx] - VaR_current, 0) / den
    lambda_from_e_taylor(e_re, lambda_max = lambda_max)
  }, numeric(1))

  lambda_grel <- combine_interval_lambdas(
    lambda_vec = lambda_vec_grel,
    interval_list = intervals,
    weighting = weighting,
    tau = tau
  )

  list(
    lambda_gree = lambda_gree,
    lambda_grel = lambda_grel,
    lambda_vec_gree = lambda_vec_gree,
    lambda_vec_grel = lambda_vec_grel,
    intervals = intervals
  )
}

# ============================================================
# 6) Strategy runners
# ============================================================

# ------------------------------------------------------------
# Single-method runners for GREE / GREL only
# ------------------------------------------------------------

run_var_e_process <- function(L_back, VaR_back, alpha_VaR,
                              interval_mode = "preset",
                              preset = "S0",
                              weighting = "size_recency",
                              tau = 150,
                              lambda_max = 0.99,
                              method = c("GREE", "GREL"),
                              var_guard = NULL) {
  method <- match.arg(method)

  T <- length(L_back)
  e_var_back <- e_values_var(L_back, VaR_back, alpha_VaR)

  M <- numeric(T + 1L)
  M[1] <- 0

  lambda_used <- numeric(T)

  for (t in 1:T) {
    lam_info <- lambda_var_at_t(
      t = t,
      L_back = L_back,
      VaR_back = VaR_back,
      alpha_VaR = alpha_VaR,
      e_var_back = e_var_back,
      interval_mode = interval_mode,
      preset = preset,
      weighting = weighting,
      tau = tau,
      lambda_max = lambda_max,
      var_guard = var_guard
    )

    lambda_t <- if (method == "GREE") lam_info$lambda_gree else lam_info$lambda_grel
    lambda_used[t] <- lambda_t

    M[t + 1L] <- update_e_process(M[t], e_var_back[t], lambda_t)
  }

  list(M = M, lambda = lambda_used, e = e_var_back)
}

run_es_e_process <- function(L_back, ES_back, VaR_back, alpha_ES,
                             interval_mode = "preset",
                             preset = "S0",
                             weighting = "size_recency",
                             tau = 150,
                             lambda_max = 0.99,
                             method = c("GREE", "GREL"),
                             var_guard = NULL) {
  method <- match.arg(method)

  T <- length(L_back)
  stopifnot(length(ES_back) == T, length(VaR_back) == T)

  e_es_back <- e_values_es(L_back, ES_back, VaR_back, alpha_ES)

  M <- numeric(T + 1L)
  M[1] <- 0

  lambda_used <- numeric(T)

  for (t in 1:T) {
    lam_info <- lambda_es_at_t(
      t = t,
      L_back = L_back,
      ES_back = ES_back,
      VaR_back = VaR_back,
      alpha_ES = alpha_ES,
      e_es_back = e_es_back,
      interval_mode = interval_mode,
      preset = preset,
      weighting = weighting,
      tau = tau,
      lambda_max = lambda_max,
      var_guard = var_guard
    )

    lambda_t <- if (method == "GREE") lam_info$lambda_gree else lam_info$lambda_grel
    lambda_used[t] <- lambda_t

    M[t + 1L] <- update_e_process(M[t], e_es_back[t], lambda_t)
  }

  list(M = M, lambda = lambda_used, e = e_es_back)
}

run_var_both_methods <- function(L_back, VaR_back, alpha_VaR,
                                 interval_mode = "preset",
                                 preset = "S0",
                                 weighting = "size_recency",
                                 tau = 150,
                                 lambda_max = 1,
                                 var_guard = NULL) {

  out_gree <- run_var_e_process(
    L_back = L_back,
    VaR_back = VaR_back,
    alpha_VaR = alpha_VaR,
    interval_mode = interval_mode,
    preset = preset,
    weighting = weighting,
    tau = tau,
    lambda_max = lambda_max,
    method = "GREE",
    var_guard = var_guard
  )

  out_grel <- run_var_e_process(
    L_back = L_back,
    VaR_back = VaR_back,
    alpha_VaR = alpha_VaR,
    interval_mode = interval_mode,
    preset = preset,
    weighting = weighting,
    tau = tau,
    lambda_max = lambda_max,
    method = "GREL",
    var_guard = var_guard
  )

  list(
    GREE = out_gree,
    GREL = out_grel,
    settings = list(
      interval_mode = interval_mode,
      preset = preset,
      weighting = weighting,
      tau = tau,
      lambda_max = lambda_max,
      var_guard = var_guard
    )
  )
}

run_es_both_methods <- function(L_back, ES_back, VaR_back, alpha_ES,
                                interval_mode = "preset",
                                preset = "S0",
                                weighting = "size_recency",
                                tau = 150,
                                lambda_max = 1,
                                var_guard = NULL) {

  out_gree <- run_es_e_process(
    L_back = L_back,
    ES_back = ES_back,
    VaR_back = VaR_back,
    alpha_ES = alpha_ES,
    interval_mode = interval_mode,
    preset = preset,
    weighting = weighting,
    tau = tau,
    lambda_max = lambda_max,
    method = "GREE",
    var_guard = var_guard
  )

  out_grel <- run_es_e_process(
    L_back = L_back,
    ES_back = ES_back,
    VaR_back = VaR_back,
    alpha_ES = alpha_ES,
    interval_mode = interval_mode,
    preset = preset,
    weighting = weighting,
    tau = tau,
    lambda_max = lambda_max,
    method = "GREL",
    var_guard = var_guard
  )

  list(
    GREE = out_gree,
    GREL = out_grel,
    settings = list(
      interval_mode = interval_mode,
      preset = preset,
      weighting = weighting,
      tau = tau,
      lambda_max = lambda_max,
      var_guard = var_guard
    )
  )
}

# ------------------------------------------------------------
# Diagnostics for COV
# The covariance-based mixture scores each interval by how "stable" the loss
# series and e-value series look locally.
# ------------------------------------------------------------

cap_series <- function(x, cap = 50) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- cap
  pmin(x, cap)
}

autocorr_score <- function(x, K = 5L, cap = Inf) {
  x <- as.numeric(x)
  if (is.finite(cap)) x <- cap_series(x, cap = cap)
  x <- x[is.finite(x)]

  W <- length(x)
  if (W <= 3) return(0)

  K_eff <- min(as.integer(K), W - 2L)
  if (K_eff <= 0L) return(0)

  ac <- stats::acf(x, lag.max = K_eff, plot = FALSE)$acf[-1]
  sum(as.numeric(ac)^2)
}

# Stable softmax mapping from diagnostic scores to a mixing weight.
alpha_from_cov_scores <- function(S_e, S_L, S_F = NULL,
                                  kappa = 5,
                                  kappa_F = 2) {
  q_e <- 1 / (1 + max(S_e, 0))
  q_L <- 1 / (1 + max(S_L, 0))

  if (is.null(S_F) || !is.finite(S_F)) {
    noise_term <- 0
  } else {
    q_F <- 1 / (1 + max(S_F, 0))
    noise_term <- (1 - q_F)
  }

  a <- kappa * q_e + kappa_F * noise_term
  b <- kappa * q_L

  m <- max(a, b)
  ea <- exp(a - m)
  eb <- exp(b - m)
  ea / (ea + eb)
}

alpha_vec_cov_at_t <- function(intervals, L_back, e_series,
                               forecast_series = NULL,
                               W_diag = 50L,
                               K = 5L,
                               kappa = 5,
                               kappa_F = 2,
                               cap_e = 50) {
  vapply(intervals, function(idx) {
    if (length(idx) <= 2L) return(0.5)

    idx2 <- tail(idx, min(length(idx), as.integer(W_diag)))

    S_L <- autocorr_score(L_back[idx2], K = K, cap = Inf)
    S_e <- autocorr_score(e_series[idx2], K = K, cap = cap_e)

    S_F <- NULL
    if (!is.null(forecast_series)) {
      f <- as.numeric(forecast_series[idx2])
      f <- f[is.finite(f)]
      if (length(f) >= 4) {
        df <- diff(f)
        if (length(df) >= 2 && any(is.finite(df))) {
          S_F <- stats::var(df, na.rm = TRUE)
          if (!is.finite(S_F)) S_F <- NULL
        }
      }
    }

    alpha_from_cov_scores(
      S_e = S_e,
      S_L = S_L,
      S_F = S_F,
      kappa = kappa,
      kappa_F = kappa_F
    )
  }, numeric(1))
}

# ------------------------------------------------------------
# GREM
# Wealth-proportional mixture of the GREE and GREL candidate lambdas.
# ------------------------------------------------------------

run_var_grem <- function(L_back, VaR_back, alpha_VaR,
                         interval_mode = "preset",
                         preset = "S0",
                         weighting = "size_recency",
                         tau = 150,
                         lambda_max = 1,
                         var_guard = NULL) {

  T <- length(L_back)
  e_var_back <- e_values_var(L_back, VaR_back, alpha_VaR)

  M_E   <- numeric(T + 1L); M_E[1]   <- 0
  M_L   <- numeric(T + 1L); M_L[1]   <- 0
  M_mix <- numeric(T + 1L); M_mix[1] <- 0

  lam_E <- numeric(T)
  lam_L <- numeric(T)
  lam_M <- numeric(T)

  for (t in 1:T) {
    info <- lambda_var_at_t(
      t = t,
      L_back = L_back,
      VaR_back = VaR_back,
      alpha_VaR = alpha_VaR,
      e_var_back = e_var_back,
      interval_mode = interval_mode,
      preset = preset,
      weighting = weighting,
      tau = tau,
      lambda_max = lambda_max,
      var_guard = var_guard
    )

    lam_E[t] <- info$lambda_gree
    lam_L[t] <- info$lambda_grel

    M_E[t + 1L] <- update_e_process(M_E[t], e_var_back[t], lam_E[t])
    M_L[t + 1L] <- update_e_process(M_L[t], e_var_back[t], lam_L[t])

    wE <- weight_from_logs(M_E[t], M_L[t])
    lam_M[t] <- wE * lam_E[t] + (1 - wE) * lam_L[t]

    M_mix[t + 1L] <- update_e_process(M_mix[t], e_var_back[t], lam_M[t])
  }

  list(
    M = M_mix,
    lambda = lam_M,
    base = list(M_E = M_E, M_L = M_L, lambda_E = lam_E, lambda_L = lam_L),
    e = e_var_back
  )
}

run_es_grem <- function(L_back, ES_back, VaR_back, alpha_ES,
                        interval_mode = "preset",
                        preset = "S0",
                        weighting = "size_recency",
                        tau = 150,
                        lambda_max = 1,
                        var_guard = NULL) {

  T <- length(L_back)
  e_es_back <- e_values_es(L_back, ES_back, VaR_back, alpha_ES)

  M_E   <- numeric(T + 1L); M_E[1]   <- 0
  M_L   <- numeric(T + 1L); M_L[1]   <- 0
  M_mix <- numeric(T + 1L); M_mix[1] <- 0

  lam_E <- numeric(T)
  lam_L <- numeric(T)
  lam_M <- numeric(T)

  for (t in 1:T) {
    info <- lambda_es_at_t(
      t = t,
      L_back = L_back,
      ES_back = ES_back,
      VaR_back = VaR_back,
      alpha_ES = alpha_ES,
      e_es_back = e_es_back,
      interval_mode = interval_mode,
      preset = preset,
      weighting = weighting,
      tau = tau,
      lambda_max = lambda_max,
      var_guard = var_guard
    )

    lam_E[t] <- info$lambda_gree
    lam_L[t] <- info$lambda_grel

    M_E[t + 1L] <- update_e_process(M_E[t], e_es_back[t], lam_E[t])
    M_L[t + 1L] <- update_e_process(M_L[t], e_es_back[t], lam_L[t])

    wE <- weight_from_logs(M_E[t], M_L[t])
    lam_M[t] <- wE * lam_E[t] + (1 - wE) * lam_L[t]

    M_mix[t + 1L] <- update_e_process(M_mix[t], e_es_back[t], lam_M[t])
  }

  list(
    M = M_mix,
    lambda = lam_M,
    base = list(M_E = M_E, M_L = M_L, lambda_E = lam_E, lambda_L = lam_L),
    e = e_es_back
  )
}

# ------------------------------------------------------------
# WEALTH
# Similar to GREM, but uses a direct wealth-adaptive mixing parameter alpha_t.
# ------------------------------------------------------------

run_var_wealth_mix <- function(L_back, VaR_back, alpha_VaR,
                               kappa = 1,
                               interval_mode = "preset",
                               preset = "S0",
                               weighting = "size_recency",
                               tau = 150,
                               lambda_max = 1,
                               var_guard = NULL) {

  T <- length(L_back)
  e_var_back <- e_values_var(L_back, VaR_back, alpha_VaR)

  M_E   <- numeric(T + 1L); M_E[1]   <- 0
  M_L   <- numeric(T + 1L); M_L[1]   <- 0
  M_mix <- numeric(T + 1L); M_mix[1] <- 0

  lam_E <- numeric(T)
  lam_L <- numeric(T)
  lam_M <- numeric(T)
  alpha_t <- numeric(T)

  for (t in 1:T) {
    info <- lambda_var_at_t(
      t = t,
      L_back = L_back,
      VaR_back = VaR_back,
      alpha_VaR = alpha_VaR,
      e_var_back = e_var_back,
      interval_mode = interval_mode,
      preset = preset,
      weighting = weighting,
      tau = tau,
      lambda_max = lambda_max,
      var_guard = var_guard
    )

    lam_E[t] <- info$lambda_gree
    lam_L[t] <- info$lambda_grel

    M_E[t + 1L] <- update_e_process(M_E[t], e_var_back[t], lam_E[t])
    M_L[t + 1L] <- update_e_process(M_L[t], e_var_back[t], lam_L[t])

    alpha_t[t] <- weight_from_logs(kappa * M_E[t], kappa * M_L[t])

    lam_M[t] <- alpha_t[t] * lam_E[t] + (1 - alpha_t[t]) * lam_L[t]
    M_mix[t + 1L] <- update_e_process(M_mix[t], e_var_back[t], lam_M[t])
  }

  list(
    M = M_mix,
    lambda = lam_M,
    alpha = alpha_t,
    base = list(M_E = M_E, M_L = M_L, lambda_E = lam_E, lambda_L = lam_L),
    e = e_var_back
  )
}

run_es_wealth_mix <- function(L_back, ES_back, VaR_back, alpha_ES,
                              kappa = 1,
                              interval_mode = "preset",
                              preset = "S0",
                              weighting = "size_recency",
                              tau = 150,
                              lambda_max = 1,
                              var_guard = NULL) {

  T <- length(L_back)
  e_es_back <- e_values_es(L_back, ES_back, VaR_back, alpha_ES)

  M_E   <- numeric(T + 1L); M_E[1]   <- 0
  M_L   <- numeric(T + 1L); M_L[1]   <- 0
  M_mix <- numeric(T + 1L); M_mix[1] <- 0

  lam_E <- numeric(T)
  lam_L <- numeric(T)
  lam_M <- numeric(T)
  alpha_t <- numeric(T)

  for (t in 1:T) {
    info <- lambda_es_at_t(
      t = t,
      L_back = L_back,
      ES_back = ES_back,
      VaR_back = VaR_back,
      alpha_ES = alpha_ES,
      e_es_back = e_es_back,
      interval_mode = interval_mode,
      preset = preset,
      weighting = weighting,
      tau = tau,
      lambda_max = lambda_max,
      var_guard = var_guard
    )

    lam_E[t] <- info$lambda_gree
    lam_L[t] <- info$lambda_grel

    M_E[t + 1L] <- update_e_process(M_E[t], e_es_back[t], lam_E[t])
    M_L[t + 1L] <- update_e_process(M_L[t], e_es_back[t], lam_L[t])

    alpha_t[t] <- weight_from_logs(kappa * M_E[t], kappa * M_L[t])

    lam_M[t] <- alpha_t[t] * lam_E[t] + (1 - alpha_t[t]) * lam_L[t]
    M_mix[t + 1L] <- update_e_process(M_mix[t], e_es_back[t], lam_M[t])
  }

  list(
    M = M_mix,
    lambda = lam_M,
    alpha = alpha_t,
    base = list(M_E = M_E, M_L = M_L, lambda_E = lam_E, lambda_L = lam_L),
    e = e_es_back
  )
}

# ------------------------------------------------------------
# COV
# Covariance-adaptive interval-wise mixture using diagnostic scores.
# ------------------------------------------------------------

run_var_cov_mix <- function(L_back, VaR_back, alpha_VaR,
                            K = 5L,
                            kappa_cov = 5,
                            W_diag = 50L,
                            K_diag = 5L,
                            kappa_F = -2,
                            cap_e = 50,
                            interval_mode = "preset",
                            preset = "S0",
                            weighting = "size_recency",
                            tau = 150,
                            lambda_max = 1,
                            var_guard = NULL) {

  T <- length(L_back)
  e_var_back <- e_values_var(L_back, VaR_back, alpha_VaR)

  M_mix <- numeric(T + 1L); M_mix[1] <- 0
  lam_M <- numeric(T)
  alpha_used <- vector("list", T)

  for (t in 1:T) {
    info <- lambda_var_at_t(
      t = t,
      L_back = L_back,
      VaR_back = VaR_back,
      alpha_VaR = alpha_VaR,
      e_var_back = e_var_back,
      interval_mode = interval_mode,
      preset = preset,
      weighting = weighting,
      tau = tau,
      lambda_max = lambda_max,
      var_guard = var_guard
    )

    intervals <- info$intervals
    if (attr(intervals, "past_end") <= 0L) {
      lam_M[t] <- 0
      M_mix[t + 1L] <- update_e_process(M_mix[t], e_var_back[t], lam_M[t])
      next
    }

    alpha_vec <- alpha_vec_cov_at_t(
      intervals = intervals,
      L_back = L_back,
      e_series = e_var_back,
      forecast_series = VaR_back,
      W_diag = W_diag,
      K = K_diag,
      kappa = kappa_cov,
      kappa_F = kappa_F,
      cap_e = cap_e
    )

    alpha_used[[t]] <- alpha_vec

    lam_vec_mix <- alpha_vec * info$lambda_vec_gree + (1 - alpha_vec) * info$lambda_vec_grel

    lam_M[t] <- combine_interval_lambdas(
      lambda_vec = lam_vec_mix,
      interval_list = intervals,
      weighting = weighting,
      tau = tau
    )

    M_mix[t + 1L] <- update_e_process(M_mix[t], e_var_back[t], lam_M[t])
  }

  list(M = M_mix, lambda = lam_M, alpha_vec = alpha_used, e = e_var_back)
}

run_es_cov_mix <- function(L_back, ES_back, VaR_back, alpha_ES,
                           K = 5L,
                           kappa_cov = 5,
                           W_diag = 50L,
                           K_diag = 5L,
                           kappa_F = 2,
                           cap_e = 100,
                           interval_mode = "preset",
                           preset = "S0",
                           weighting = "size_recency",
                           tau = 150,
                           lambda_max = 1,
                           var_guard = NULL) {

  T <- length(L_back)
  e_es_back <- e_values_es(L_back, ES_back, VaR_back, alpha_ES)

  M_mix <- numeric(T + 1L); M_mix[1] <- 0
  lam_M <- numeric(T)
  alpha_used <- vector("list", T)

  for (t in 1:T) {
    info <- lambda_es_at_t(
      t = t,
      L_back = L_back,
      ES_back = ES_back,
      VaR_back = VaR_back,
      alpha_ES = alpha_ES,
      e_es_back = e_es_back,
      interval_mode = interval_mode,
      preset = preset,
      weighting = weighting,
      tau = tau,
      lambda_max = lambda_max,
      var_guard = var_guard
    )

    intervals <- info$intervals
    if (attr(intervals, "past_end") <= 0L) {
      lam_M[t] <- 0
      M_mix[t + 1L] <- update_e_process(M_mix[t], e_es_back[t], lam_M[t])
      next
    }

    alpha_vec <- alpha_vec_cov_at_t(
      intervals = intervals,
      L_back = L_back,
      e_series = e_es_back,
      forecast_series = (ES_back - VaR_back),
      W_diag = W_diag,
      K = K_diag,
      kappa = kappa_cov,
      kappa_F = kappa_F,
      cap_e = cap_e
    )

    alpha_used[[t]] <- alpha_vec

    lam_vec_mix <- alpha_vec * info$lambda_vec_gree + (1 - alpha_vec) * info$lambda_vec_grel

    lam_M[t] <- combine_interval_lambdas(
      lambda_vec = lam_vec_mix,
      interval_list = intervals,
      weighting = weighting,
      tau = tau
    )

    M_mix[t + 1L] <- update_e_process(M_mix[t], e_es_back[t], lam_M[t])
  }

  list(M = M_mix, lambda = lam_M, alpha_vec = alpha_used, e = e_es_back)
}

# ============================================================
# 7) One replication
# ============================================================

run_one_rep <- function(world = c("A", "B", "A7", "B7", "C7", "D"),
                        risk  = c("VaR", "ES"),
                        b_star = NULL,
                        f = 1.0,
                        step_f_low = 0.90,
                        step_f_high = 1.10,
                        step_center = NULL,
                        step_width  = NULL,
                        step_apply_to = "both",
                        distort_mode_var = c("none", "scale", "smooth_step"),
                        distort_mode_es  = c("none", "both", "es_only", "smooth_step"),
                        interval_mode = "preset",
                        preset = "S0",
                        weighting = "size_recency",
                        tau = 150,
                        lambda_max = 0.99,
                        wealth_kappa = 4,
                        cov_K = 5L,
                        cov_kappa = 5,
                        cov_W_diag  = 200L,
                        cov_K_diag  = 10L,
                        cov_kappa_F = 2,
                        cov_cap_e_var = 50,
                        cov_cap_e_es  = 100,
                        q_type = 8,
                        forecast_mode = c("direct_world", "oracle_empirical", "oracle_true", "ml_norm", "ml_t"),
                        refit_every = 50L,
                        oracle_mc_N = 200000L,
                        segments = NULL,
                        ex7_under_factor = 0.90,
                        ex7_theta = 0.01,
                        ex7_noise_set = c(0, as.vector(sapply(1:5, function(i) c(-i, i))) / 10),
                        var_guard = NULL) {

  world <- match.arg(world)
  risk  <- match.arg(risk)
  distort_mode_var <- match.arg(distort_mode_var)
  distort_mode_es  <- match.arg(distort_mode_es)
  forecast_mode <- match.arg(forecast_mode)
  forecast_mode <- as.character(forecast_mode)[1]

  # Example-7 worlds use alpha = 0.95 in the direct setup.
  if (world %in% c("A7", "B7", "C7", "D")) {
    alpha_VaR <- 0.95
    alpha_ES  <- 0.95
  }

  sim <- switch(
    world,
    "A"  = simulate_world_A(T_total),
    "B"  = simulate_world_B(T_total, b_star = b_star),
    "A7" = simulate_world_A7(T_total, alpha = alpha_VaR),
    "B7" = simulate_world_B7(T_total, alpha = alpha_VaR, theta = ex7_theta),
    "C7" = simulate_world_C7(T_total, alpha = alpha_VaR, eps_set = ex7_noise_set),
    "D"  = simulate_world_D(
      T_total,
      segments = segments,
      alpha = alpha_VaR,
      theta = ex7_theta,
      eps_set = ex7_noise_set
    )
  )

  sl <- get_backtest_slice(sim, T_train = T_train, T_backtest = T_backtest)
  L_back <- sl$L_back

  out <- list()

  direct_fc_available <- !is.null(sl$VaR_back) && (risk == "VaR" || !is.null(sl$ES_back))
  direct_fc <- direct_fc_available && (forecast_mode == "direct_world")

  if (!direct_fc) {
    if (world == "A") {
      fc <- build_forecasts_worldA(
        sim = sim,
        forecast_mode = forecast_mode,
        T_train = T_train,
        T_backtest = T_backtest,
        alpha_VaR = alpha_VaR,
        alpha_ES = alpha_ES,
        q_type = q_type,
        refit_every = refit_every,
        oracle_mc_N = oracle_mc_N
      )

    } else if (world == "B") {
      fcB <- make_ec44_forecasts_worldB(
        L_full = sim$L,
        alpha_VaR = alpha_VaR,
        alpha_ES  = alpha_ES,
        presample_length = T_train,
        backtest_length  = T_backtest,
        q_type = q_type
      )
      fc <- list(VaR = fcB$VaR, ES = fcB$ES, mode = "ec44_qml")

    } else if (world == "D") {
      if (forecast_mode == "oracle_true") {
        cV <- const_norm(alpha_VaR)
        cE <- const_norm(alpha_ES)

        idx_bt <- (T_train + 1L):(T_train + T_backtest)
        VaR <- sim$mu[idx_bt] + sim$sig[idx_bt] * cV$q
        ES  <- sim$mu[idx_bt] + sim$sig[idx_bt] * cE$es

        fc <- list(VaR = VaR, ES = ES, mode = "oracle_true")

      } else if (forecast_mode == "oracle_empirical") {
        fc0 <- make_empirical_baseline_forecasts(
          sim = sim,
          alpha_VaR = alpha_VaR,
          alpha_ES = alpha_ES,
          T_train = T_train,
          T_backtest = T_backtest,
          q_type = q_type
        )
        fc <- list(VaR = fc0$VaR, ES = fc0$ES, mode = "oracle_empirical")

      } else if (forecast_mode == "direct_world") {
        fc <- list(VaR = sl$VaR_back, ES = sl$ES_back, mode = "direct_world")

      } else {
        fc0 <- make_empirical_baseline_forecasts(
          sim = sim,
          alpha_VaR = alpha_VaR,
          alpha_ES = alpha_ES,
          T_train = T_train,
          T_backtest = T_backtest,
          q_type = q_type
        )
        fc <- list(VaR = fc0$VaR, ES = fc0$ES, mode = "oracle_empirical")
      }

    } else {
      fc <- make_empirical_baseline_forecasts(
        sim = sim,
        alpha_VaR = alpha_VaR,
        alpha_ES = alpha_ES,
        T_train = T_train,
        T_backtest = T_backtest,
        q_type = q_type
      )
      fc$mode <- "oracle_empirical"
    }

    if (risk == "VaR") {
      VaR0 <- fc$VaR

      d <- distort_var_forecast(
        VaR_base = VaR0,
        f = f,
        mode = distort_mode_var,
        f_low = step_f_low,
        f_high = step_f_high,
        center = step_center,
        width  = step_width
      )
      VaR <- d$VaR

      out$dist_mult  <- d$mult
      out$dist_label <- d$label

      both <- run_var_both_methods(
        L_back = L_back,
        VaR_back = VaR,
        alpha_VaR = alpha_VaR,
        interval_mode = interval_mode,
        preset = preset,
        weighting = weighting,
        tau = tau,
        lambda_max = lambda_max,
        var_guard = var_guard
      )

      out$GREE <- both$GREE
      out$GREL <- both$GREL
      out$GREM <- run_var_grem(
        L_back = L_back,
        VaR_back = VaR,
        alpha_VaR = alpha_VaR,
        interval_mode = interval_mode,
        preset = preset,
        weighting = weighting,
        tau = tau,
        lambda_max = lambda_max,
        var_guard = var_guard
      )
      out$WEALTH <- run_var_wealth_mix(
        L_back = L_back,
        VaR_back = VaR,
        alpha_VaR = alpha_VaR,
        kappa = wealth_kappa,
        interval_mode = interval_mode,
        preset = preset,
        weighting = weighting,
        tau = tau,
        lambda_max = lambda_max,
        var_guard = var_guard
      )
      out$COV <- run_var_cov_mix(
        L_back = L_back,
        VaR_back = VaR,
        alpha_VaR = alpha_VaR,
        K = cov_K,
        kappa_cov = cov_kappa,
        W_diag = cov_W_diag,
        K_diag = cov_K_diag,
        kappa_F = cov_kappa_F,
        cap_e = cov_cap_e_var,
        interval_mode = interval_mode,
        preset = preset,
        weighting = weighting,
        tau = tau,
        lambda_max = lambda_max,
        var_guard = var_guard
      )

      out$series <- list(L = L_back, VaR = VaR, VaR_base = VaR0)

    } else {
      ES0  <- fc$ES
      VaR0 <- fc$VaR

      d <- distort_es_forecasts(
        ES_base = ES0,
        VaR_base = VaR0,
        f = f,
        mode = distort_mode_es,
        f_low = step_f_low,
        f_high = step_f_high,
        center = step_center,
        width  = step_width,
        apply_to = step_apply_to
      )
      ES  <- d$ES
      VaR <- d$VaR

      out$dist_mult  <- d$mult
      out$dist_label <- d$label

      both <- run_es_both_methods(
        L_back = L_back,
        ES_back = ES,
        VaR_back = VaR,
        alpha_ES = alpha_ES,
        interval_mode = interval_mode,
        preset = preset,
        weighting = weighting,
        tau = tau,
        lambda_max = lambda_max,
        var_guard = var_guard
      )

      out$GREE <- both$GREE
      out$GREL <- both$GREL
      out$GREM <- run_es_grem(
        L_back = L_back,
        ES_back = ES,
        VaR_back = VaR,
        alpha_ES = alpha_ES,
        interval_mode = interval_mode,
        preset = preset,
        weighting = weighting,
        tau = tau,
        lambda_max = lambda_max,
        var_guard = var_guard
      )
      out$WEALTH <- run_es_wealth_mix(
        L_back = L_back,
        ES_back = ES,
        VaR_back = VaR,
        alpha_ES = alpha_ES,
        kappa = wealth_kappa,
        interval_mode = interval_mode,
        preset = preset,
        weighting = weighting,
        tau = tau,
        lambda_max = lambda_max,
        var_guard = var_guard
      )
      out$COV <- run_es_cov_mix(
        L_back = L_back,
        ES_back = ES,
        VaR_back = VaR,
        alpha_ES = alpha_ES,
        K = cov_K,
        kappa_cov = cov_kappa,
        W_diag = cov_W_diag,
        K_diag = cov_K_diag,
        kappa_F = cov_kappa_F,
        cap_e = cov_cap_e_es,
        interval_mode = interval_mode,
        preset = preset,
        weighting = weighting,
        tau = tau,
        lambda_max = lambda_max,
        var_guard = var_guard
      )

      out$series <- list(L = L_back, ES = ES, VaR = VaR, ES_base = ES0, VaR_base = VaR0)
    }

  } else {
    if (risk == "VaR") {
      VaR0 <- sl$VaR_back

      d <- distort_var_forecast(
        VaR_base = VaR0,
        f = f,
        mode = distort_mode_var,
        f_low = step_f_low,
        f_high = step_f_high,
        center = step_center,
        width  = step_width
      )
      VaR <- d$VaR

      out$dist_mult  <- d$mult
      out$dist_label <- d$label

      both <- run_var_both_methods(
        L_back = L_back,
        VaR_back = VaR,
        alpha_VaR = alpha_VaR,
        interval_mode = interval_mode,
        preset = preset,
        weighting = weighting,
        tau = tau,
        lambda_max = lambda_max,
        var_guard = var_guard
      )

      out$GREE <- both$GREE
      out$GREL <- both$GREL
      out$GREM <- run_var_grem(
        L_back = L_back,
        VaR_back = VaR,
        alpha_VaR = alpha_VaR,
        interval_mode = interval_mode,
        preset = preset,
        weighting = weighting,
        tau = tau,
        lambda_max = lambda_max,
        var_guard = var_guard
      )
      out$WEALTH <- run_var_wealth_mix(
        L_back = L_back,
        VaR_back = VaR,
        alpha_VaR = alpha_VaR,
        kappa = wealth_kappa,
        interval_mode = interval_mode,
        preset = preset,
        weighting = weighting,
        tau = tau,
        lambda_max = lambda_max,
        var_guard = var_guard
      )
      out$COV <- run_var_cov_mix(
        L_back = L_back,
        VaR_back = VaR,
        alpha_VaR = alpha_VaR,
        K = cov_K,
        kappa_cov = cov_kappa,
        W_diag = cov_W_diag,
        K_diag = cov_K_diag,
        kappa_F = cov_kappa_F,
        cap_e = cov_cap_e_var,
        interval_mode = interval_mode,
        preset = preset,
        weighting = weighting,
        tau = tau,
        lambda_max = lambda_max,
        var_guard = var_guard
      )

      out$series <- list(L = L_back, VaR = VaR, VaR_base = VaR0)

    } else {
      VaR0 <- sl$VaR_back
      ES0  <- sl$ES_back

      d <- distort_es_forecasts(
        ES_base = ES0,
        VaR_base = VaR0,
        f = f,
        mode = distort_mode_es,
        f_low = step_f_low,
        f_high = step_f_high,
        center = step_center,
        width  = step_width,
        apply_to = step_apply_to
      )
      ES  <- d$ES
      VaR <- d$VaR

      out$dist_mult  <- d$mult
      out$dist_label <- d$label

      both <- run_es_both_methods(
        L_back = L_back,
        ES_back = ES,
        VaR_back = VaR,
        alpha_ES = alpha_ES,
        interval_mode = interval_mode,
        preset = preset,
        weighting = weighting,
        tau = tau,
        lambda_max = lambda_max,
        var_guard = var_guard
      )

      out$GREE <- both$GREE
      out$GREL <- both$GREL
      out$GREM <- run_es_grem(
        L_back = L_back,
        ES_back = ES,
        VaR_back = VaR,
        alpha_ES = alpha_ES,
        interval_mode = interval_mode,
        preset = preset,
        weighting = weighting,
        tau = tau,
        lambda_max = lambda_max,
        var_guard = var_guard
      )
      out$WEALTH <- run_es_wealth_mix(
        L_back = L_back,
        ES_back = ES,
        VaR_back = VaR,
        alpha_ES = alpha_ES,
        kappa = wealth_kappa,
        interval_mode = interval_mode,
        preset = preset,
        weighting = weighting,
        tau = tau,
        lambda_max = lambda_max,
        var_guard = var_guard
      )
      out$COV <- run_es_cov_mix(
        L_back = L_back,
        ES_back = ES,
        VaR_back = VaR,
        alpha_ES = alpha_ES,
        K = cov_K,
        kappa_cov = cov_kappa,
        W_diag = cov_W_diag,
        K_diag = cov_K_diag,
        kappa_F = cov_kappa_F,
        cap_e = cov_cap_e_es,
        interval_mode = interval_mode,
        preset = preset,
        weighting = weighting,
        tau = tau,
        lambda_max = lambda_max,
        var_guard = var_guard
      )

      out$series <- list(L = L_back, ES = ES, VaR = VaR, ES_base = ES0, VaR_base = VaR0)
    }
  }

  out$meta <- sim$meta
  if (!is.null(sim$regime)) out$regime_full <- sim$regime
  out$forecast_mode <- if (!direct_fc) fc$mode else "direct_world"
  out$refit_every   <- if (!direct_fc && !is.null(fc$refit_every)) fc$refit_every else NA_integer_

  out
}

# ============================================================
# 8) Repeated runs and detection summaries
# ============================================================

detect_time <- function(M, c = 5) {
  thr <- log(c)
  hit <- which(M >= thr)
  if (length(hit) == 0) return(Inf)
  hit[1] - 1L
}

run_many <- function(R = 1000,
                     seed0 = 1,
                     world = "A",
                     risk  = "ES",
                     b_star = NULL,
                     f = 0.9,
                     step_f_low = 0.90,
                     step_f_high = 1.10,
                     step_center = NULL,
                     step_width  = NULL,
                     step_apply_to = "both",
                     distort_mode_var = "scale",
                     distort_mode_es  = "both",
                     interval_mode = "preset",
                     preset = "S0",
                     weighting = "size_recency",
                     tau = 100,
                     lambda_max = 0.99,
                     wealth_kappa = 4,
                     cov_K = 5L,
                     cov_kappa = 5,
                     cov_W_diag  = 200L,
                     cov_K_diag  = 10L,
                     cov_kappa_F = 2,
                     cov_cap_e_var = 50,
                     cov_cap_e_es  = 100,
                     forecast_mode = c("direct_world", "oracle_empirical", "oracle_true", "ml_norm", "ml_t"),
                     refit_every = 50L,
                     oracle_mc_N = 200000L,
                     segments = NULL,
                     ex7_under_factor = 0.90,
                     ex7_theta = 0.01,
                     ex7_noise_set = c(0, as.vector(sapply(1:5, function(i) c(-i, i))) / 10),
                     var_guard = NULL) {
  forecast_mode <- as.character(forecast_mode)[1]

  methods <- c("GREE", "GREL", "GREM", "WEALTH", "COV")
  T <- T_backtest

  sum_VaR <- rep(0, T)
  sum_ES  <- rep(0, T)
  n_fc    <- 0L

  M_store <- lapply(methods, function(.) matrix(NA_real_, nrow = R, ncol = T + 1L))
  names(M_store) <- methods

  det_store <- lapply(methods, function(.) {
    matrix(
      NA_real_,
      nrow = R,
      ncol = 3,
      dimnames = list(NULL, c("c2", "c5", "c10"))
    )
  })
  names(det_store) <- methods

  for (r in 1:R) {
    set.seed(seed0 + r)

    rep_out <- run_one_rep(
      world = world,
      risk = risk,
      b_star = b_star,
      f = f,
      step_f_low = step_f_low,
      step_f_high = step_f_high,
      step_center = step_center,
      step_width  = step_width,
      step_apply_to = step_apply_to,
      distort_mode_var = distort_mode_var,
      distort_mode_es  = distort_mode_es,
      interval_mode = interval_mode,
      preset = preset,
      weighting = weighting,
      tau = tau,
      lambda_max = lambda_max,
      wealth_kappa = wealth_kappa,
      cov_K = cov_K,
      cov_kappa = cov_kappa,
      cov_W_diag  = cov_W_diag,
      cov_K_diag  = cov_K_diag,
      cov_kappa_F = cov_kappa_F,
      cov_cap_e_var = cov_cap_e_var,
      cov_cap_e_es  = cov_cap_e_es,
      forecast_mode = forecast_mode,
      refit_every = refit_every,
      oracle_mc_N = oracle_mc_N,
      segments = segments,
      ex7_under_factor = ex7_under_factor,
      ex7_theta = ex7_theta,
      ex7_noise_set = ex7_noise_set,
      var_guard = var_guard
    )

    VaR_r <- rep_out$series$VaR
    ES_r  <- if (risk == "ES") rep_out$series$ES else NULL

    ok_var <- is.numeric(VaR_r) && length(VaR_r) == T && all(is.finite(VaR_r))
    ok_es  <- (risk != "ES") || (is.numeric(ES_r) && length(ES_r) == T && all(is.finite(ES_r)))

    if (ok_var && ok_es) {
      sum_VaR <- sum_VaR + VaR_r
      if (risk == "ES") sum_ES <- sum_ES + ES_r
      n_fc <- n_fc + 1L
    }

    for (m in methods) {
      M_store[[m]][r, ] <- rep_out[[m]]$M
      det_store[[m]][r, "c2"]  <- detect_time(rep_out[[m]]$M, c = 2)
      det_store[[m]][r, "c5"]  <- detect_time(rep_out[[m]]$M, c = 5)
      det_store[[m]][r, "c10"] <- detect_time(rep_out[[m]]$M, c = 10)
    }
  }

  forecast_avg <- NULL
  forecast_scalar <- NULL

  if (n_fc > 0L) {
    VaR_avg_path <- sum_VaR / n_fc
    ES_avg_path  <- if (risk == "ES") (sum_ES / n_fc) else NULL

    forecast_avg <- list(
      n_used = n_fc,
      VaR_avg = VaR_avg_path,
      ES_avg  = ES_avg_path
    )

    forecast_scalar <- list(
      n_used = n_fc,
      T = length(VaR_avg_path),
      VaR_mean = mean(VaR_avg_path, na.rm = TRUE),
      ES_mean  = if (!is.null(ES_avg_path)) mean(ES_avg_path, na.rm = TRUE) else NA_real_
    )
  }

  list(
    M_store = M_store,
    det_store = det_store,
    forecast_avg = forecast_avg,
    forecast_scalar = forecast_scalar,
    settings = list(
      world = world,
      risk = risk,
      f = f,
      interval_mode = interval_mode,
      preset = preset,
      weighting = weighting,
      tau = tau,
      lambda_max = lambda_max,
      wealth_kappa = wealth_kappa,
      cov_K = cov_K,
      cov_kappa = cov_kappa,
      forecast_mode = forecast_mode,
      refit_every = refit_every,
      oracle_mc_N = oracle_mc_N,
      segments = segments,
      ex7_under_factor = ex7_under_factor,
      ex7_theta = ex7_theta,
      ex7_noise_set = ex7_noise_set
    )
  )
}

# ============================================================
# 9) Summaries and plotting helpers
# ============================================================

summarize_mean_paths <- function(M_store) {
  methods <- names(M_store)
  T1 <- ncol(M_store[[1]])
  tvec <- 0:(T1 - 1)

  bind_rows(lapply(methods, function(m) {
    mat <- M_store[[m]]
    mat <- pmin(mat, 100)

    tibble(
      t = tvec,
      method = m,
      mean_logE = colMeans(mat, na.rm = TRUE),
      q10 = apply(mat, 2, quantile, probs = 0.10, na.rm = TRUE),
      q90 = apply(mat, 2, quantile, probs = 0.90, na.rm = TRUE)
    )
  }))
}

plot_avg_log_eprocess <- function(sum_df, title = NULL, ribbon = FALSE) {
  p <- ggplot(sum_df, aes(x = t, y = mean_logE, color = method)) +
    geom_line(linewidth = 1) +
    labs(x = "backtest time", y = "average log e-process", title = title) +
    theme_bw()

  if (ribbon) {
    p <- p +
      geom_ribbon(aes(ymin = q10, ymax = q90, fill = method), alpha = 0.15, color = NA) +
      guides(fill = "none")
  }

  p
}

summarize_detection <- function(det_store, T) {
  methods <- names(det_store)

  dplyr::bind_rows(lapply(methods, function(m) {
    mat <- det_store[[m]]

    tibble::tibble(
      method = m,
      det2  = mean(is.finite(mat[, "c2"]),  na.rm = TRUE),
      det5  = mean(is.finite(mat[, "c5"]),  na.rm = TRUE),
      det10 = mean(is.finite(mat[, "c10"]), na.rm = TRUE),
      arl2  = if (any(is.finite(mat[, "c2"])))  mean(mat[is.finite(mat[, "c2"]),  "c2"],  na.rm = TRUE) else Inf,
      arl5  = if (any(is.finite(mat[, "c5"])))  mean(mat[is.finite(mat[, "c5"]),  "c5"],  na.rm = TRUE) else Inf,
      arl10 = if (any(is.finite(mat[, "c10"]))) mean(mat[is.finite(mat[, "c10"]), "c10"], na.rm = TRUE) else Inf
    )
  }))
}

# Plot one simulated loss path together with its reported forecast(s).
plot_one_path_loss_forecast <- function(one_rep,
                                        risk = c("VaR", "ES"),
                                        title = NULL,
                                        show_base = TRUE,
                                        lw = 0.30,
                                        alpha_forecast = 0.95,
                                        alpha_base = 0.45,
                                        alpha_loss = 0.55) {
  risk <- match.arg(risk)

  L <- one_rep$series$L
  t <- seq_along(L)

  if (risk == "VaR") {
    df <- tibble::tibble(t = t, Loss = L, VaR = one_rep$series$VaR)
    if (show_base && !is.null(one_rep$series$VaR_base)) {
      df$VaR_base <- one_rep$series$VaR_base
    }
  } else {
    df <- tibble::tibble(
      t = t,
      Loss = L,
      ES = one_rep$series$ES,
      VaR = one_rep$series$VaR
    )
    if (show_base) {
      if (!is.null(one_rep$series$ES_base))  df$ES_base  <- one_rep$series$ES_base
      if (!is.null(one_rep$series$VaR_base)) df$VaR_base <- one_rep$series$VaR_base
    }
  }

  df_long <- df |>
    tidyr::pivot_longer(-t, names_to = "series", values_to = "value") |>
    dplyr::mutate(
      is_base = grepl("_base$", series),
      is_loss = series == "Loss",
      alpha_line = dplyr::case_when(
        is_loss ~ alpha_loss,
        is_base ~ alpha_base,
        TRUE    ~ alpha_forecast
      )
    )

  ggplot2::ggplot(df_long, ggplot2::aes(x = t, y = value, color = series, alpha = alpha_line)) +
    ggplot2::geom_line(linewidth = lw, lineend = "round") +
    ggplot2::labs(x = "backtest time", y = "", title = title, color = "") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::scale_alpha_identity()
}

# Plot the distortion multiplier used in a smooth-step or scale experiment.
plot_game_multiplier <- function(one_rep, title = NULL) {
  df <- tibble::tibble(t = 1:length(one_rep$dist_mult), mult = one_rep$dist_mult)

  ggplot2::ggplot(df, ggplot2::aes(t, mult)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "backtest time", y = "distortion multiplier", title = title)
}

# ============================================================
# 10) Guard presets
# ============================================================

# Strict guard: reacts quickly to disagreement across interval-wise lambdas.
var_guard_strict <- list(
  diff_thresh = 0.04,
  sd_thresh   = 0.03,
  compare_q       = 1L,
  mean_gap_thresh = 0.015,
  slope_thresh    = 0.004,
  min_keep = 1L,
  max_keep = 3L,
  prefix_rule = "diff",
  prefix_diff_thresh = 0.03
)

# Smooth guard: looser version that usually keeps more intervals.
var_guard_smooth <- list(
  diff_thresh = 0.07,
  sd_thresh   = 0.05,
  compare_q       = 2L,
  mean_gap_thresh = 0.03,
  slope_thresh    = 0.008,
  min_keep = 2L,
  max_keep = 6L,
  prefix_rule = "sd",
  prefix_sd_thresh = 0.03
)
