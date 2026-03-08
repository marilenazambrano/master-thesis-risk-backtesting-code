###############################################################################
# Tail risk e-backtesting
#
# This script implements the univariate tail-risk simulation study used in the
# thesis. 
#
# Included forecast modes:
#   (1) oracle_true      : uses the true latent mu/sigma and known t degrees
#                          of freedom
#   (2) oracle_empirical : resamples presample standardized residuals
#   (3) ml_t             : rolling AR(1)-GARCH(1,1) fit with standardized
#                          Student-t innovations
#
# Included tail-risk tests:
#   (1) ES
#   (2) TailSD
#   (3) MedianShortfall
#   (4) RVaR
#
# Sequential aggregation:
#   - componentwise GREM-style e-processes
#   - Bonferroni and mixture combinations
#   - additional c-threshold summaries for c in {2, 5, 10}
###############################################################################

suppressPackageStartupMessages({
  library(rugarch)
})

# ----------------------------
# 0) Global settings
# ----------------------------
set.seed(123)

T_train    <- 500
T_backtest <- 500
T_total    <- T_train + T_backtest

# Tail levels (loss right tail)
alpha_tail      <- 0.975
alpha_rvar_low  <- 0.95
alpha_rvar_high <- 0.975

# Global test level
# Used for the delta-based Bonferroni / mixture summaries. The thesis tables are
# reported mainly for c-thresholds, which are returned in addition.
delta_global <- 0.05

# Betting settings (GREM)
J_lambda   <- 500L
lambda_max <- 0.99

# Default distortion grid
f_set_default <- c(0.9, 1.0, 1.1)
apply_to_set  <- c("none", "r_only", "both")

# ============================================================
# 1) Simulation: univariate losses L_t with AR(1)-GARCH(1,1)
#    standardized-t innovations (Var=1) with df_innov
# ============================================================

simulate_world_tail <- function(T_total, df_innov = 4, burnin = 1000L) {
  if (!is.finite(df_innov) || df_innov <= 2) {
    stop("df_innov must be finite and > 2.")
  }

  T_all <- as.integer(T_total + burnin)

  L   <- numeric(T_all)
  mu  <- numeric(T_all)
  sig <- numeric(T_all)

  # GARCH(1,1) on the loss scale
  g_omega <- 0.01
  g_alpha <- 0.10
  g_beta  <- 0.85

  # AR(1) mean for losses
  mu0 <- -0.05
  phi <- 0.30

  # Draw standardized t innovations with unit variance.
  rstd_t <- function(n, df) {
    s <- sqrt((df - 2) / df)
    s * stats::rt(n, df = df)
  }

  sig[1] <- sqrt(g_omega / (1 - g_alpha - g_beta))
  mu[1]  <- mu0
  z1     <- rstd_t(1, df_innov)
  L[1]   <- mu[1] + sig[1] * z1

  for (t in 2:T_all) {
    mu[t] <- mu0 + phi * L[t - 1]

    z_prev <- (L[t - 1] - mu[t - 1]) / pmax(sig[t - 1], 1e-8)
    sig[t] <- sqrt(pmax(
      g_omega + g_alpha * sig[t - 1]^2 * z_prev^2 + g_beta * sig[t - 1]^2,
      0
    ))

    zt   <- rstd_t(1, df_innov)
    L[t] <- mu[t] + sig[t] * zt
  }

  idx <- (burnin + 1L):T_all
  list(L = L[idx], mu = mu[idx], sig = sig[idx], df_innov = df_innov)
}

get_backtest_slice_tail <- function(sim, T_train, T_backtest) {
  idx_bt <- (T_train + 1L):(T_train + T_backtest)
  list(
    idx_bt = idx_bt,
    L_bt   = sim$L[idx_bt],
    mu_bt  = sim$mu[idx_bt],
    sig_bt = sim$sig[idx_bt]
  )
}

# ============================================================
# 2) ML forecasting: rolling refit AR(1)-GARCH(1,1) with Student-t
# ============================================================

fit_ar1_garch11_t <- function(x) {
  spec <- rugarch::ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model     = list(armaOrder = c(1, 0), include.mean = TRUE),
    distribution.model = "std"
  )
  rugarch::ugarchfit(spec = spec, data = x, solver = "hybrid")
}

build_ml_mu_sig_t <- function(series_full, T_train, T_backtest, refit_every = 10L) {
  refit_every <- as.integer(refit_every)

  mu_f  <- rep(NA_real_, T_backtest)
  sig_f <- rep(NA_real_, T_backtest)
  nu_f  <- rep(NA_real_, T_backtest)

  fits <- vector("list", length = ceiling(T_backtest / refit_every))
  fit_idx <- 0L

  k <- 1L
  while (k <= T_backtest) {
    t_prev_full <- T_train + (k - 1L)
    win_start <- t_prev_full - T_train + 1L
    win_end   <- t_prev_full
    x_win <- series_full[win_start:win_end]

    fit <- tryCatch(fit_ar1_garch11_t(x_win), error = function(e) NULL)
    if (is.null(fit)) {
      stop("ugarchfit failed (Student-t). Try a different seed or refit_every.")
    }

    fit_idx <- fit_idx + 1L
    fits[[fit_idx]] <- fit

    cf <- coef(fit)
    mu0 <- if ("mu" %in% names(cf)) as.numeric(cf[["mu"]]) else 0
    ar1_name <- intersect(names(cf), c("ar1", "ar[1]", "ar"))
    ar1 <- if (length(ar1_name)) as.numeric(cf[[ar1_name[1]]]) else 0
    omg <- if ("omega"  %in% names(cf)) as.numeric(cf[["omega"]])  else NA_real_
    a1  <- if ("alpha1" %in% names(cf)) as.numeric(cf[["alpha1"]]) else NA_real_
    b1  <- if ("beta1"  %in% names(cf)) as.numeric(cf[["beta1"]])  else NA_real_
    nu_block <- if ("shape" %in% names(cf)) as.numeric(cf[["shape"]]) else NA_real_

    sig_series <- as.numeric(rugarch::sigma(fit))
    mu_series  <- as.numeric(rugarch::fitted(fit))

    prev <- series_full[t_prev_full]
    sig_prev <- sig_series[length(sig_series)]
    mu_prev  <- mu_series[length(mu_series)]
    res_prev <- prev - mu_prev

    block_len <- min(refit_every, T_backtest - k + 1L)

    for (j in seq_len(block_len)) {
      t_full <- t_prev_full + 1L

      mu1  <- mu0 + ar1 * prev
      sig1 <- sqrt(pmax(omg + a1 * (res_prev^2) + b1 * (sig_prev^2), 0))

      mu_f[k]  <- mu1
      sig_f[k] <- sig1
      nu_f[k]  <- nu_block

      prev     <- series_full[t_full]
      mu_prev  <- mu1
      sig_prev <- sig1
      res_prev <- prev - mu_prev

      t_prev_full <- t_full
      k <- k + 1L
      if (k > T_backtest) break
    }
  }

  list(mu = mu_f, sig = sig_f, nu = nu_f, fits = fits[seq_len(fit_idx)])
}

# ============================================================
# 3) Helpers: analytic VaR / ES 
# ============================================================

q_std_t <- function(p, df) {
  s <- sqrt((df - 2) / df)
  s * stats::qt(p, df = df)
}

es_std_t_righttail <- function(alpha, df) {
  # For T ~ t_df (non-standardized), right-tail ES is
  # E[T | T > q] = ((df + q^2) / (df - 1)) * f(q) / (1 - alpha).
  # We then multiply by s = sqrt((df - 2) / df) to standardize variance to one.
  q  <- stats::qt(alpha, df = df)
  fq <- stats::dt(q, df = df)
  ES_T <- ((df + q^2) / (df - 1)) * (fq / (1 - alpha))
  s <- sqrt((df - 2) / df)
  s * ES_T
}

# If M solves P(L <= M | L > VaR_alpha) = 0.5, then F(M) = (1 + alpha) / 2.
cond_median_from_Finv <- function(alpha, Finv) {
  Finv((1 + alpha) / 2)
}

# Conditional draw from a right tail using inverse CDF sampling.
cond_draw_from_Finv <- function(N, alpha, Finv) {
  u <- stats::runif(N, min = alpha, max = 1)
  Finv(u)
}

mean_or_na <- function(x) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

col_means_or_na <- function(mat) {
  out <- vapply(seq_len(ncol(mat)), function(j) mean_or_na(mat[, j]), numeric(1))
  stats::setNames(out, colnames(mat))
}

normalize_apply_to_tail <- function(apply_to) {
  apply_to <- as.character(apply_to)[1]

  if (identical(apply_to, "all")) {
    return("both")
  }

  if (identical(apply_to, "aux_only")) {
    stop("'aux_only' is not part of the thesis tail-risk setup. Use 'none', 'r_only', or 'both'.")
  }

  if (!apply_to %in% apply_to_set) {
    stop("apply_to must be one of: ", paste(apply_to_set, collapse = ", "))
  }

  apply_to
}

# ============================================================
# 4) Forecast construction 
# ============================================================

needs_for_test <- function(test = c("ES", "TailSD", "RVaR", "MedianShortfall")) {
  test <- match.arg(test)

  if (test == "ES") {
    return(list(var_alpha = TRUE, es = TRUE, tailsd = FALSE, rvar = FALSE, med = FALSE))
  }
  if (test == "TailSD") {
    return(list(var_alpha = TRUE, es = TRUE, tailsd = TRUE, rvar = FALSE, med = FALSE))
  }
  if (test == "RVaR") {
    return(list(var_alpha = FALSE, es = FALSE, tailsd = FALSE, rvar = TRUE, med = FALSE))
  }

  list(var_alpha = TRUE, es = FALSE, tailsd = FALSE, rvar = FALSE, med = TRUE)
}

forecast_tail_oracle_true <- function(sim, T_train, T_backtest,
                                      test,
                                      alpha_tail, alpha_rvar_low, alpha_rvar_high,
                                      N_cond = 3000L) {
  need <- needs_for_test(test)
  sl <- get_backtest_slice_tail(sim, T_train, T_backtest)
  mu  <- sl$mu_bt
  sig <- sl$sig_bt

  df <- sim$df_innov
  if (!is.finite(df) || df <= 2) {
    stop("oracle_true requires df_innov > 2.")
  }

  out <- list(
    mu = mu,
    sig = sig,
    df = df,
    mode = "oracle_true",
    test = test,
    params = list(
      alpha_tail = alpha_tail,
      alpha_low = alpha_rvar_low,
      alpha_high = alpha_rvar_high,
      N_cond = as.integer(N_cond)
    )
  )

  if (need$var_alpha || need$es || need$tailsd || need$med) {
    qA <- q_std_t(alpha_tail, df)
    out$VaR_a <- mu + sig * qA
  }

  if (need$es || need$tailsd) {
    esA <- es_std_t_righttail(alpha_tail, df)
    out$ES_a <- mu + sig * esA
  }

  if (need$tailsd) {
    N_cond <- as.integer(N_cond)
    FinvA <- function(u) q_std_t(u, df)
    Zc <- cond_draw_from_Finv(N_cond, alpha_tail, FinvA)

    TailSD <- numeric(T_backtest)
    for (t in seq_len(T_backtest)) {
      Lc <- mu[t] + sig[t] * Zc
      TailSD[t] <- sqrt(mean((Lc - out$ES_a[t])^2))
    }
    out$TailSD <- TailSD
  }

  if (need$med) {
    Finv <- function(u) q_std_t(u, df)
    medZ <- cond_median_from_Finv(alpha_tail, Finv)
    out$MedShortfall <- mu + sig * medZ
  }

  if (need$rvar) {
    u_grid <- seq(alpha_rvar_low, alpha_rvar_high, length.out = 51)
    q_grid <- q_std_t(u_grid, df)
    out$VaR_low  <- mu + sig * q_std_t(alpha_rvar_low, df)
    out$VaR_high <- mu + sig * q_std_t(alpha_rvar_high, df)
    out$RVaR     <- mu + sig * mean(q_grid)
  }

  out
}

forecast_tail_oracle_empirical <- function(sim, T_train, T_backtest,
                                           test,
                                           alpha_tail, alpha_rvar_low, alpha_rvar_high,
                                           N_mc = 5000L) {
  need <- needs_for_test(test)
  sl <- get_backtest_slice_tail(sim, T_train, T_backtest)
  mu  <- sl$mu_bt
  sig <- sl$sig_bt

  # Presample standardized residuals using the true latent states.
  idx_pre <- seq_len(as.integer(T_train))
  z_hat <- (sim$L[idx_pre] - sim$mu[idx_pre]) / pmax(sim$sig[idx_pre], 1e-8)
  z_hat <- z_hat[is.finite(z_hat)]
  if (length(z_hat) < 50) {
    stop("oracle_empirical: too few presample residuals.")
  }

  N_mc <- as.integer(N_mc)
  out <- list(
    mu = mu,
    sig = sig,
    mode = "oracle_empirical",
    test = test,
    params = list(
      alpha_tail = alpha_tail,
      alpha_low = alpha_rvar_low,
      alpha_high = alpha_rvar_high,
      N_mc = N_mc
    )
  )

  if (need$var_alpha || need$es || need$tailsd || need$med) {
    out$VaR_a <- rep(NA_real_, T_backtest)
  }
  if (need$es || need$tailsd) out$ES_a <- rep(NA_real_, T_backtest)
  if (need$tailsd) out$TailSD <- rep(NA_real_, T_backtest)
  if (need$med) out$MedShortfall <- rep(NA_real_, T_backtest)
  if (need$rvar) {
    out$VaR_low  <- rep(NA_real_, T_backtest)
    out$VaR_high <- rep(NA_real_, T_backtest)
    out$RVaR     <- rep(NA_real_, T_backtest)
  }

  for (t in seq_len(T_backtest)) {
    z <- sample(z_hat, size = N_mc, replace = TRUE)
    Ls <- mu[t] + sig[t] * z

    if (need$var_alpha || need$es || need$tailsd || need$med) {
      out$VaR_a[t] <- as.numeric(stats::quantile(Ls, probs = alpha_tail, type = 8))
    }

    if (need$es || need$tailsd) {
      tail <- Ls[Ls > out$VaR_a[t]]
      if (length(tail) < 30) tail <- Ls
      out$ES_a[t] <- mean(tail)
    }

    if (need$tailsd) {
      tail <- Ls[Ls > out$VaR_a[t]]
      if (length(tail) < 30) tail <- Ls
      out$TailSD[t] <- sqrt(mean((tail - out$ES_a[t])^2))
      if (!is.finite(out$TailSD[t]) || out$TailSD[t] <= 0) {
        out$TailSD[t] <- NA_real_
      }
    }

    if (need$med) {
      tail <- Ls[Ls > out$VaR_a[t]]
      if (length(tail) < 30) tail <- Ls
      out$MedShortfall[t] <- as.numeric(stats::median(tail))
    }

    if (need$rvar) {
      out$VaR_low[t]  <- as.numeric(stats::quantile(Ls, probs = alpha_rvar_low, type = 8))
      out$VaR_high[t] <- as.numeric(stats::quantile(Ls, probs = alpha_rvar_high, type = 8))
      u_grid <- seq(alpha_rvar_low, alpha_rvar_high, length.out = 15)
      qs <- as.numeric(stats::quantile(Ls, probs = u_grid, type = 8))
      out$RVaR[t] <- mean(qs)
    }
  }

  out
}

forecast_tail_ml_t <- function(sim, T_train, T_backtest,
                               test,
                               alpha_tail, alpha_rvar_low, alpha_rvar_high,
                               N_cond = 3000L,
                               refit_every = 10L) {
  need <- needs_for_test(test)

  st <- build_ml_mu_sig_t(sim$L, T_train, T_backtest, refit_every = refit_every)
  mu <- st$mu
  sig <- st$sig
  nu <- st$nu

  out <- list(
    mu = mu,
    sig = sig,
    nu = nu,
    mode = "ml_t",
    test = test,
    params = list(
      alpha_tail = alpha_tail,
      alpha_low = alpha_rvar_low,
      alpha_high = alpha_rvar_high,
      N_cond = as.integer(N_cond),
      refit_every = as.integer(refit_every)
    )
  )

  Tbt <- T_backtest
  block_id <- ceiling(seq_len(Tbt) / as.integer(refit_every))
  blocks <- sort(unique(block_id))

  if (need$var_alpha || need$es || need$tailsd || need$med) {
    out$VaR_a <- rep(NA_real_, Tbt)
  }
  if (need$es || need$tailsd) out$ES_a <- rep(NA_real_, Tbt)
  if (need$tailsd) out$TailSD <- rep(NA_real_, Tbt)
  if (need$med) out$MedShortfall <- rep(NA_real_, Tbt)
  if (need$rvar) {
    out$VaR_low  <- rep(NA_real_, Tbt)
    out$VaR_high <- rep(NA_real_, Tbt)
    out$RVaR     <- rep(NA_real_, Tbt)
  }

  N_cond <- as.integer(N_cond)

  for (b in blocks) {
    idx_b <- which(block_id == b)

    # Within a refit block the fitted shape parameter is typically constant.
    nu_b <- suppressWarnings(stats::median(nu[idx_b], na.rm = TRUE))
    if (!is.finite(nu_b) || nu_b <= 2) nu_b <- 8

    Finv <- function(u) q_std_t(u, nu_b)

    qA   <- q_std_t(alpha_tail, nu_b)
    esA  <- es_std_t_righttail(alpha_tail, nu_b)
    medZ <- q_std_t((1 + alpha_tail) / 2, nu_b)
    Zc   <- cond_draw_from_Finv(N_cond, alpha_tail, Finv)

    u_grid_r <- seq(alpha_rvar_low, alpha_rvar_high, length.out = 51)
    q_grid_r <- q_std_t(u_grid_r, nu_b)
    qL <- q_std_t(alpha_rvar_low, nu_b)
    qH <- q_std_t(alpha_rvar_high, nu_b)

    for (t in idx_b) {
      if (!is.finite(mu[t]) || !is.finite(sig[t]) || sig[t] <= 0) next

      if (need$var_alpha || need$es || need$tailsd || need$med) {
        out$VaR_a[t] <- mu[t] + sig[t] * qA
      }
      if (need$es || need$tailsd) {
        out$ES_a[t] <- mu[t] + sig[t] * esA
      }
      if (need$med) {
        out$MedShortfall[t] <- mu[t] + sig[t] * medZ
      }
      if (need$tailsd) {
        Lc <- mu[t] + sig[t] * Zc
        out$TailSD[t] <- sqrt(mean((Lc - out$ES_a[t])^2))
      }
      if (need$rvar) {
        out$VaR_low[t]  <- mu[t] + sig[t] * qL
        out$VaR_high[t] <- mu[t] + sig[t] * qH
        out$RVaR[t]     <- mu[t] + sig[t] * mean(q_grid_r)
      }
    }
  }

  out
}

forecast_tail_dispatch <- function(sim,
                                   mode,
                                   T_train, T_backtest,
                                   test,
                                   alpha_tail, alpha_rvar_low, alpha_rvar_high,
                                   N_mc = 5000L,
                                   N_cond = 3000L,
                                   refit_every = 10L) {
  mode <- match.arg(mode, choices = c("oracle_true", "oracle_empirical", "ml_t"))
  test <- match.arg(test, choices = c("ES", "TailSD", "RVaR", "MedianShortfall"))

  if (mode == "oracle_true") {
    return(forecast_tail_oracle_true(
      sim, T_train, T_backtest, test,
      alpha_tail, alpha_rvar_low, alpha_rvar_high,
      N_cond = N_cond
    ))
  }

  if (mode == "oracle_empirical") {
    return(forecast_tail_oracle_empirical(
      sim, T_train, T_backtest, test,
      alpha_tail, alpha_rvar_low, alpha_rvar_high,
      N_mc = N_mc
    ))
  }

  forecast_tail_ml_t(
    sim, T_train, T_backtest, test,
    alpha_tail, alpha_rvar_low, alpha_rvar_high,
    N_cond = N_cond,
    refit_every = refit_every
  )
}

# ============================================================
# 5) Distortion map (none / r_only / both)
# ============================================================

scale_abs <- function(x, f) {
  x + (f - 1) * abs(x)
}

apply_distortion_tail <- function(fc, f = 1.0,
                                  apply_to = "none",
                                  test = c("ES", "TailSD", "RVaR", "MedianShortfall")) {
  apply_to <- normalize_apply_to_tail(apply_to)
  test <- match.arg(test)

  if (!is.finite(f) || f <= 0) {
    stop("f must be > 0.")
  }

  out <- fc
  out$dist <- list(f = f, apply_to = apply_to, test = test)

  if (apply_to == "none") return(out)

  scale_aux <- identical(apply_to, "both")
  scale_r   <- apply_to %in% c("r_only", "both")

  if (test == "ES") {
    if (scale_aux && !is.null(out$VaR_a)) out$VaR_a <- scale_abs(out$VaR_a, f)
    if (scale_r   && !is.null(out$ES_a))  out$ES_a  <- scale_abs(out$ES_a, f)
    return(out)
  }

  if (test == "TailSD") {
    if (scale_aux) {
      if (!is.null(out$VaR_a)) out$VaR_a <- scale_abs(out$VaR_a, f)
      if (!is.null(out$ES_a))  out$ES_a  <- scale_abs(out$ES_a, f)
    }
    if (scale_r && !is.null(out$TailSD)) {
      out$TailSD <- scale_abs(out$TailSD, f)
    }
    return(out)
  }

  if (test == "MedianShortfall") {
    if (scale_aux && !is.null(out$VaR_a)) {
      out$VaR_a <- scale_abs(out$VaR_a, f)
    }
    if (scale_r && !is.null(out$MedShortfall)) {
      out$MedShortfall <- scale_abs(out$MedShortfall, f)
    }
    return(out)
  }

  # RVaR
  if (scale_aux) {
    if (!is.null(out$VaR_low))  out$VaR_low  <- scale_abs(out$VaR_low,  f)
    if (!is.null(out$VaR_high)) out$VaR_high <- scale_abs(out$VaR_high, f)
  }
  if (scale_r && !is.null(out$RVaR)) {
    out$RVaR <- scale_abs(out$RVaR, f)
  }

  out
}

# ============================================================
# 6) One-step e-values for tail components
# ============================================================

# Two one-sided VaR guard components used inside the multicomponent tests.
e_var_under <- function(L, z, alpha) {
  as.numeric(L > z) / (1 - alpha)
}

e_var_over <- function(L, z, alpha) {
  as.numeric(L <= z) / alpha
}

# ES e-value (canonical, right tail)
e_es <- function(L, ES, VaR, alpha) {
  den <- (1 - alpha) * (ES - VaR)
  if (!is.finite(den) || den <= 0) return(Inf)
  pmax(L - VaR, 0) / den
}

# TailSD e-value: quadratic block using VaR and ES auxiliaries.
e_tailsd_quad <- function(L, r, VaR, ES, alpha) {
  if (!is.finite(r) || r <= 0) return(Inf)
  den <- (1 - alpha) * (r^2)
  if (!is.finite(den) || den <= 0) return(Inf)
  num <- ((L - ES) * as.numeric(L > VaR))^2
  num / den
}

# Median shortfall conditional-median inequality e-value.
e_median_shortfall <- function(L, r, VaR, alpha) {
  den <- 0.5 * (1 - alpha)
  if (!is.finite(den) || den <= 0) return(Inf)
  as.numeric(L > max(r, VaR)) / den
}

# RVaR e-value 
S_gamma <- function(gamma, z, x) {
  I <- as.numeric(x <= z)
  (I - gamma) * z - I * x
}

V3_rvar <- function(alpha, beta, z1, z2, r, x) {
  r + (S_gamma(beta, z2, x) - S_gamma(alpha, z1, x)) / (beta - alpha)
}

M_rvar_closed <- function(alpha, beta, z1, z2, r) {
  D_left <- (1 - beta) * z2 - (1 - alpha) * z1
  r + D_left / (beta - alpha)
}

e_rvar_q <- function(L, alpha, beta, z1, z2, r) {
  if (!is.finite(alpha) || !is.finite(beta) || alpha <= 0 || beta >= 1 || alpha >= beta) {
    return(Inf)
  }
  if (!is.finite(z1) || !is.finite(z2) || !is.finite(r)) {
    return(Inf)
  }
  if (z2 <= z1) return(Inf)

  M <- M_rvar_closed(alpha, beta, z1, z2, r)
  if (!is.finite(M) || M <= 0) return(1)

  V3x <- V3_rvar(alpha, beta, z1, z2, r, L)
  pmax(1 - V3x / M, 0)
}

# ============================================================
# 7) GREL 
# ============================================================

past_e_grel_tail_component <- function(component, idx0, t, L, fc,
                                       alpha_tail, alpha_low, alpha_high) {
  Lp <- L[idx0]

  if (component == "VaR_under") {
    zt <- fc$VaR_a[t]
    return(as.numeric(Lp > zt) / (1 - alpha_tail))
  }
  if (component == "VaR_over") {
    zt <- fc$VaR_a[t]
    return(as.numeric(Lp <= zt) / alpha_tail)
  }

  if (component == "ES") {
    return(sapply(Lp, function(x) {
      e_es(x, ES = fc$ES_a[t], VaR = fc$VaR_a[t], alpha = alpha_tail)
    }))
  }

  if (component == "TailSD") {
    return(sapply(Lp, function(x) {
      e_tailsd_quad(
        x,
        r = fc$TailSD[t],
        VaR = fc$VaR_a[t],
        ES = fc$ES_a[t],
        alpha = alpha_tail
      )
    }))
  }

  if (component == "MedSF") {
    return(sapply(Lp, function(x) {
      e_median_shortfall(
        x,
        r = fc$MedShortfall[t],
        VaR = fc$VaR_a[t],
        alpha = alpha_tail
      )
    }))
  }

  if (component == "VaRlow_under") {
    zt <- fc$VaR_low[t]
    return(as.numeric(Lp > zt) / (1 - alpha_low))
  }
  if (component == "VaRlow_over") {
    zt <- fc$VaR_low[t]
    return(as.numeric(Lp <= zt) / alpha_low)
  }
  if (component == "VaRhigh_under") {
    zt <- fc$VaR_high[t]
    return(as.numeric(Lp > zt) / (1 - alpha_high))
  }
  if (component == "VaRhigh_over") {
    zt <- fc$VaR_high[t]
    return(as.numeric(Lp <= zt) / alpha_high)
  }
  if (component == "RVaR") {
    return(sapply(Lp, function(x) {
      e_rvar_q(
        x,
        alpha = alpha_low,
        beta = alpha_high,
        z1 = fc$VaR_low[t],
        z2 = fc$VaR_high[t],
        r = fc$RVaR[t]
      )
    }))
  }

  stop("Unknown component: ", component)
}

# ============================================================
# 8) GREM core (S0): GREE/GREL mixture
# ============================================================

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

update_log_wealth <- function(M_prev, e_t, lambda_t) {
  if (!is.finite(M_prev)) return(Inf)
  if (!is.finite(lambda_t) || lambda_t < 0) lambda_t <- 0
  if (lambda_t == 0) return(M_prev)
  if (!is.finite(e_t)) return(Inf)

  inc <- 1 - lambda_t + lambda_t * e_t
  if (!is.finite(inc)) return(Inf)
  if (inc <= 0) return(-Inf)
  M_prev + log(inc)
}

weight_from_logs <- function(logA, logB) {
  if (is.infinite(logA) && logA < 0 && is.infinite(logB) && logB < 0) return(0.5)
  if (is.infinite(logA) && logA < 0) return(0)
  if (is.infinite(logB) && logB < 0) return(1)
  1 / (1 + exp(logB - logA))
}

run_grem_S0_gree_grel_tail <- function(component, e_series, L, fc,
                                       alpha_tail, alpha_low, alpha_high,
                                       J_lambda = 500L, lambda_max = 0.99) {
  T <- length(e_series)
  J <- as.integer(J_lambda)
  if (J <= 0) stop("J_lambda must be positive.")

  M_E   <- numeric(T + 1L); M_E[1]   <- 0
  M_L   <- numeric(T + 1L); M_L[1]   <- 0
  M_mix <- numeric(T + 1L); M_mix[1] <- 0

  lam_E <- numeric(T)
  lam_L <- numeric(T)
  lam_M <- numeric(T)

  for (t in seq_len(T)) {
    past_end <- t - 1L

    if (past_end <= 0L) {
      lam_E[t] <- 0
      lam_L[t] <- 0
    } else {
      idx0 <- max(1L, past_end - J + 1L):past_end

      past_GREE <- e_series[idx0]
      past_GREL <- past_e_grel_tail_component(
        component, idx0, t, L, fc,
        alpha_tail, alpha_low, alpha_high
      )

      lam_E[t] <- lambda_from_e_taylor(past_GREE, lambda_max = lambda_max)
      lam_L[t] <- lambda_from_e_taylor(past_GREL, lambda_max = lambda_max)
    }

    M_E[t + 1L] <- update_log_wealth(M_E[t], e_series[t], lam_E[t])
    M_L[t + 1L] <- update_log_wealth(M_L[t], e_series[t], lam_L[t])

    wE <- weight_from_logs(M_E[t], M_L[t])
    lam_M[t] <- wE * lam_E[t] + (1 - wE) * lam_L[t]
    M_mix[t + 1L] <- update_log_wealth(M_mix[t], e_series[t], lam_M[t])
  }

  W <- exp(pmin(M_mix, 700))
  list(
    logW = M_mix,
    W = W,
    lambda = lam_M,
    lambda_GREE = lam_E,
    lambda_GREL = lam_L,
    logW_GREE = M_E,
    logW_GREL = M_L,
    e = e_series,
    component = component
  )
}

# ============================================================
# 9) Combine components within one test
# ============================================================

logsumexp <- function(v) {
  m <- max(v)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(v - m)))
}

combine_bonferroni <- function(logW_list, delta = 0.05, weights = NULL) {
  comp_names <- names(logW_list)
  K <- length(comp_names)
  if (is.null(weights)) weights <- rep(1 / K, K)
  weights <- weights / sum(weights)
  delta_i <- delta * weights

  hit_t <- rep(Inf, K)
  for (k in seq_len(K)) {
    thr <- log(1 / delta_i[k])
    hit <- which(logW_list[[k]] >= thr)
    if (length(hit)) hit_t[k] <- hit[1] - 1L
  }

  t_star <- min(hit_t)
  which_min <- which(hit_t == t_star)
  list(
    detected = is.finite(t_star),
    t = t_star,
    who = if (is.finite(t_star)) comp_names[which_min] else character(0),
    per_component = data.frame(
      component = comp_names,
      weight = weights,
      delta_i = delta_i,
      detect_time = hit_t
    )
  )
}

combine_mixture <- function(logW_list, delta = 0.05, weights = NULL) {
  comp_names <- names(logW_list)
  K <- length(comp_names)
  if (is.null(weights)) weights <- rep(1 / K, K)
  weights <- weights / sum(weights)

  T <- length(logW_list[[1]]) - 1L
  logWmix <- rep(-Inf, T + 1L)
  for (t in seq_len(T + 1L)) {
    logs <- sapply(seq_len(K), function(k) log(weights[k]) + logW_list[[k]][t])
    logWmix[t] <- logsumexp(logs)
  }

  thr <- log(1 / delta)
  hit <- which(logWmix >= thr)
  t_star <- if (length(hit)) hit[1] - 1L else Inf

  list(detected = is.finite(t_star), t = t_star, logWmix = logWmix)
}

combine_bonferroni_c <- function(logW_list, c = 10, weights = NULL) {
  comp_names <- names(logW_list)
  K <- length(comp_names)
  if (is.null(weights)) weights <- rep(1 / K, K)
  weights <- weights / sum(weights)

  thr_k <- (1 / weights) * log(c)

  hit_t <- rep(Inf, K)
  for (k in seq_len(K)) {
    hit <- which(logW_list[[k]] >= thr_k[k])
    if (length(hit)) hit_t[k] <- hit[1] - 1L
  }

  t_star <- min(hit_t)
  which_min <- which(hit_t == t_star)

  list(
    detected = is.finite(t_star),
    t = t_star,
    who = if (is.finite(t_star)) comp_names[which_min] else character(0),
    per_component = data.frame(
      component = comp_names,
      weight = weights,
      thr_log = thr_k,
      detect_time = hit_t
    ),
    c = c
  )
}

combine_mixture_c <- function(logW_list, c = 10, weights = NULL) {
  comp_names <- names(logW_list)
  K <- length(comp_names)
  if (is.null(weights)) weights <- rep(1 / K, K)
  weights <- weights / sum(weights)

  T <- length(logW_list[[1]]) - 1L
  logWmix <- rep(-Inf, T + 1L)
  for (t in seq_len(T + 1L)) {
    logs <- sapply(seq_len(K), function(k) log(weights[k]) + logW_list[[k]][t])
    logWmix[t] <- logsumexp(logs)
  }

  thr <- log(c)
  hit <- which(logWmix >= thr)
  t_star <- if (length(hit)) hit[1] - 1L else Inf

  list(detected = is.finite(t_star), t = t_star, logWmix = logWmix, c = c)
}

# ============================================================
# 10) One replication for one tail-risk test
# ============================================================

components_for_test_tail <- function(test = c("ES", "TailSD", "RVaR", "MedianShortfall")) {
  test <- match.arg(test)

  if (test == "ES") {
    return(c("VaR_under", "VaR_over", "ES"))
  }
  if (test == "TailSD") {
    return(c("VaR_under", "VaR_over", "ES", "TailSD"))
  }
  if (test == "MedianShortfall") {
    return(c("VaR_under", "VaR_over", "MedSF"))
  }

  # RVaR needs both quantiles correct plus the RVaR target component.
  c("VaRlow_under", "VaRlow_over", "VaRhigh_under", "VaRhigh_over", "RVaR")
}

run_one_rep_tail <- function(test = c("ES", "TailSD", "RVaR", "MedianShortfall"),
                             mode = c("oracle_true", "oracle_empirical", "ml_t"),
                             f = 1.0,
                             apply_to = "none",
                             delta = 0.05,
                             N_mc = 5000L,
                             N_cond = 3000L,
                             refit_every = 10L,
                             weights = NULL,
                             seed = 1) {
  test <- match.arg(test)
  mode <- match.arg(mode)
  apply_to <- normalize_apply_to_tail(apply_to)

  set.seed(seed)

  sim <- simulate_world_tail(T_total, df_innov = 4)
  sl  <- get_backtest_slice_tail(sim, T_train, T_backtest)
  L   <- sl$L_bt

  fc <- forecast_tail_dispatch(
    sim,
    mode = mode,
    T_train = T_train,
    T_backtest = T_backtest,
    test = test,
    alpha_tail = alpha_tail,
    alpha_rvar_low = alpha_rvar_low,
    alpha_rvar_high = alpha_rvar_high,
    N_mc = N_mc,
    N_cond = N_cond,
    refit_every = refit_every
  )

  fc <- apply_distortion_tail(fc, f = f, apply_to = apply_to, test = test)

  comps <- components_for_test_tail(test)
  e_series <- list()

  if ("VaR_under" %in% comps) {
    e_series$VaR_under <- sapply(seq_len(T_backtest), function(t) {
      e_var_under(L[t], fc$VaR_a[t], alpha_tail)
    })
  }
  if ("VaR_over" %in% comps) {
    e_series$VaR_over <- sapply(seq_len(T_backtest), function(t) {
      e_var_over(L[t], fc$VaR_a[t], alpha_tail)
    })
  }
  if ("ES" %in% comps) {
    e_series$ES <- sapply(seq_len(T_backtest), function(t) {
      e_es(L[t], ES = fc$ES_a[t], VaR = fc$VaR_a[t], alpha = alpha_tail)
    })
  }
  if ("TailSD" %in% comps) {
    e_series$TailSD <- sapply(seq_len(T_backtest), function(t) {
      e_tailsd_quad(
        L[t],
        r = fc$TailSD[t],
        VaR = fc$VaR_a[t],
        ES = fc$ES_a[t],
        alpha = alpha_tail
      )
    })
  }
  if ("MedSF" %in% comps) {
    e_series$MedSF <- sapply(seq_len(T_backtest), function(t) {
      e_median_shortfall(
        L[t],
        r = fc$MedShortfall[t],
        VaR = fc$VaR_a[t],
        alpha = alpha_tail
      )
    })
  }

  if ("VaRlow_under" %in% comps) {
    e_series$VaRlow_under <- sapply(seq_len(T_backtest), function(t) {
      e_var_under(L[t], fc$VaR_low[t], alpha_rvar_low)
    })
  }
  if ("VaRlow_over" %in% comps) {
    e_series$VaRlow_over <- sapply(seq_len(T_backtest), function(t) {
      e_var_over(L[t], fc$VaR_low[t], alpha_rvar_low)
    })
  }
  if ("VaRhigh_under" %in% comps) {
    e_series$VaRhigh_under <- sapply(seq_len(T_backtest), function(t) {
      e_var_under(L[t], fc$VaR_high[t], alpha_rvar_high)
    })
  }
  if ("VaRhigh_over" %in% comps) {
    e_series$VaRhigh_over <- sapply(seq_len(T_backtest), function(t) {
      e_var_over(L[t], fc$VaR_high[t], alpha_rvar_high)
    })
  }
  if ("RVaR" %in% comps) {
    e_series$RVaR <- sapply(seq_len(T_backtest), function(t) {
      e_rvar_q(
        L[t],
        alpha = alpha_rvar_low,
        beta = alpha_rvar_high,
        z1 = fc$VaR_low[t],
        z2 = fc$VaR_high[t],
        r = fc$RVaR[t]
      )
    })
  }

  proc <- lapply(names(e_series), function(nm) {
    run_grem_S0_gree_grel_tail(
      component = nm,
      e_series = e_series[[nm]],
      L = L,
      fc = fc,
      alpha_tail = alpha_tail,
      alpha_low = alpha_rvar_low,
      alpha_high = alpha_rvar_high,
      J_lambda = J_lambda,
      lambda_max = lambda_max
    )
  })
  names(proc) <- names(e_series)
  logW_list <- lapply(proc, function(p) p$logW)

  comp_names <- names(logW_list)
  K <- length(comp_names)
  if (is.null(weights)) {
    w <- rep(1 / K, K)
    names(w) <- comp_names
  } else {
    if (is.null(names(weights))) {
      stop("weights must be a named numeric vector matching component names.")
    }
    w <- weights[comp_names]
    if (any(!is.finite(w))) stop("weights missing or non-finite for some components.")
    if (any(w < 0)) stop("weights must be nonnegative.")
    if (sum(w) <= 0) stop("weights sum must be > 0.")
    w <- w / sum(w)
  }

  bonf <- combine_bonferroni(logW_list, delta = delta, weights = w)
  mix  <- combine_mixture(logW_list, delta = delta, weights = w)

  bonf_c2  <- combine_bonferroni_c(logW_list, c = 2,  weights = w)
  bonf_c5  <- combine_bonferroni_c(logW_list, c = 5,  weights = w)
  bonf_c10 <- combine_bonferroni_c(logW_list, c = 10, weights = w)

  mix_c2   <- combine_mixture_c(logW_list, c = 2,  weights = w)
  mix_c5   <- combine_mixture_c(logW_list, c = 5,  weights = w)
  mix_c10  <- combine_mixture_c(logW_list, c = 10, weights = w)

  avg_fc <- list(
    avg_VaR_a    = mean_or_na(fc$VaR_a),
    avg_ES_a     = mean_or_na(fc$ES_a),
    avg_TailSD   = mean_or_na(fc$TailSD),
    avg_MedSF    = mean_or_na(fc$MedShortfall),
    avg_VaR_low  = mean_or_na(fc$VaR_low),
    avg_VaR_high = mean_or_na(fc$VaR_high),
    avg_RVaR     = mean_or_na(fc$RVaR)
  )

  list(
    meta = list(
      test = test,
      mode = mode,
      f = f,
      apply_to = apply_to,
      delta = delta,
      N_mc = N_mc,
      N_cond = N_cond,
      refit_every = refit_every
    ),
    forecasts = fc,
    avg_forecasts = avg_fc,
    components = list(e = e_series, proc = proc, logW = logW_list, weights = w),
    results = list(
      bonferroni = bonf,
      mixture = mix,
      bonf_c = list(c2 = bonf_c2, c5 = bonf_c5, c10 = bonf_c10),
      mix_c = list(c2 = mix_c2, c5 = mix_c5, c10 = mix_c10)
    )
  )
}

# ============================================================
# 11) Many repetitions + summaries
# ============================================================

summarize_detection <- function(times) {
  hit <- is.finite(times)
  data.frame(
    detect_pct = 100 * mean(hit),
    n_detect = sum(hit),
    mean_time = if (sum(hit) > 0) mean(times[hit]) else NA_real_,
    median_time = if (sum(hit) > 0) stats::median(times[hit]) else NA_real_
  )
}

run_many_tail <- function(R = 200,
                          seed0 = 1,
                          test = c("ES", "TailSD", "RVaR", "MedianShortfall"),
                          mode = c("oracle_true", "oracle_empirical", "ml_t"),
                          f_set = f_set_default,
                          apply_to = "r_only",
                          delta = delta_global,
                          N_mc = 3000L,
                          N_cond = 3000L,
                          refit_every = 10L,
                          weights = NULL,
                          verbose = TRUE) {
  test <- match.arg(test)
  mode <- match.arg(mode)
  apply_to <- normalize_apply_to_tail(apply_to)

  f_set <- as.numeric(f_set)
  out <- vector("list", length(f_set))
  names(out) <- paste0("f_", f_set)

  for (k in seq_along(f_set)) {
    f <- f_set[k]

    t_bonf <- rep(Inf, R)
    t_mix  <- rep(Inf, R)

    t_bonf_c2  <- rep(Inf, R)
    t_bonf_c5  <- rep(Inf, R)
    t_bonf_c10 <- rep(Inf, R)
    t_mix_c2   <- rep(Inf, R)
    t_mix_c5   <- rep(Inf, R)
    t_mix_c10  <- rep(Inf, R)

    who_bonf <- vector("list", R)

    avg1 <- matrix(NA_real_, nrow = R, ncol = 7)
    colnames(avg1) <- c("VaR_a", "ES_a", "TailSD", "MedSF", "VaR_low", "VaR_high", "RVaR")

    for (r in seq_len(R)) {
      rep_out <- run_one_rep_tail(
        test = test,
        mode = mode,
        f = f,
        apply_to = apply_to,
        delta = delta,
        N_mc = N_mc,
        N_cond = N_cond,
        refit_every = refit_every,
        weights = weights,
        seed = seed0 + r
      )

      t_bonf[r] <- rep_out$results$bonferroni$t
      t_mix[r]  <- rep_out$results$mixture$t

      t_bonf_c2[r]  <- rep_out$results$bonf_c$c2$t
      t_bonf_c5[r]  <- rep_out$results$bonf_c$c5$t
      t_bonf_c10[r] <- rep_out$results$bonf_c$c10$t

      t_mix_c2[r]   <- rep_out$results$mix_c$c2$t
      t_mix_c5[r]   <- rep_out$results$mix_c$c5$t
      t_mix_c10[r]  <- rep_out$results$mix_c$c10$t

      who_bonf[[r]] <- rep_out$results$bonferroni$who

      af <- rep_out$avg_forecasts
      avg1[r, ] <- c(
        af$avg_VaR_a,
        af$avg_ES_a,
        af$avg_TailSD,
        af$avg_MedSF,
        af$avg_VaR_low,
        af$avg_VaR_high,
        af$avg_RVaR
      )
    }

    avg_tbl <- as.data.frame(as.list(col_means_or_na(avg1)))

    out[[k]] <- list(
      settings = list(
        R = R,
        test = test,
        mode = mode,
        f = f,
        apply_to = apply_to,
        delta = delta,
        N_mc = N_mc,
        N_cond = N_cond,
        refit_every = refit_every
      ),
      bonf = list(times = t_bonf, summary = summarize_detection(t_bonf), who = who_bonf),
      mix  = list(times = t_mix,  summary = summarize_detection(t_mix)),
      bonf_c = list(
        c2  = list(times = t_bonf_c2,  summary = summarize_detection(t_bonf_c2)),
        c5  = list(times = t_bonf_c5,  summary = summarize_detection(t_bonf_c5)),
        c10 = list(times = t_bonf_c10, summary = summarize_detection(t_bonf_c10))
      ),
      mix_c = list(
        c2  = list(times = t_mix_c2,  summary = summarize_detection(t_mix_c2)),
        c5  = list(times = t_mix_c5,  summary = summarize_detection(t_mix_c5)),
        c10 = list(times = t_mix_c10, summary = summarize_detection(t_mix_c10))
      ),
      avg_forecasts = avg_tbl
    )

    if (isTRUE(verbose)) {
      cat("\n=== TEST =", test, "| mode =", mode, "| f =", f, "| apply_to =", apply_to, "===\n")
      cat("Bonferroni (delta):\n"); print(out[[k]]$bonf$summary)
      cat("Mixture (delta):\n"); print(out[[k]]$mix$summary)
      cat("Avg forecasts:\n"); print(out[[k]]$avg_forecasts)
      cat("Bonferroni (c = 2):\n"); print(out[[k]]$bonf_c$c2$summary)
      cat("Bonferroni (c = 5):\n"); print(out[[k]]$bonf_c$c5$summary)
      cat("Bonferroni (c = 10):\n"); print(out[[k]]$bonf_c$c10$summary)
      cat("Mixture (c = 2):\n"); print(out[[k]]$mix_c$c2$summary)
      cat("Mixture (c = 5):\n"); print(out[[k]]$mix_c$c5$summary)
      cat("Mixture (c = 10):\n"); print(out[[k]]$mix_c$c10$summary)
    }
  }

  out
}

# ============================================================
# 12) Grid run
# ============================================================

default_tail_scenarios <- function() {
  list(
    list(name = "SIZE_f1_none",       apply_to = "none",   f_set = c(1.0)),
    list(name = "ALT_f09_f11_r_only", apply_to = "r_only", f_set = c(0.9, 1.1)),
    list(name = "ALT_f09_f11_both",   apply_to = "both",   f_set = c(0.9, 1.1))
  )
}

run_tail_thesis_grid <- function(R = 1000L,
                                 seed_base = 10000L,
                                 tests_grid = c("ES", "TailSD", "RVaR", "MedianShortfall"),
                                 modes_grid = c("oracle_true", "ml_t", "oracle_empirical"),
                                 scenarios = NULL,
                                 delta = delta_global,
                                 N_mc = 3000L,
                                 N_cond = 3000L,
                                 refit_every = 10L,
                                 weights = NULL,
                                 verbose = TRUE,
                                 save_path = NULL) {
  if (is.null(scenarios)) {
    scenarios <- default_tail_scenarios()
  }

  results_all <- list()
  run_id <- 0L

  for (test in tests_grid) {
    test <- match.arg(test, choices = c("ES", "TailSD", "RVaR", "MedianShortfall"))

    for (mode in modes_grid) {
      mode <- match.arg(mode, choices = c("oracle_true", "oracle_empirical", "ml_t"))

      for (sc in scenarios) {
        run_id <- run_id + 1L
        key <- paste(test, mode, sc$name, sep = "__")

        if (isTRUE(verbose)) {
          cat("\n\n========================================================\n")
          cat("RUN", run_id, ":", key, "\n")
          cat(
            "R =", R,
            "| test =", test,
            "| mode =", mode,
            "| apply_to =", sc$apply_to,
            "| f_set =", paste(sc$f_set, collapse = ", "),
            "\n"
          )
          cat("========================================================\n")
        }

        results_all[[key]] <- run_many_tail(
          R = R,
          seed0 = seed_base + 100000L * run_id,
          test = test,
          mode = mode,
          f_set = sc$f_set,
          apply_to = sc$apply_to,
          delta = delta,
          N_mc = N_mc,
          N_cond = N_cond,
          refit_every = refit_every,
          weights = weights,
          verbose = verbose
        )
      }
    }
  }

  if (!is.null(save_path)) {
    saveRDS(results_all, file = save_path)
    if (isTRUE(verbose)) {
      cat("\nSaved results to:", save_path, "\n")
    }
  }

  results_all
}
