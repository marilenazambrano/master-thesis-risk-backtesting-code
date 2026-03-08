###############################################################################
# World A simulation + forecasting + canonical e-values + GREM(S0)
#
# This script:
# - simulates World A (AR(1) + GARCH(1,1) with skewed t innovations),
# - builds forecasts under different modes (oracle_true, oracle_empirical,
#   ml_norm, ml_t),
# - optionally distorts forecasts by a multiplicative scaling factor f,
# - computes canonical e-values for several risk measures,
# - runs GREM(S0) with the GREE/GREL mixture,
# - repeats the experiment R times and summarizes detection performance.
#
# Performance note:
# - For the "new risks" (Entropic / Expectile / OCE / EVaR / Distortion_ESgrid),
#   Monte Carlo sampling is done inside the backtest loop t = 1, ..., T_backtest.
# - For large R (for example R = 1000), it is better to run one mode and one
#   risk at a time and keep N_mc moderate.
###############################################################################

suppressPackageStartupMessages({
  library(rugarch)
  library(fGarch)   # rsstd / qsstd
  library(dplyr)
})

set.seed(123)

alpha_VaR <- 0.975
alpha_ES  <- 0.975

df_true   <- 5
skew_true <- 1.5

T_train    <- 500
T_backtest <- 500
T_total    <- T_train + T_backtest


# ============================================================
# 1) WORLD A
#    AR(1) mean + GARCH(1,1) volatility + skewed t innovations
# ============================================================

simulate_world_A <- function(T_total, df = df_true, skew = skew_true, burnin = 1000L) {
  T_all <- as.integer(T_total + burnin)

  L   <- numeric(T_all)   # losses
  mu  <- numeric(T_all)   # conditional mean
  sig <- numeric(T_all)   # conditional sd
  Z   <- numeric(T_all)   # innovations

  # GARCH(1,1)
  g_omega <- 0.01
  g_alpha <- 0.10
  g_beta  <- 0.85

  # AR(1) mean
  mu0 <- -0.05
  phi <-  0.30

  sig[1] <- sqrt(g_omega / (1 - g_alpha - g_beta))
  mu[1]  <- mu0
  Z[1]   <- fGarch::rsstd(1, mean = 0, sd = 1, nu = df, xi = skew)
  L[1]   <- mu[1] + sig[1] * Z[1]

  for (t in 2:T_all) {
    mu[t]  <- mu0 + phi * L[t - 1]
    sig[t] <- sqrt(g_omega + g_alpha * sig[t - 1]^2 * Z[t - 1]^2 + g_beta * sig[t - 1]^2)
    Z[t]   <- fGarch::rsstd(1, mean = 0, sd = 1, nu = df, xi = skew)
    L[t]   <- mu[t] + sig[t] * Z[t]
  }

  # drop burn-in
  idx <- (burnin + 1L):T_all
  list(L = L[idx], mu = mu[idx], sig = sig[idx], Z = Z[idx])
}

get_backtest_slice <- function(sim, T_train, T_backtest) {
  idx_bt <- (T_train + 1L):(T_train + T_backtest)
  list(
    idx_bt   = idx_bt,
    L_back   = sim$L[idx_bt],
    mu_back  = sim$mu[idx_bt],
    sig_back = sim$sig[idx_bt]
  )
}


# ============================================================
# 2) ORACLE CONSTANTS 
#    Used for oracle_true VaR/ES constants under skew-t innovations
# ============================================================

.oracle_cache <- new.env(parent = emptyenv())

oracle_true_constants_sstd <- function(alpha, df, skew, N_mc = 200000L) {
  key <- paste0("a=", alpha, "|df=", df, "|sk=", skew, "|N=", N_mc)
  if (exists(key, envir = .oracle_cache, inherits = FALSE)) {
    return(get(key, envir = .oracle_cache, inherits = FALSE))
  }
  q <- as.numeric(fGarch::qsstd(alpha, mean = 0, sd = 1, nu = df, xi = skew))
  Z <- fGarch::rsstd(as.integer(N_mc), mean = 0, sd = 1, nu = df, xi = skew)
  tailZ <- Z[Z > q]
  es <- if (length(tailZ) == 0L) q else mean(tailZ)
  out <- list(q = q, es = es, N_mc = N_mc)
  assign(key, out, envir = .oracle_cache)
  out
}

# Gaussian constants 
const_norm <- function(alpha) {
  q <- qnorm(alpha)
  es <- dnorm(q) / (1 - alpha)
  list(q = q, es = es)
}

# Standardized t with var = 1: Z = s * T, s = sqrt((nu - 2) / nu)
const_std <- function(alpha, nu) {
  nu <- as.numeric(nu)
  if (!is.finite(nu) || nu <= 2) {
    # fallback
    nu2 <- max(nu, 2.1)
    s <- sqrt((nu2 - 2) / nu2)
    qT <- qt(alpha, df = nu2)
    esT <- (dt(qT, df = nu2) / (1 - alpha)) * ((nu2 + qT^2) / (nu2 - 1))
    return(list(q = s * qT, es = s * esT, nu = nu2, note = "nu<=2 fallback"))
  }
  s <- sqrt((nu - 2) / nu)
  qT <- qt(alpha, df = nu)
  esT <- (dt(qT, df = nu) / (1 - alpha)) * ((nu + qT^2) / (nu - 1))
  list(q = s * qT, es = s * esT, nu = nu)
}


# ============================================================
# 3) ML FITTING (rolling refits)
#    Fits AR(1)-GARCH(1,1) with dist in {norm, std}
# ============================================================

fit_ar1_garch11 <- function(x, dist = c("norm", "std")) {
  dist <- match.arg(dist)
  spec <- rugarch::ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model     = list(armaOrder = c(1, 0), include.mean = TRUE),
    distribution.model = dist
  )
  rugarch::ugarchfit(spec = spec, data = x, solver = "hybrid")
}

# Build mu_{t|t-1} and sigma_{t|t-1} forecasts for backtest times under ML
# rolling refits. Also returns nu_used[t] for dist = "std".
build_ml_mu_sig_worldA <- function(L_full, T_train, T_backtest,
                                   dist = c("norm", "std"),
                                   refit_every = 50L) {
  dist <- match.arg(dist)
  refit_every <- as.integer(refit_every)
  if (refit_every <= 0L) stop("refit_every must be positive.")

  mu_f  <- rep(NA_real_, T_backtest)
  sig_f <- rep(NA_real_, T_backtest)
  nu_f  <- rep(NA_real_, T_backtest)

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
    if (is.null(fit)) break

    cf <- coef(fit)
    mu0 <- if ("mu" %in% names(cf)) as.numeric(cf[["mu"]]) else 0
    ar1_name <- intersect(names(cf), c("ar1", "ar[1]", "ar"))
    ar1 <- if (length(ar1_name)) as.numeric(cf[[ar1_name[1]]]) else 0
    omg <- if ("omega"  %in% names(cf)) as.numeric(cf[["omega"]])  else NA_real_
    a1  <- if ("alpha1" %in% names(cf)) as.numeric(cf[["alpha1"]]) else NA_real_
    b1  <- if ("beta1"  %in% names(cf)) as.numeric(cf[["beta1"]])  else NA_real_

    nu_block <- NA_real_
    if (dist == "std" && ("shape" %in% names(cf))) nu_block <- as.numeric(cf[["shape"]])

    sig_series <- as.numeric(rugarch::sigma(fit))
    mu_series  <- as.numeric(rugarch::fitted(fit))

    L_prev   <- L_full[t_prev_full]
    sig_prev <- sig_series[length(sig_series)]
    mu_prev  <- mu_series[length(mu_series)]
    res_prev <- L_prev - mu_prev

    block_len <- min(refit_every, T_backtest - k + 1L)

    for (j in 1:block_len) {
      t_full <- t_prev_full + 1L

      mu1  <- mu0 + ar1 * L_prev
      sig1 <- sqrt(pmax(omg + a1 * (res_prev^2) + b1 * (sig_prev^2), 0))

      mu_f[k]  <- mu1
      sig_f[k] <- sig1
      nu_f[k]  <- nu_block

      # update with realized observation
      L_prev   <- L_full[t_full]
      mu_prev  <- mu1
      sig_prev <- sig1
      res_prev <- L_prev - mu_prev

      t_prev_full <- t_full
      k <- k + 1L
      if (k > T_backtest) break
    }
  }

  list(mu = mu_f, sig = sig_f, nu = nu_f, dist = dist, refit_every = refit_every)
}

# Oracle empirical: presample standardized residuals from truth
# (simulation knows mu and sig)
estimate_Zhat_presample <- function(sim, presample_length) {
  idx <- 1:as.integer(presample_length)
  Z_hat <- (sim$L[idx] - sim$mu[idx]) / sim$sig[idx]
  Z_hat[is.finite(Z_hat)]
}


# ============================================================
# 4) RISK-FORECASTS
#    For VaR/ES: analytic constants
#    For other risks: MC predictive sampling per t
# ============================================================

# Draw innovations depending on forecast mode
draw_Z <- function(N, mode, Z_hat = NULL, nu = NA_real_) {
  N <- as.integer(N)
  if (N <= 0L) stop("N must be positive.")

  if (mode == "oracle_true") {
    return(fGarch::rsstd(N, mean = 0, sd = 1, nu = df_true, xi = skew_true))
  }
  if (mode == "oracle_empirical") {
    if (is.null(Z_hat) || length(Z_hat) < 10) stop("Need Z_hat for oracle_empirical.")
    return(sample(Z_hat, size = N, replace = TRUE))
  }
  if (mode == "ml_norm") {
    return(rnorm(N))
  }
  if (mode == "ml_t") {
    nu2 <- as.numeric(nu)
    if (!is.finite(nu2) || nu2 <= 2) nu2 <- 8
    s <- sqrt((nu2 - 2) / nu2)
    return(s * rt(N, df = nu2))
  }
  stop("Unknown mode in draw_Z(): ", mode)
}

# Expectile identification function V(r, x)
V_expectile <- function(r, x, tau) {
  (1 - tau) * (r - x) * (x <= r) + tau * (r - x) * (x > r)
}

# Sample expectile: solve mean[V(r, X)] = 0
expectile_sample <- function(X, tau) {
  X <- as.numeric(X)
  X <- X[is.finite(X)]
  if (length(X) < 5) return(NA_real_)

  f <- function(r) mean(V_expectile(r, X, tau))
  lo <- min(X)
  hi <- max(X)

  flo <- f(lo)
  fhi <- f(hi)
  if (!is.finite(flo) || !is.finite(fhi)) return(NA_real_)
  if (flo == 0) return(lo)
  if (fhi == 0) return(hi)

  # If no sign change, expand interval; otherwise fallback to median
  if (flo * fhi > 0) {
    span <- hi - lo
    for (k in 1:6) {
      lo2 <- lo - 0.5 * span
      hi2 <- hi + 0.5 * span
      flo <- f(lo2)
      fhi <- f(hi2)
      if (is.finite(flo) && is.finite(fhi) && flo * fhi <= 0) {
        lo <- lo2
        hi <- hi2
        break
      }
      span <- 2 * span
    }
    if (flo * fhi > 0) return(stats::median(X))
  }

  uniroot(f, lower = lo, upper = hi, tol = 1e-10)$root
}

phi_p_shortfall <- function(t, p) (pmax(t, 0)^p) / p

phi_huber <- function(t, c) {
  ifelse(t <= 0, 0,
         ifelse(t <= c, 0.5 * t^2,
                c * t - 0.5 * c^2))
}

# OCE: optimize m -> m + mean(phi(X - m))
oce_fit_from_sample <- function(X, phi_fun, m_lo = NULL, m_hi = NULL) {
  X <- as.numeric(X)
  X <- X[is.finite(X)]
  if (length(X) < 10) return(list(z = NA_real_, r = NA_real_))

  if (is.null(m_lo)) m_lo <- stats::quantile(X, 0.01, type = 8) - stats::sd(X)
  if (is.null(m_hi)) m_hi <- stats::quantile(X, 0.99, type = 8) + stats::sd(X)

  obj <- function(m) m + mean(phi_fun(X - m))
  z <- optimize(obj, interval = c(m_lo, m_hi), tol = 1e-8)$minimum
  r <- z + mean(phi_fun(X - z))
  list(z = z, r = r)
}

# EVaR: minimize over t > 0:
# (1 / t) * (log(mean(exp(tX))) - log(1 - alpha))
evar_fit_from_sample <- function(X, alpha, t_lo = 1e-4, t_hi = 50) {
  X <- as.numeric(X)
  X <- X[is.finite(X)]
  if (length(X) < 10) return(list(t_star = NA_real_, r = NA_real_))

  log1ma <- log(1 - alpha)

  # stable log-mean-exp
  lme <- function(v) {
    m <- max(v)
    m + log(mean(exp(v - m)))
  }

  obj <- function(t) {
    t <- as.numeric(t)
    if (!is.finite(t) || t <= 0) return(Inf)
    val <- (lme(t * X) - log1ma) / t
    if (!is.finite(val)) Inf else val
  }

  opt <- optimize(obj, interval = c(t_lo, t_hi), tol = 1e-6)
  t_star <- opt$minimum
  r <- opt$objective
  list(t_star = t_star, r = r)
}

# Distortion ES-grid: compute VaR_u and ES_u from sample X
es_from_sample <- function(X, u) {
  q <- as.numeric(stats::quantile(X, probs = u, type = 8))
  tail <- X[X > q]
  es <- if (length(tail) == 0) q else mean(tail)
  list(q = q, es = es)
}

# Build forecast arrays for a given risk type
build_risk_forecasts_worldA <- function(sim,
                                        forecast_mode = c("oracle_empirical", "oracle_true", "ml_norm", "ml_t"),
                                        risk_type = c("VaR", "ES", "Entropic", "Expectile", "OCE_pShortfall", "OCE_Huber", "EVaR", "Distortion_ESgrid"),
                                        T_train, T_backtest,
                                        alpha_VaR, alpha_ES,
                                        N_mc = 20000L,
                                        a_ent = 4,
                                        tau_ex = 0.975,
                                        a_lower = -10,
                                        p_short = 2,
                                        huber_c = 1.5,
                                        alpha_evar = 0.975,
                                        u_vec = c(0.90, 0.95, 0.975),
                                        w_vec = c(0.2, 0.3, 0.5),
                                        refit_every = 50L,
                                        oracle_mc_N = 200000L,
                                        q_type = 8) {

  forecast_mode <- match.arg(forecast_mode)
  risk_type <- match.arg(risk_type)

  # backtest slice indices
  sl <- get_backtest_slice(sim, T_train, T_backtest)
  idx_bt <- sl$idx_bt

  # Build mu_{t|t-1}, sig_{t|t-1} and possibly nu_t under the chosen mode
  if (forecast_mode %in% c("oracle_true", "oracle_empirical")) {
    mu_f  <- sim$mu[idx_bt]
    sig_f <- sim$sig[idx_bt]
    nu_f  <- rep(NA_real_, T_backtest)
  } else if (forecast_mode == "ml_norm") {
    st <- build_ml_mu_sig_worldA(sim$L, T_train, T_backtest, dist = "norm", refit_every = refit_every)
    mu_f <- st$mu
    sig_f <- st$sig
    nu_f <- st$nu
  } else {
    st <- build_ml_mu_sig_worldA(sim$L, T_train, T_backtest, dist = "std", refit_every = refit_every)
    mu_f <- st$mu
    sig_f <- st$sig
    nu_f <- st$nu
  }

  if (any(!is.finite(mu_f)) || any(!is.finite(sig_f))) {
    stop("Forecast state (mu_f/sig_f) contains NA/Inf. Likely ML fit failed. Investigate refit_every or data.")
  }

  # Oracle empirical Z_hat 
  Z_hat <- NULL
  if (forecast_mode == "oracle_empirical") {
    Z_hat <- estimate_Zhat_presample(sim, presample_length = T_train)
  }

  # ---- VaR / ES forecasts are analytic ----
  if (risk_type %in% c("VaR", "ES")) {
    VaR <- numeric(T_backtest)
    ES  <- numeric(T_backtest)

    if (forecast_mode == "oracle_true") {
      qV <- as.numeric(fGarch::qsstd(alpha_VaR, mean = 0, sd = 1, nu = df_true, xi = skew_true))
      constES <- oracle_true_constants_sstd(alpha_ES, df_true, skew_true, N_mc = oracle_mc_N)
      esE <- constES$es
      VaR <- mu_f + sig_f * qV
      ES  <- mu_f + sig_f * esE

    } else if (forecast_mode == "oracle_empirical") {
      q_hat <- as.numeric(stats::quantile(Z_hat, probs = alpha_VaR, type = q_type))
      q_es  <- as.numeric(stats::quantile(Z_hat, probs = alpha_ES, type = q_type))
      tailZ <- Z_hat[Z_hat > q_es]
      es_hat <- if (length(tailZ) == 0) q_es else mean(tailZ)

      VaR <- mu_f + sig_f * q_hat
      ES  <- mu_f + sig_f * es_hat

    } else if (forecast_mode == "ml_norm") {
      cV <- const_norm(alpha_VaR)
      cE <- const_norm(alpha_ES)
      VaR <- mu_f + sig_f * cV$q
      ES  <- mu_f + sig_f * cE$es

    } else {  # ml_t
      VaR <- numeric(T_backtest)
      ES  <- numeric(T_backtest)
      for (t in 1:T_backtest) {
        nu_t <- nu_f[t]
        cV <- const_std(alpha_VaR, nu_t)
        cE <- const_std(alpha_ES, nu_t)
        VaR[t] <- mu_f[t] + sig_f[t] * cV$q
        ES[t]  <- mu_f[t] + sig_f[t] * cE$es
      }
    }

    return(list(
      type = risk_type,
      mode = forecast_mode,
      mu = mu_f,
      sig = sig_f,
      nu = nu_f,
      VaR = VaR,
      ES  = ES
    ))
  }

  # ---- For new risks: MC predictive sampling at each t ----
  r  <- rep(NA_real_, T_backtest)
  z  <- rep(NA_real_, T_backtest)     # for OCE
  t_star <- rep(NA_real_, T_backtest) # for EVaR
  z_mat <- NULL                       # for distortion ES-grid

  if (risk_type == "Distortion_ESgrid") {
    u_vec <- as.numeric(u_vec)
    w_vec <- as.numeric(w_vec)
    if (length(u_vec) != length(w_vec)) stop("u_vec and w_vec must have same length.")
    if (any(u_vec <= 0 | u_vec >= 1)) stop("u_vec must be in (0,1).")
    sw <- sum(w_vec)
    if (!is.finite(sw) || sw <= 0) stop("w_vec must sum to positive.")
    w_vec <- w_vec / sw
    z_mat <- matrix(NA_real_, nrow = T_backtest, ncol = length(u_vec))
    colnames(z_mat) <- paste0("u", u_vec)
  }

  for (t in 1:T_backtest) {
    Z <- draw_Z(N = N_mc, mode = forecast_mode, Z_hat = Z_hat, nu = nu_f[t])
    X <- mu_f[t] + sig_f[t] * Z

    if (risk_type == "Entropic") {
      a <- as.numeric(a_ent)
      if (!is.finite(a) || a <= 0) stop("a_ent must be > 0.")
      # r_t = (1 / a) log E[exp(aX)]
      m <- max(a * X)
      r[t] <- (m + log(mean(exp(a * X - m)))) / a

    } else if (risk_type == "Expectile") {
      r[t] <- expectile_sample(X, tau = tau_ex)

    } else if (risk_type == "OCE_pShortfall") {
      fit <- oce_fit_from_sample(X, phi_fun = function(u) phi_p_shortfall(u, p = p_short))
      z[t] <- fit$z
      r[t] <- fit$r

    } else if (risk_type == "OCE_Huber") {
      fit <- oce_fit_from_sample(X, phi_fun = function(u) phi_huber(u, c = huber_c))
      z[t] <- fit$z
      r[t] <- fit$r

    } else if (risk_type == "EVaR") {
      fit <- evar_fit_from_sample(X, alpha = alpha_evar, t_lo = 1e-4, t_hi = 50)
      t_star[t] <- fit$t_star
      r[t] <- fit$r

    } else if (risk_type == "Distortion_ESgrid") {
      es_vec <- numeric(length(u_vec))
      for (j in seq_along(u_vec)) {
        tmp <- es_from_sample(X, u = u_vec[j])
        z_mat[t, j] <- tmp$q
        es_vec[j] <- tmp$es
      }
      r[t] <- sum(w_vec * es_vec)

    } else {
      stop("Unknown risk_type in MC part: ", risk_type)
    }
  }

  out <- list(
    type = risk_type,
    mode = forecast_mode,
    mu = mu_f,
    sig = sig_f,
    nu = nu_f,
    r = r
  )
  if (risk_type %in% c("OCE_pShortfall", "OCE_Huber")) out$z <- z
  if (risk_type == "EVaR") out$t_star <- t_star
  if (risk_type == "Distortion_ESgrid") {
    out$u_vec <- as.numeric(u_vec)
    out$w_vec <- as.numeric(w_vec)
    out$z_mat <- z_mat
  }
  out
}


# ============================================================
# 5) DISTORTION (simple scaling; apply_to = "all" or "r_only")
# ============================================================

scale_abs <- function(x, f) x + (f - 1) * abs(x)   # f = 0.9 down, f = 1.1 up

distort_forecast_scale <- function(forecast_obj, f = 1.0,
                                   apply_to = c("all", "r_only")) {
  apply_to <- match.arg(apply_to)
  stopifnot(is.finite(f), f > 0)

  out <- forecast_obj
  out$dist_f <- f
  out$dist_apply_to <- apply_to

  type <- out$type

  if (type == "VaR") {
    out$VaR <- scale_abs(out$VaR, f)
    return(out)
  }

  if (type == "ES") {
    out$ES <- scale_abs(out$ES, f)
    if (apply_to == "all") out$VaR <- scale_abs(out$VaR, f)
    return(out)
  }

  if (type %in% c("Entropic", "Expectile", "EVaR")) {
    out$r <- scale_abs(out$r, f)
    return(out)
  }
  if (type %in% c("OCE_pShortfall", "OCE_Huber")) {
    out$r <- scale_abs(out$r, f)
    if (apply_to == "all") out$z <- scale_abs(out$z, f)
    return(out)
  }

  if (type == "Distortion_ESgrid") {
    out$r <- scale_abs(out$r, f)
    if (apply_to == "all") out$z_mat <- scale_abs(out$z_mat, f)
    return(out)
  }

  stop("Unknown risk type in distort_forecast_scale(): ", type)
}


# ============================================================
# 6) E-VALUES 
# ============================================================

# VaR canonical e-value (indicator exceedance)
e_var <- function(L, VaR, alpha) (1 / (1 - alpha)) * as.numeric(L > VaR)

# ES canonical e-value (requires ES > VaR)
e_es <- function(L, ES, VaR, alpha) {
  den <- (1 - alpha) * (ES - VaR)
  if (!is.finite(den) || den <= 0) return(Inf)
  pmax(L - VaR, 0) / den
}

# Expectile canonical e-value:
# e = 1 - V(r, L) / ((1 - tau) (r - a_lower))
e_expectile <- function(L, r, tau, a_lower) {
  den <- (1 - tau) * (r - a_lower)
  if (!is.finite(den) || den <= 0) return(Inf)
  1 - V_expectile(r, L, tau) / den
}

# Entropic canonical e-value:
# e = exp(a(L - r))
e_entropic <- function(L, r, a_ent) {
  a <- as.numeric(a_ent)
  if (!is.finite(a) || a <= 0) return(Inf)
  val <- a * (L - r)
  if (val > 700) return(Inf)  # avoid overflow
  exp(val)
}

# OCE e-value (phi >= 0 version):
# e = phi(L - z) / (r - z)
e_oce <- function(L, r, z, phi_fun) {
  den <- r - z
  if (!is.finite(den) || den <= 0) return(Inf)
  num <- phi_fun(L - z)
  if (!is.finite(num) || num < 0) num <- Inf
  num / den
}

# EVaR e-value:
# e = exp(t * (L - r)) / (1 - alpha)
e_evar <- function(L, r, t_star, alpha) {
  if (!is.finite(t_star) || t_star <= 0) return(Inf)
  if (!is.finite(r)) return(Inf)
  val <- t_star * (L - r)
  if (val > 700) return(Inf)
  exp(val) / (1 - alpha)
}

# Distortion ES-grid e-value:
# e = [sum w_i (L - z_i)_+ / (1 - u_i)] / [r - sum w_i z_i]
e_dist_esgrid <- function(L, r, z_vec, u_vec, w_vec) {
  u_vec <- as.numeric(u_vec)
  w_vec <- as.numeric(w_vec)
  z_vec <- as.numeric(z_vec)
  if (length(u_vec) != length(w_vec) || length(z_vec) != length(u_vec)) return(Inf)

  den <- r - sum(w_vec * z_vec)
  if (!is.finite(den) || den <= 0) return(Inf)

  num <- sum(w_vec * (pmax(L - z_vec, 0) / (1 - u_vec)))
  if (!is.finite(num) || num < 0) return(Inf)

  num / den
}

# Dispatch realized e-value at time t (uses global alpha_VaR / alpha_ES)
e_at <- function(L, t, fc,
                 tau_ex, a_lower, a_ent, p_short, huber_c, alpha_evar) {
  type <- fc$type
  if (type == "VaR") return(e_var(L, fc$VaR[t], alpha = alpha_VaR))
  if (type == "ES")  return(e_es(L, fc$ES[t], fc$VaR[t], alpha = alpha_ES))
  if (type == "Expectile") return(e_expectile(L, fc$r[t], tau = tau_ex, a_lower = a_lower))
  if (type == "Entropic")  return(e_entropic(L, fc$r[t], a_ent = a_ent))
  if (type == "OCE_pShortfall") {
    return(e_oce(L, fc$r[t], fc$z[t], phi_fun = function(u) phi_p_shortfall(u, p = p_short)))
  }
  if (type == "OCE_Huber") {
    return(e_oce(L, fc$r[t], fc$z[t], phi_fun = function(u) phi_huber(u, c = huber_c)))
  }
  if (type == "EVaR") {
    return(e_evar(L, fc$r[t], fc$t_star[t], alpha = alpha_evar))
  }
  if (type == "Distortion_ESgrid") {
    return(e_dist_esgrid(L, fc$r[t], fc$z_mat[t, ], fc$u_vec, fc$w_vec))
  }
  stop("Unknown fc$type in e_at(): ", type)
}

# Optional helper: recompute e(L_i, forecast_at_t) for a set of past losses
e_re_at <- function(L_vec, t, fc,
                    tau_ex, a_lower, a_ent, p_short, huber_c, alpha_evar) {
  sapply(L_vec, function(Li) e_at(Li, t = t, fc = fc,
                                  tau_ex = tau_ex, a_lower = a_lower,
                                  a_ent = a_ent, p_short = p_short, huber_c = huber_c,
                                  alpha_evar = alpha_evar))
}


# ============================================================
# 7) GREE / GREL lambdas on S0 + GREM mixture
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

weight_from_logs <- function(logA, logB) {
  if (is.infinite(logA) && logA < 0 && is.infinite(logB) && logB < 0) return(0.5)
  if (is.infinite(logA) && logA < 0) return(0)
  if (is.infinite(logB) && logB < 0) return(1)
  1 / (1 + exp(logB - logA))
}

run_grem_S0 <- function(L_back, fc,
                        J_lambda = 500L,
                        lambda_max = 0.99,
                        tau_ex = 0.975,
                        a_lower = -10,
                        a_ent = 4,
                        p_short = 2,
                        huber_c = 1.5,
                        alpha_evar = 0.975) {

  T <- length(L_back)
  J <- as.integer(J_lambda)
  if (J <= 0) stop("J_lambda must be positive.")

  # realized e_t series under the time-t forecast
  e_series <- numeric(T)
  for (t in 1:T) {
    e_series[t] <- e_at(L_back[t], t, fc,
                        tau_ex = tau_ex, a_lower = a_lower,
                        a_ent = a_ent, p_short = p_short, huber_c = huber_c,
                        alpha_evar = alpha_evar)
  }

  # log-wealth processes
  M_E   <- numeric(T + 1L)
  M_E[1] <- 0
  M_L   <- numeric(T + 1L)
  M_L[1] <- 0
  M_mix <- numeric(T + 1L)
  M_mix[1] <- 0

  lam_E <- numeric(T)
  lam_L <- numeric(T)
  lam_M <- numeric(T)

  for (t in 1:T) {
    past_end <- t - 1L
    if (past_end <= 0L) {
      lam_E[t] <- 0
      lam_L[t] <- 0
    } else {
      idx0 <- max(1L, past_end - J + 1L):past_end

      # GREE: lambda from past realized e-values
      lam_E[t] <- lambda_from_e_taylor(e_series[idx0], lambda_max = lambda_max)

      # GREL: recompute e-values of past losses under CURRENT forecast at time t
      e_re <- sapply(idx0, function(i) e_at(L_back[i], t, fc,
                                            tau_ex = tau_ex, a_lower = a_lower,
                                            a_ent = a_ent, p_short = p_short, huber_c = huber_c,
                                            alpha_evar = alpha_evar))
      lam_L[t] <- lambda_from_e_taylor(e_re, lambda_max = lambda_max)
    }

    # update base wealths with realized e_t
    M_E[t + 1L] <- update_log_wealth(M_E[t], e_series[t], lam_E[t])
    M_L[t + 1L] <- update_log_wealth(M_L[t], e_series[t], lam_L[t])

    # GREM mixture
    wE <- weight_from_logs(M_E[t], M_L[t])
    lam_M[t] <- wE * lam_E[t] + (1 - wE) * lam_L[t]
    M_mix[t + 1L] <- update_log_wealth(M_mix[t], e_series[t], lam_M[t])
  }

  list(
    M = M_mix,
    lambda = lam_M,
    base = list(M_E = M_E, M_L = M_L, lambda_E = lam_E, lambda_L = lam_L),
    e = e_series
  )
}


# ============================================================
# 8) ONE REPLICATION
# ============================================================

run_one_rep_A <- function(risk_type = c("VaR", "ES", "Entropic", "Expectile", "OCE_pShortfall", "OCE_Huber", "EVaR", "Distortion_ESgrid"),
                          forecast_mode = c("oracle_empirical", "oracle_true", "ml_norm", "ml_t"),
                          f = 1.0,
                          dist_apply_to = c("all", "r_only"),
                          J_lambda = 500L,
                          lambda_max = 0.99,
                          N_mc = 20000L,
                          a_ent = 4,
                          tau_ex = 0.975,
                          a_lower = -10,
                          p_short = 2,
                          huber_c = 1.5,
                          alpha_evar = 0.975,
                          u_vec = c(0.90, 0.95, 0.975),
                          w_vec = c(0.2, 0.3, 0.5),
                          refit_every = 50L,
                          oracle_mc_N = 200000L,
                          q_type = 8) {

  risk_type <- match.arg(risk_type)
  forecast_mode <- match.arg(forecast_mode)
  dist_apply_to <- match.arg(dist_apply_to)

  sim <- simulate_world_A(T_total)
  sl <- get_backtest_slice(sim, T_train, T_backtest)
  L_back <- sl$L_back

  fc <- build_risk_forecasts_worldA(
    sim = sim,
    forecast_mode = forecast_mode,
    risk_type = risk_type,
    T_train = T_train,
    T_backtest = T_backtest,
    alpha_VaR = alpha_VaR,
    alpha_ES  = alpha_ES,
    N_mc = N_mc,
    a_ent = a_ent,
    tau_ex = tau_ex,
    a_lower = a_lower,
    p_short = p_short,
    huber_c = huber_c,
    alpha_evar = alpha_evar,
    u_vec = u_vec,
    w_vec = w_vec,
    refit_every = refit_every,
    oracle_mc_N = oracle_mc_N,
    q_type = q_type
  )

  # apply distortion scaling
  fc <- distort_forecast_scale(fc, f = f, apply_to = dist_apply_to)

  # run S0-only GREM
  grem <- run_grem_S0(
    L_back = L_back, fc = fc,
    J_lambda = J_lambda,
    lambda_max = lambda_max,
    tau_ex = tau_ex,
    a_lower = a_lower,
    a_ent = a_ent,
    p_short = p_short,
    huber_c = huber_c,
    alpha_evar = alpha_evar
  )

  list(
    GREM = grem,
    series = list(L = L_back, fc = fc),
    meta = list(
      world = "A",
      risk_type = risk_type,
      forecast_mode = forecast_mode,
      f = f,
      dist_apply_to = dist_apply_to,
      J_lambda = J_lambda,
      lambda_max = lambda_max,
      N_mc = N_mc
    )
  )
}


# ============================================================
# 9) MANY REPS + DETECTION SUMMARY
# ============================================================

detect_time <- function(M, c = 5) {
  thr <- log(c)
  hit <- which(M >= thr)
  if (length(hit) == 0) return(Inf)
  hit[1] - 1
}

summarize_detection <- function(det_store, T_max = T_backtest) {
  thr_names <- colnames(det_store)   # "c2", "c5", "c10"
  thr_vals  <- as.numeric(sub("^c", "", thr_names))

  res <- lapply(seq_along(thr_names), function(j) {
    v <- det_store[, j]
    hit <- is.finite(v)
    n_hit <- sum(hit)

    v_trunc <- pmin(v, T_max)   # truncate non-detections at horizon
    arl <- mean(v_trunc)

    data.frame(
      threshold   = thr_vals[j],
      detect_pct  = 100 * mean(hit),
      arl         = arl,
      n_detect    = n_hit,
      mean_time   = if (n_hit > 0) mean(v[hit]) else NA_real_,
      median_time = if (n_hit > 0) median(v[hit]) else NA_real_
    )
  })

  do.call(rbind, res)
}

print_det <- function(det_df, title = "") {
  cat("\n====================\n")
  cat(title, "\n")
  cat("====================\n")
  print(det_df)
}

run_many_A <- function(R = 1000,
                       seed0 = 1,
                       risk_type = "ES",
                       forecast_mode = "ml_t",
                       f_set = c(0.9, 1.0, 1.1),
                       dist_apply_to = "r_only",
                       J_lambda = 500L,
                       lambda_max = 0.99,
                       N_mc = 20000L,
                       a_ent = 4,
                       tau_ex = 0.975,
                       a_lower = -10,
                       p_short = 2,
                       huber_c = 1.5,
                       alpha_evar = 0.975,
                       u_vec = c(0.90, 0.95, 0.975),
                       w_vec = c(0.2, 0.3, 0.5),
                       refit_every = 50L,
                       oracle_mc_N = 200000L,
                       q_type = 8) {

  f_set <- as.numeric(f_set)
  out <- vector("list", length(f_set))
  names(out) <- paste0("f_", f_set)

  for (k in seq_along(f_set)) {
    f <- f_set[k]

    M_store <- matrix(NA_real_, nrow = R, ncol = T_backtest + 1L)
    det_store <- matrix(NA_real_, nrow = R, ncol = 3,
                        dimnames = list(NULL, c("c2", "c5", "c10")))

    # To compute average forecasts across replications
    sum_VaR <- NULL
    sum_ES  <- NULL
    sum_r   <- NULL
    sum_z   <- NULL
    sum_t_star <- NULL
    sum_z_mat  <- NULL
    u_keep <- NULL
    w_keep <- NULL
    fc_type_keep <- NULL

    for (r in 1:R) {
      set.seed(seed0 + r)

      rep_out <- run_one_rep_A(
        risk_type = risk_type,
        forecast_mode = forecast_mode,
        f = f,
        dist_apply_to = dist_apply_to,
        J_lambda = J_lambda,
        lambda_max = lambda_max,
        N_mc = N_mc,
        a_ent = a_ent,
        tau_ex = tau_ex,
        a_lower = a_lower,
        p_short = p_short,
        huber_c = huber_c,
        alpha_evar = alpha_evar,
        u_vec = u_vec,
        w_vec = w_vec,
        refit_every = refit_every,
        oracle_mc_N = oracle_mc_N,
        q_type = q_type
      )
      fc <- rep_out$series$fc

      # initialize sums on first replication
      if (r == 1) {
        fc_type_keep <- fc$type

        if (fc$type %in% c("VaR", "ES")) {
          sum_VaR <- rep(0, length(fc$VaR))
          sum_ES  <- rep(0, length(fc$ES))

        } else if (fc$type == "Distortion_ESgrid") {
          sum_r <- rep(0, length(fc$r))
          sum_z_mat <- 0 * fc$z_mat
          u_keep <- fc$u_vec
          w_keep <- fc$w_vec

        } else {
          sum_r <- rep(0, length(fc$r))
          if (!is.null(fc$z))      sum_z      <- rep(0, length(fc$z))
          if (!is.null(fc$t_star)) sum_t_star <- rep(0, length(fc$t_star))
        }
      }

      # add to sums
      if (fc$type %in% c("VaR", "ES")) {
        sum_VaR <- sum_VaR + fc$VaR
        sum_ES  <- sum_ES  + fc$ES

      } else if (fc$type == "Distortion_ESgrid") {
        sum_r     <- sum_r + fc$r
        sum_z_mat <- sum_z_mat + fc$z_mat

      } else {
        sum_r <- sum_r + fc$r
        if (!is.null(fc$z))      sum_z      <- sum_z + fc$z
        if (!is.null(fc$t_star)) sum_t_star <- sum_t_star + fc$t_star
      }

      # store paths and detection times
      M <- rep_out$GREM$M
      M_store[r, ] <- M
      det_store[r, "c2"]  <- detect_time(M, c = 2)
      det_store[r, "c5"]  <- detect_time(M, c = 5)
      det_store[r, "c10"] <- detect_time(M, c = 10)
    }

    # average forecasts across replications (optional)
    forecast_avg <- NULL

    if (!is.null(sum_ES)) {
      forecast_avg <- list(
        type = "VaR/ES",
        VaR_avg = sum_VaR / R,
        ES_avg  = sum_ES  / R
      )

    } else if (!is.null(sum_z_mat)) {
      forecast_avg <- list(
        type = "Distortion_ESgrid",
        r_avg = sum_r / R,
        z_mat_avg = sum_z_mat / R,
        u_vec = u_keep,
        w_vec = w_keep
      )

    } else if (!is.null(sum_r)) {
      forecast_avg <- list(
        type = fc_type_keep,
        r_avg = sum_r / R
      )
      if (!is.null(sum_z))      forecast_avg$z_avg      <- sum_z / R
      if (!is.null(sum_t_star)) forecast_avg$t_star_avg <- sum_t_star / R
    }

    det_summary <- summarize_detection(det_store)
    print(det_summary)

    out[[k]] <- list(
      M_store = M_store,
      det_store = det_store,
      det_summary = det_summary,
      forecast_avg = forecast_avg,
      settings = list(
        risk_type = risk_type,
        forecast_mode = forecast_mode,
        f = f,
        dist_apply_to = dist_apply_to,
        J_lambda = J_lambda,
        lambda_max = lambda_max,
        N_mc = N_mc
      )
    )
  }

  out
}


# ============================================================
# 10) GRID RUN WRAPPER (World A)
#     Loops over risk_types × apply_to × forecast_modes × f_set
# ============================================================

run_grid_worldA_allrisks <- function(
    R = 1000,
    seed0 = 123,
    risk_types = c("VaR", "ES", "Entropic", "Expectile", "OCE_pShortfall", "OCE_Huber", "EVaR", "Distortion_ESgrid"),
    forecast_modes = c("oracle_true", "oracle_empirical", "ml_t", "ml_norm"),
    f_set = c(0.9, 1.0, 1.1),
    dist_apply_default = "r_only",  # fallback if risk not in apply_map
    N_mc = 2000,
    # risk params
    a_ent = 4,
    tau_ex = 0.975,
    a_lower = -10,
    p_short = 2,
    huber_c = 1.5,
    alpha_evar = 0.975,
    u_vec = c(0.90, 0.95, 0.975),
    w_vec = c(0.2, 0.3, 0.5),
    oracle_mc_N = 200000L,
    q_type = 8,
    verbose_print = TRUE
) {

  # which apply_to variants to run per risk type
  apply_map <- list(
    VaR = c("all"),                       # scaling VaR is always the same notion
    ES  = c("r_only", "all"),            # ES-only scaling vs scaling both VaR and ES
    Entropic = c("r_only"),
    Expectile = c("r_only"),
    EVaR = c("r_only"),
    OCE_pShortfall = c("r_only", "all"), # scale r only vs scale (r, z)
    OCE_Huber = c("r_only", "all"),
    Distortion_ESgrid = c("r_only", "all")
  )

  all_results <- list()
  det_rows <- list()
  idx <- 1L

  for (risk in risk_types) {
    apply_vec <- apply_map[[risk]]
    if (is.null(apply_vec)) apply_vec <- dist_apply_default

    for (app in apply_vec) {
      for (fm in forecast_modes) {
        cat("\n\n########################################################\n")
        cat("### World A GRID:", "RISK=", risk, " | mode=", fm,
            " | dist_apply_to=", app,
            " | f_set=", paste(f_set, collapse = ","), "\n", sep = "")
        cat("########################################################\n")

        res <- run_many_A(
          R = R,
          seed0 = seed0,
          risk_type = risk,
          forecast_mode = fm,
          f_set = f_set,
          dist_apply_to = app,
          N_mc = N_mc,
          a_ent = a_ent,
          tau_ex = tau_ex,
          a_lower = a_lower,
          p_short = p_short,
          huber_c = huber_c,
          alpha_evar = alpha_evar,
          u_vec = u_vec,
          w_vec = w_vec,
          oracle_mc_N = oracle_mc_N,
          q_type = q_type
        )

        tag <- paste0(risk, "__", app, "__", fm)
        all_results[[tag]] <- res

        for (nm in names(res)) {
          det_df <- res[[nm]]$det_summary

          if (verbose_print) {
            print_det(det_df, title = paste0("WorldA | ", risk, " | ", app, " | ", nm, " | ", fm))
          }

          f_val <- as.numeric(sub("^f_", "", nm))
          det_df2 <- det_df
          det_df2$risk_type <- risk
          det_df2$dist_apply_to <- app
          det_df2$f <- f_val
          det_df2$forecast_mode <- fm

          det_rows[[idx]] <- det_df2
          idx <- idx + 1L
        }
      }
    }
  }

  det_tbl <- dplyr::bind_rows(det_rows) %>%
    dplyr::select(forecast_mode, risk_type, dist_apply_to, f,
                  threshold, detect_pct, arl, n_detect, mean_time, median_time) %>%
    dplyr::arrange(forecast_mode, risk_type, dist_apply_to, f, threshold)

  invisible(list(summary = det_tbl, raw = all_results))
}


