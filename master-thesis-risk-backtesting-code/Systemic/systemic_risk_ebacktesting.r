###############################################################################
# Systemic risk e-backtesting
#
# Settings:
# - Simulate correlated AR(1)-GARCH(1,1) returns (X_ret, Y_ret) with t innovations
# - Backtest on LOSSES: X = -X_ret, Y = -Y_ret
# - Forecast modes:
#     (1) oracle_true      : uses true latent mu/sigma and true t copula params
#     (2) oracle_empirical : uses presample standardized residual pairs (paired bootstrap)
#     (3) ml_t             : rolling AR(1)-GARCH(1,1) fit with Student-t innovations
# - Sequential aggregation:
#     Components (VaR_under, VaR_over, and one of MES / CoVaR / CoES)
#     Each component uses a GREM adaptive betting fraction lambda_t
#     Combine components via Bonferroni or mixture
#
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

beta_X  <- 0.8   # VaR_X level (tail event on X-losses: X > VaR_X_beta)
alpha_Y <- 0.8   # conditional tail level for Y given X tail (CoVaR/CoES)
a_fixed <- -5    # lower truncation for MES payoff (denominator safety)
delta_global <- 0.2

# GREM settings
J_lambda   <- 500L
lambda_max <- 0.99

# distortion grid defaults
f_set_default <- c(0.9, 1.0, 1.1)
apply_to_set  <- c("none", "r_only", "aux_only", "all")

# ============================================================
# 1) SIMULATION: systemic (X_ret, Y_ret) with correlated t innovations
# ============================================================

simulate_systemic_world <- function(T_total, rho = 0.6, burnin = 1000L, df_innov = 4) {
  # Returns are simulated; later we backtest losses (minus returns).
  if (!is.finite(df_innov) || df_innov <= 2) stop("df_innov must be finite and > 2.")

  T_all <- as.integer(T_total + burnin)

  X_ret <- numeric(T_all); muX <- numeric(T_all); sigX <- numeric(T_all)
  Y_ret <- numeric(T_all); muY <- numeric(T_all); sigY <- numeric(T_all)

  # X GARCH 
  ox <- 0.01; ax <- 0.10; bx <- 0.85
  mu0x <- -0.05; phix <- 0.25

  # Y GARCH 
  oy <- 0.015; ay <- 0.12; by <- 0.80
  mu0y <- -0.03; phiy <- 0.20

  sigX[1] <- sqrt(ox / (1 - ax - bx))
  sigY[1] <- sqrt(oy / (1 - ay - by))
  muX[1]  <- mu0x
  muY[1]  <- mu0y

  # Draw one correlated standardized t pair
  draw_corr_t_pair <- function() {
    u <- rnorm(2)
    z <- c(u[1], rho * u[1] + sqrt(1 - rho^2) * u[2])

    w <- sqrt(df_innov / stats::rchisq(1, df = df_innov)) # common scale
    s <- sqrt((df_innov - 2) / df_innov)                  # makes Var=1
    s * z * w
  }

  z <- draw_corr_t_pair()
  X_ret[1] <- muX[1] + sigX[1] * z[1]
  Y_ret[1] <- muY[1] + sigY[1] * z[2]

  for (t in 2:T_all) {
    # AR(1) means on returns
    muX[t] <- mu0x + phix * X_ret[t - 1]
    muY[t] <- mu0y + phiy * Y_ret[t - 1]

    # standardized residuals 
    zx_prev <- (X_ret[t - 1] - muX[t - 1]) / pmax(sigX[t - 1], 1e-8)
    zy_prev <- (Y_ret[t - 1] - muY[t - 1]) / pmax(sigY[t - 1], 1e-8)

    sigX[t] <- sqrt(pmax(ox + ax * sigX[t - 1]^2 * zx_prev^2 + bx * sigX[t - 1]^2, 0))
    sigY[t] <- sqrt(pmax(oy + ay * sigY[t - 1]^2 * zy_prev^2 + by * sigY[t - 1]^2, 0))

    z <- draw_corr_t_pair()
    X_ret[t] <- muX[t] + sigX[t] * z[1]
    Y_ret[t] <- muY[t] + sigY[t] * z[2]
  }

  idx <- (burnin + 1L):T_all
  list(
    X_ret = X_ret[idx], Y_ret = Y_ret[idx],
    muX = muX[idx], sigX = sigX[idx],
    muY = muY[idx], sigY = sigY[idx],
    rho = rho,
    df_innov = df_innov
  )
}

get_backtest_slice <- function(sim, T_train, T_backtest) {
  idx_bt <- (T_train + 1L):(T_train + T_backtest)
  list(
    idx_bt = idx_bt,
    X_bt_ret = sim$X_ret[idx_bt],
    Y_bt_ret = sim$Y_ret[idx_bt],
    muX_bt = sim$muX[idx_bt], sigX_bt = sim$sigX[idx_bt],
    muY_bt = sim$muY[idx_bt], sigY_bt = sim$sigY[idx_bt]
  )
}

# ============================================================
# 2) ML forecasting: rolling refit AR(1)-GARCH(1,1) with Student-t
# ============================================================

fit_ar1_garch11_t <- function(x) {
  # Student-t innovations 
  spec <- rugarch::ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model     = list(armaOrder = c(1, 0), include.mean = TRUE),
    distribution.model = "std"
  )
  rugarch::ugarchfit(spec = spec, data = x, solver = "hybrid")
}

build_ml_mu_sig_t <- function(series_full, T_train, T_backtest, refit_every = 10L) {
  # Rolling re-fit every 'refit_every' steps, always using a window of length T_train.
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
    if (is.null(fit)) stop("ugarchfit failed (Student-t). Try different seed/refit_every.")

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

    for (j in 1:block_len) {
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

  list(mu = mu_f, sig = sig_f, nu = nu_f, fits = fits[1:fit_idx], refit_every = refit_every)
}

rho_hat_from_fits <- function(fitX, fitY) {
  # Estimate correlation of standardized residuals (paired) inside the current fit block.
  epsX <- as.numeric(residuals(fitX, standardize = TRUE))
  epsY <- as.numeric(residuals(fitY, standardize = TRUE))
  m <- min(length(epsX), length(epsY))
  epsX <- epsX[(length(epsX) - m + 1L):length(epsX)]
  epsY <- epsY[(length(epsY) - m + 1L):length(epsY)]
  ok <- is.finite(epsX) & is.finite(epsY)
  if (sum(ok) < 50) return(0)
  as.numeric(cor(epsX[ok], epsY[ok]))
}

draw_corr_t_pair_given_rho <- function(N, rho_hat, df) {
  # Draw N correlated standardized t innovations with correlation rho_hat.
  N <- as.integer(N)
  rho_hat <- max(min(as.numeric(rho_hat), 0.999), -0.999)
  if (!is.finite(df) || df <= 2) df <- 8

  u1 <- rnorm(N)
  u2 <- rnorm(N)
  zx <- u1
  zy <- rho_hat * u1 + sqrt(1 - rho_hat^2) * u2

  w <- sqrt(df / stats::rchisq(N, df = df))
  s <- sqrt((df - 2) / df)

  list(zx = s * zx * w, zy = s * zy * w)
}

# ============================================================
# 3) Empirical presample residual pairs (oracle_empirical)
# ============================================================

estimate_presample_Z_pairs_loss <- function(sim, presample_length) {
  # We work with LOSSES L = -R. Hence standardized loss innovation is:
  # Z_loss = -(R - mu_R)/sigma_R.
  idx <- 1:as.integer(presample_length)
  Zx <- -(sim$X_ret[idx] - sim$muX[idx]) / sim$sigX[idx]
  Zy <- -(sim$Y_ret[idx] - sim$muY[idx]) / sim$sigY[idx]
  ok <- is.finite(Zx) & is.finite(Zy)
  cbind(Zx[ok], Zy[ok])
}

# For oracle_true conditional sampling: sample (Zx, Zy) given Zx in right tail.
draw_pair_std_t_loss <- function(N, rho, df) {
  # Returns LOSS innovations directly (right tail means large loss).
  u1 <- rnorm(N)
  u2 <- rnorm(N)
  z1 <- u1
  z2 <- rho * u1 + sqrt(1 - rho^2) * u2
  w  <- sqrt(df / stats::rchisq(N, df = df))
  s  <- sqrt((df - 2) / df)
  list(zx = -(s * z1 * w), zy = -(s * z2 * w))
}

sample_tail_pairs_loss <- function(N, beta, rho, df) {
  # Accept-reject sampling for (Zx, Zy) conditional on Zx > q_beta (right tail).
  qx <- sqrt((df - 2) / df) * stats::qt(beta, df = df)  # threshold for LOSS innovation
  zx <- zy <- numeric(0)
  while (length(zx) < N) {
    m <- ceiling((N - length(zx)) * 1.3 / (1 - beta))
    Z <- draw_pair_std_t_loss(m, rho = rho, df = df)
    keep <- which(Z$zx > qx)
    if (length(keep)) {
      zx <- c(zx, Z$zx[keep])
      zy <- c(zy, Z$zy[keep])
    }
  }
  list(zx = zx[1:N], zy = zy[1:N], qx = qx)
}

# ============================================================
# 4) Forecasts 
# ============================================================

forecast_oracle_true <- function(sim, T_train, T_backtest, beta_X, alpha_Y,
                                 test = c("MES", "CoVaR", "CoES"),
                                 N_mc_cond = 5000L) {
  test <- match.arg(test)

  sl <- get_backtest_slice(sim, T_train, T_backtest)

  # Convert RETURN moments to LOSS moments: L = -R
  muX <- -sl$muX_bt
  sigX <- sl$sigX_bt
  muY <- -sl$muY_bt
  sigY <- sl$sigY_bt

  rho <- sim$rho
  df_innov <- sim$df_innov
  if (!is.finite(df_innov) || df_innov <= 2) stop("oracle_true requires df_innov finite > 2.")

  # Closed form VaR for LOSS under standardized t:
  qx <- sqrt((df_innov - 2) / df_innov) * stats::qt(beta_X, df = df_innov)
  VaR_X_beta <- muX + sigX * qx

  aY_lower <- a_fixed
  MES <- CoVaR <- CoES <- rep(NA_real_, T_backtest)

  N_mc_cond <- as.integer(N_mc_cond)

  if (test == "MES") {
    MES <- numeric(T_backtest)
    for (t in 1:T_backtest) {
      Z <- sample_tail_pairs_loss(N_mc_cond, beta_X, rho = rho, df = df_innov)
      Yc <- muY[t] + sigY[t] * Z$zy
      MES[t] <- mean(pmax(Yc, a_fixed))
    }
    return(list(
      muX = muX, sigX = sigX, muY = muY, sigY = sigY, rho = rho,
      VaR_X_beta = VaR_X_beta, MES = MES, CoVaR = CoVaR, CoES = CoES,
      aY_lower = a_fixed, mode = "oracle_true", N_mc_cond = N_mc_cond
    ))
  }

  for (t in 1:T_backtest) {
    Z <- sample_tail_pairs_loss(N_mc_cond, beta_X, rho = rho, df = df_innov)
    Yc <- muY[t] + sigY[t] * Z$zy

    CoVaR[t] <- as.numeric(quantile(Yc, probs = alpha_Y, type = 8))
    if (test == "CoES") {
      idx <- which(Yc > CoVaR[t])
      if (length(idx) < 10) idx <- seq_along(Yc)
      CoES[t] <- mean(Yc[idx])
    }
  }

  list(
    muX = muX, sigX = sigX, muY = muY, sigY = sigY, rho = rho,
    VaR_X_beta = VaR_X_beta, MES = MES, CoVaR = CoVaR, CoES = CoES,
    aY_lower = aY_lower, mode = "oracle_true", N_mc_cond = N_mc_cond
  )
}

forecast_oracle_empirical <- function(sim, T_train, T_backtest, beta_X, alpha_Y,
                                      test = c("MES", "CoVaR", "CoES"),
                                      N_mc = 5000L) {
  test <- match.arg(test)
  sl <- get_backtest_slice(sim, T_train, T_backtest)

  muX <- -sl$muX_bt
  sigX <- sl$sigX_bt
  muY <- -sl$muY_bt
  sigY <- sl$sigY_bt

  Zhat_pair <- estimate_presample_Z_pairs_loss(sim, presample_length = T_train)
  if (nrow(Zhat_pair) < 50) stop("oracle_empirical: too few presample paired residuals.")

  N_mc <- as.integer(N_mc)
  VaR_X_beta <- rep(NA_real_, T_backtest)
  MES <- CoVaR <- CoES <- rep(NA_real_, T_backtest)

  aY_lower <- a_fixed

  for (t in 1:T_backtest) {
    ii <- sample.int(nrow(Zhat_pair), size = N_mc, replace = TRUE)
    zx <- Zhat_pair[ii, 1]
    zy <- Zhat_pair[ii, 2]

    Xs <- muX[t] + sigX[t] * zx
    Ys <- muY[t] + sigY[t] * zy

    VaR_X_beta[t] <- as.numeric(quantile(Xs, probs = beta_X, type = 8))

    idx_tailX <- which(Xs > VaR_X_beta[t])
    if (length(idx_tailX) < 30) idx_tailX <- seq_along(Xs)
    Y_tailX <- Ys[idx_tailX]

    if (test == "MES") {
      MES[t] <- mean(pmax(Y_tailX, a_fixed))
    } else if (test == "CoVaR") {
      CoVaR[t] <- as.numeric(quantile(Y_tailX, probs = alpha_Y, type = 8))
    } else {
      CoVaR[t] <- as.numeric(quantile(Y_tailX, probs = alpha_Y, type = 8))
      idx2 <- which(Y_tailX > CoVaR[t])
      if (length(idx2) < 10) idx2 <- seq_along(Y_tailX)
      CoES[t] <- mean(Y_tailX[idx2])
    }
  }

  list(
    muX = muX, sigX = sigX, muY = muY, sigY = sigY, rho = NA_real_,
    VaR_X_beta = VaR_X_beta, MES = MES, CoVaR = CoVaR, CoES = CoES,
    aY_lower = aY_lower, mode = "oracle_empirical", N_mc = N_mc
  )
}

forecast_ml_t <- function(sim, T_train, T_backtest, beta_X, alpha_Y,
                          test = c("MES", "CoVaR", "CoES"),
                          N_mc = 5000L,
                          refit_every = 10L) {
  test <- match.arg(test)

  stX <- build_ml_mu_sig_t(sim$X_ret, T_train, T_backtest, refit_every = refit_every)
  stY <- build_ml_mu_sig_t(sim$Y_ret, T_train, T_backtest, refit_every = refit_every)

  # Convert return forecasts to loss forecasts
  muX <- -stX$mu
  sigX <- stX$sig
  nuX <- stX$nu
  muY <- -stY$mu
  sigY <- stY$sig
  nuY <- stY$nu

  VaR_X_beta <- rep(NA_real_, T_backtest)
  MES <- CoVaR <- CoES <- rep(NA_real_, T_backtest)

  aY_lower <- a_fixed
  N_mc <- as.integer(N_mc)

  block_id <- ceiling(seq_len(T_backtest) / as.integer(refit_every))
  blocks <- sort(unique(block_id))

  for (b in blocks) {
    idx_b <- which(block_id == b)
    fitX <- stX$fits[[b]]
    fitY <- stY$fits[[b]]
    rho_hat <- rho_hat_from_fits(fitX, fitY)

    for (t in idx_b) {
      # Use a conservative/common df for the joint sampling step
      nu_common <- suppressWarnings(min(c(nuX[t], nuY[t]), na.rm = TRUE))
      Zp <- draw_corr_t_pair_given_rho(N_mc, rho_hat = rho_hat, df = nu_common)

      
      Xs <- muX[t] + sigX[t] * (-Zp$zx)
      Ys <- muY[t] + sigY[t] * (-Zp$zy)

      VaR_X_beta[t] <- as.numeric(quantile(Xs, probs = beta_X, type = 8))

      idx_tailX <- which(Xs > VaR_X_beta[t])
      if (length(idx_tailX) < 30) idx_tailX <- seq_along(Xs)
      Y_tailX <- Ys[idx_tailX]

      if (test == "MES") {
        MES[t] <- mean(pmax(Y_tailX, a_fixed))
      } else if (test == "CoVaR") {
        CoVaR[t] <- as.numeric(quantile(Y_tailX, probs = alpha_Y, type = 8))
      } else {
        CoVaR[t] <- as.numeric(quantile(Y_tailX, probs = alpha_Y, type = 8))
        idx2 <- which(Y_tailX > CoVaR[t])
        if (length(idx2) < 10) idx2 <- seq_along(Y_tailX)
        CoES[t] <- mean(Y_tailX[idx2])
      }
    }
  }

  list(
    muX = muX, sigX = sigX, muY = muY, sigY = sigY, rho = NA_real_,
    VaR_X_beta = VaR_X_beta, MES = MES, CoVaR = CoVaR, CoES = CoES,
    aY_lower = aY_lower, mode = "ml_t", N_mc = N_mc, refit_every = refit_every
  )
}

forecast_dispatch <- function(sim, mode, T_train, T_backtest, beta_X, alpha_Y,
                              test = c("MES", "CoVaR", "CoES"),
                              N_mc = 5000L, refit_every = 10L, N_mc_cond = 5000L) {
  test <- match.arg(test)
  mode <- match.arg(mode, choices = c("oracle_true", "oracle_empirical", "ml_t"))

  if (mode == "oracle_true") {
    return(forecast_oracle_true(
      sim, T_train, T_backtest, beta_X, alpha_Y,
      test = test, N_mc_cond = N_mc_cond
    ))
  }
  if (mode == "oracle_empirical") {
    return(forecast_oracle_empirical(
      sim, T_train, T_backtest, beta_X, alpha_Y,
      test = test, N_mc = N_mc
    ))
  }
  forecast_ml_t(
    sim, T_train, T_backtest, beta_X, alpha_Y,
    test = test, N_mc = N_mc, refit_every = refit_every
  )
}

# ============================================================
# 5) Distortion
# ============================================================

scale_abs <- function(x, f) x + (f - 1) * abs(x)   # f=0.9 down, f=1.1 up

apply_distortion <- function(fc, f = 1.0,
                             apply_to = c("none", "r_only", "aux_only", "all"),
                             test = c("MES", "CoVaR", "CoES")) {
  apply_to <- match.arg(apply_to)
  test <- match.arg(test)
  out <- fc
  out$dist <- list(f = f, apply_to = apply_to, test = test)

  if (!is.finite(f) || f <= 0) stop("f must be > 0")
  if (apply_to == "none") return(out)

  scale_aux <- apply_to %in% c("aux_only", "all")
  scale_r   <- apply_to %in% c("r_only", "all")

  # Auxiliaries:
  # - always VaR_X_beta
  # - for CoES test, CoVaR is also auxiliary
  if (scale_aux) {
    out$VaR_X_beta <- scale_abs(out$VaR_X_beta, f)
    if (test == "CoES") out$CoVaR <- scale_abs(out$CoVaR, f)
  }

  # Reported target r:
  if (scale_r) {
    if (test == "MES")   out$MES   <- scale_abs(out$MES, f)
    if (test == "CoVaR") out$CoVaR <- scale_abs(out$CoVaR, f)
    if (test == "CoES")  out$CoES  <- scale_abs(out$CoES, f)
  }

  
  out
}

# ============================================================
# 6) One-step e-values 
# ============================================================

# VaR correctness components (two-sided via two one-sided e-values)
e_var_under <- function(X, z, beta) as.numeric(X > z) / (1 - beta)  # detects z too low
e_var_over  <- function(X, z, beta) as.numeric(X <= z) / beta       # detects z too high

# CoVaR (conditional quantile) e-value:
#   e = 1 + 1{X>zX} (alpha - 1{Y <= r})
e_covar <- function(X, Y, zX, r, alpha) {
  pmax(0, 1 + as.numeric(X > zX) * (alpha - as.numeric(Y <= r)))
}

# MES inequality e-value:
#   e = 1 + 1{X>zX} (Y - r)/(r - a)
# If r <= a, denominator is nonpositive -> we epsilon-correct r.
e_mes <- function(X, Y, zX, r, aY_lower, eps_den = 1e-6) {
  if (!is.finite(r)) return(Inf)
  if (r <= aY_lower + eps_den) r <- aY_lower + eps_den
  den <- (r - aY_lower)
  pmax(0, 1 + as.numeric(X > zX) * (Y - r) / den)
}

# CoES inequality e-value:
#   e = 1{X>zX} (Y - CoVaR)_+ / ( (1-beta)(1-alpha)(r - CoVaR) )
# If r <= CoVaR, we set epsilon-corrected r; if r is nonfinite -> Inf.
e_coes <- function(X, Y, zX, CoVaR, r, beta, alpha, eps_den = 1e-6) {
  if (!is.finite(r) || !is.finite(CoVaR)) return(Inf)
  if (r <= CoVaR + eps_den) r <- CoVaR + eps_den
  den <- (1 - beta) * (1 - alpha) * (r - CoVaR)
  as.numeric(X > zX) * pmax(Y - CoVaR, 0) / den
}

# GREL
past_e_recalibrated <- function(component, idx0, t, X, Y, fc, beta_X, alpha_Y) {
  Xp <- X[idx0]
  Yp <- Y[idx0]

  if (component == "VaR_under") {
    zt <- fc$VaR_X_beta[t]
    return(as.numeric(Xp > zt) / (1 - beta_X))
  }
  if (component == "VaR_over") {
    zt <- fc$VaR_X_beta[t]
    return(as.numeric(Xp <= zt) / beta_X)
  }
  if (component == "MES") {
    zt <- fc$VaR_X_beta[t]
    rt <- fc$MES[t]
    if (!is.finite(rt)) return(rep(Inf, length(idx0)))
    if (rt <= fc$aY_lower + 1e-6) rt <- fc$aY_lower + 1e-6
    den <- (rt - fc$aY_lower)
    return(pmax(0, 1 + as.numeric(Xp > zt) * (Yp - rt) / den))
  }
  if (component == "CoVaR") {
    zt <- fc$VaR_X_beta[t]
    rt <- fc$CoVaR[t]
    return(pmax(0, 1 + as.numeric(Xp > zt) * (alpha_Y - as.numeric(Yp <= rt))))
  }
  if (component == "CoES") {
    zt  <- fc$VaR_X_beta[t]
    cvt <- fc$CoVaR[t]
    rt  <- fc$CoES[t]
    if (!is.finite(rt) || !is.finite(cvt)) return(rep(Inf, length(idx0)))
    if (rt <= cvt + 1e-6) rt <- cvt + 1e-6
    den <- (1 - beta_X) * (1 - alpha_Y) * (rt - cvt)
    return(as.numeric(Xp > zt) * pmax(Yp - cvt, 0) / den)
  }

  stop("Unknown component: ", component)
}

# ============================================================
# 7) GREM (S0): lambda via Taylor 
# ============================================================

lambda_from_e_taylor <- function(e_vec, lambda_max = 1) {
  # Taylor lambda approx: lam = mean(e-1)/mean((e-1)^2), clipped to [0, lambda_max].
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
  # Update log-wealth for multiplicative factor: (1-lam) + lam*e_t
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
  # Numerically stable weight for mixing two experts by their log-wealths.
  if (is.infinite(logA) && logA < 0 && is.infinite(logB) && logB < 0) return(0.5)
  if (is.infinite(logA) && logA < 0) return(0)
  if (is.infinite(logB) && logB < 0) return(1)
  1 / (1 + exp(logB - logA))
}

run_component_grem_S0 <- function(component, e_series, X, Y, fc,
                                  beta_X, alpha_Y,
                                  J_lambda = 500L, lambda_max = 0.99) {


  T <- length(e_series)
  J <- as.integer(J_lambda)
  if (J <= 0) stop("J_lambda must be positive.")

  logW_real <- numeric(T + 1L)
  logW_real[1] <- 0
  logW_recal <- numeric(T + 1L)
  logW_recal[1] <- 0
  logW_mix <- numeric(T + 1L)
  logW_mix[1] <- 0

  lam_real <- numeric(T)
  lam_recal <- numeric(T)
  lam_mix <- numeric(T)

  for (t in 1:T) {
    past_end <- t - 1L

    if (past_end <= 0L) {
      lam_real[t] <- 0
      lam_recal[t] <- 0
    } else {
      idx0 <- max(1L, past_end - J + 1L):past_end

      past_real <- e_series[idx0]
      past_recal <- past_e_recalibrated(component, idx0, t, X, Y, fc, beta_X, alpha_Y)

      lam_real[t]  <- lambda_from_e_taylor(past_real,  lambda_max = lambda_max)
      lam_recal[t] <- lambda_from_e_taylor(past_recal, lambda_max = lambda_max)
    }

   
    logW_real[t + 1L]  <- update_log_wealth(logW_real[t],  e_series[t], lam_real[t])
    logW_recal[t + 1L] <- update_log_wealth(logW_recal[t], e_series[t], lam_recal[t])

    w <- weight_from_logs(logW_real[t], logW_recal[t])
    lam_mix[t] <- w * lam_real[t] + (1 - w) * lam_recal[t]
    logW_mix[t + 1L] <- update_log_wealth(logW_mix[t], e_series[t], lam_mix[t])
  }

  W_mix <- exp(pmin(logW_mix, 700))
  list(
    logW = logW_mix, W = W_mix, lambda = lam_mix,
    lambda_real = lam_real, lambda_recal = lam_recal,
    logW_real = logW_real, logW_recal = logW_recal,
    e = e_series, component = component
  )
}

# ============================================================
# 8) Combine components (Bonferroni + Mixture)
# ============================================================

logsumexp <- function(v) {
  m <- max(v)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(v - m)))
}

combine_bonferroni_delta <- function(logW_list, delta = 0.05, weights = NULL) {
  comp_names <- names(logW_list)
  K <- length(comp_names)
  if (is.null(weights)) weights <- rep(1 / K, K)
  weights <- weights / sum(weights)

  # Bonferroni split of alpha: delta_i = delta * w_i
  delta_i <- delta * weights

  hit_t <- rep(Inf, K)
  for (k in 1:K) {
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
      component = comp_names, weight = weights, delta_i = delta_i, detect_time = hit_t
    )
  )
}

combine_mixture_delta <- function(logW_list, delta = 0.05, weights = NULL) {
  comp_names <- names(logW_list)
  K <- length(comp_names)
  if (is.null(weights)) weights <- rep(1 / K, K)
  weights <- weights / sum(weights)

  T <- length(logW_list[[1]]) - 1L
  logWmix <- rep(-Inf, T + 1L)

  for (t in 1:(T + 1L)) {
    logs <- sapply(1:K, function(k) log(weights[k]) + logW_list[[k]][t])
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

  # Threshold split of c across components:
  # require W_k >= c^(1 / w_k)  <=> logW_k >= (1 / w_k) log c
  thr_k <- (1 / weights) * log(c)

  hit_t <- rep(Inf, K)
  for (k in 1:K) {
    hit <- which(logW_list[[k]] >= thr_k[k])
    if (length(hit)) hit_t[k] <- hit[1] - 1L
  }

  t_star <- min(hit_t)
  which_min <- which(hit_t == t_star)

  list(
    detected = is.finite(t_star),
    t = t_star,
    who = if (is.finite(t_star)) comp_names[which_min] else character(0),
    per_component = data.frame(component = comp_names, weight = weights, thr_log = thr_k, detect_time = hit_t),
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

  for (t in 1:(T + 1L)) {
    logs <- sapply(1:K, function(k) log(weights[k]) + logW_list[[k]][t])
    logWmix[t] <- logsumexp(logs)
  }

  thr <- log(c)
  hit <- which(logWmix >= thr)
  t_star <- if (length(hit)) hit[1] - 1L else Inf

  list(detected = is.finite(t_star), t = t_star, logWmix = logWmix, c = c)
}

# ============================================================
# 9) One replication 
# ============================================================

components_for_test <- function(test = c("MES", "CoVaR", "CoES")) {
  test <- match.arg(test)
  base <- c("VaR_under", "VaR_over")
  if (test == "MES")   return(c(base, "MES"))
  if (test == "CoVaR") return(c(base, "CoVaR"))
  c(base, "CoES")
}

run_one_replication <- function(test = c("MES", "CoVaR", "CoES"),
                                mode = c("oracle_true", "oracle_empirical", "ml_t"),
                                f = 1.0,
                                apply_to = c("none", "r_only", "aux_only", "all"),
                                delta = 0.05,
                                N_mc = 5000L,
                                N_mc_cond = 5000L,
                                refit_every = 10L,
                                seed = 1,
                                sim_rho = 0.6,
                                sim_df  = 4) {
  test <- match.arg(test)
  mode <- match.arg(mode)
  apply_to <- match.arg(apply_to)

  set.seed(seed)

  sim <- simulate_systemic_world(T_total, rho = sim_rho, df_innov = sim_df)
  sl  <- get_backtest_slice(sim, T_train, T_backtest)

  # Backtest is on LOSSES
  X <- -sl$X_bt_ret
  Y <- -sl$Y_bt_ret

  # For MES we cap Y from below 
  Y_cap  <- pmax(Y, a_fixed)
  Y_used <- if (test == "MES") Y_cap else Y

  fc <- forecast_dispatch(
    sim, mode = mode, T_train = T_train, T_backtest = T_backtest,
    beta_X = beta_X, alpha_Y = alpha_Y, test = test,
    N_mc = N_mc, refit_every = refit_every, N_mc_cond = N_mc_cond
  )

  fc <- apply_distortion(fc, f = f, apply_to = apply_to, test = test)

  # Build component e-series
  comps <- components_for_test(test)
  e_series <- list()

  if ("VaR_under" %in% comps) {
    e_series$VaR_under <- sapply(1:T_backtest, function(t)
      e_var_under(X[t], fc$VaR_X_beta[t], beta = beta_X)
    )
  }
  if ("VaR_over" %in% comps) {
    e_series$VaR_over <- sapply(1:T_backtest, function(t)
      e_var_over(X[t], fc$VaR_X_beta[t], beta = beta_X)
    )
  }

  if (test == "MES") {
    fc$aY_lower <- a_fixed
    e_series$MES <- sapply(1:T_backtest, function(t)
      e_mes(X[t], Y_cap[t], zX = fc$VaR_X_beta[t], r = fc$MES[t], aY_lower = fc$aY_lower)
    )
  }

  if (test == "CoVaR") {
    e_series$CoVaR <- sapply(1:T_backtest, function(t)
      e_covar(X[t], Y[t], zX = fc$VaR_X_beta[t], r = fc$CoVaR[t], alpha = alpha_Y)
    )
  }

  if (test == "CoES") {
    e_series$CoES <- sapply(1:T_backtest, function(t)
      e_coes(X[t], Y[t], zX = fc$VaR_X_beta[t], CoVaR = fc$CoVaR[t],
             r = fc$CoES[t], beta = beta_X, alpha = alpha_Y)
    )
  }

  # Run GREM per component
  proc <- lapply(names(e_series), function(nm) {
    run_component_grem_S0(
      component = nm,
      e_series  = e_series[[nm]],
      X = X, Y = Y_used, fc = fc,
      beta_X = beta_X, alpha_Y = alpha_Y,
      J_lambda = J_lambda, lambda_max = lambda_max
    )
  })
  names(proc) <- names(e_series)

  logW_list <- lapply(proc, function(p) p$logW)

  # Equal weights across components 
  K <- length(logW_list)
  w <- rep(1 / K, K)
  names(w) <- names(logW_list)

  bonf <- combine_bonferroni_delta(logW_list, delta = delta, weights = w)
  mix  <- combine_mixture_delta(logW_list,  delta = delta, weights = w)

  bonf_c2  <- combine_bonferroni_c(logW_list, c = 2,  weights = w)
  bonf_c5  <- combine_bonferroni_c(logW_list, c = 5,  weights = w)
  bonf_c10 <- combine_bonferroni_c(logW_list, c = 10, weights = w)

  mix_c2   <- combine_mixture_c(logW_list, c = 2,  weights = w)
  mix_c5   <- combine_mixture_c(logW_list, c = 5,  weights = w)
  mix_c10  <- combine_mixture_c(logW_list, c = 10, weights = w)

  avg_fc <- list(
    avg_VaR_X = mean(fc$VaR_X_beta, na.rm = TRUE),
    avg_MES   = if (!all(is.na(fc$MES)))   mean(fc$MES, na.rm = TRUE)   else NA_real_,
    avg_CoVaR = if (!all(is.na(fc$CoVaR))) mean(fc$CoVaR, na.rm = TRUE) else NA_real_,
    avg_CoES  = if (!all(is.na(fc$CoES)))  mean(fc$CoES, na.rm = TRUE)  else NA_real_
  )

  list(
    meta = list(
      test = test, mode = mode, f = f, apply_to = apply_to, delta = delta,
      N_mc = N_mc, N_mc_cond = N_mc_cond, refit_every = refit_every,
      sim_rho = sim_rho, sim_df = sim_df
    ),
    forecasts = fc,
    avg_forecasts = avg_fc,
    components = list(e = e_series, proc = proc, logW = logW_list, weights = w),
    results = list(
      bonferroni = bonf,
      mixture    = mix,
      bonf_c = list(c2 = bonf_c2, c5 = bonf_c5, c10 = bonf_c10),
      mix_c  = list(c2 = mix_c2,  c5 = mix_c5,  c10 = mix_c10)
    )
  )
}

# ============================================================
# 10) Many reps + summaries
# ============================================================

summarize_detection <- function(times) {
  hit <- is.finite(times)
  data.frame(
    detect_pct = 100 * mean(hit),
    n_detect = sum(hit),
    mean_time = if (sum(hit) > 0) mean(times[hit]) else NA_real_,
    median_time = if (sum(hit) > 0) median(times[hit]) else NA_real_
  )
}

run_many_replications <- function(R = 200,
                                  seed0 = 1,
                                  test = c("MES", "CoVaR", "CoES"),
                                  mode = c("oracle_true", "oracle_empirical", "ml_t"),
                                  f_set = f_set_default,
                                  apply_to = "r_only",
                                  delta = delta_global,
                                  N_mc = 3000L,
                                  N_mc_cond = 3000L,
                                  refit_every = 10L,
                                  sim_rho = 0.6,
                                  sim_df  = 4) {
  test <- match.arg(test)
  mode <- match.arg(mode)
  apply_to <- match.arg(apply_to, choices = apply_to_set)

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

    avg_VaR <- avg_r <- avg_aux <- rep(NA_real_, R)

    for (r in 1:R) {
      rep_out <- run_one_replication(
        test = test,
        mode = mode,
        f = f,
        apply_to = apply_to,
        delta = delta,
        N_mc = N_mc,
        N_mc_cond = N_mc_cond,
        refit_every = refit_every,
        seed = seed0 + r,
        sim_rho = sim_rho,
        sim_df  = sim_df
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

      avg_VaR[r] <- rep_out$avg_forecasts$avg_VaR_X
      if (test == "MES")   avg_r[r] <- rep_out$avg_forecasts$avg_MES
      if (test == "CoVaR") avg_r[r] <- rep_out$avg_forecasts$avg_CoVaR
      if (test == "CoES")  avg_r[r] <- rep_out$avg_forecasts$avg_CoES
      if (test == "CoES")  avg_aux[r] <- rep_out$avg_forecasts$avg_CoVaR
    }

    avg_tbl <- data.frame(
      avg_VaR_X = mean(avg_VaR, na.rm = TRUE),
      avg_r     = mean(avg_r, na.rm = TRUE),
      avg_aux_CoVaR = if (test == "CoES") mean(avg_aux, na.rm = TRUE) else NA_real_
    )

    out[[k]] <- list(
      settings = list(
        R = R, test = test, mode = mode, f = f, apply_to = apply_to, delta = delta,
        N_mc = N_mc, N_mc_cond = N_mc_cond, refit_every = refit_every,
        sim_rho = sim_rho, sim_df = sim_df
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

    cat("\n=== TEST =", test, "| mode =", mode, "| f =", f, "| apply_to =", apply_to, "===\n")
    cat("Bonferroni (delta):\n"); print(out[[k]]$bonf$summary)
    cat("Mixture (delta):\n");    print(out[[k]]$mix$summary)
    cat("Avg forecasts:\n");      print(out[[k]]$avg_forecasts)
    cat("Bonferroni (c=2/5/10):\n")
    print(out[[k]]$bonf_c$c2$summary); print(out[[k]]$bonf_c$c5$summary); print(out[[k]]$bonf_c$c10$summary)
    cat("Mixture (c=2/5/10):\n")
    print(out[[k]]$mix_c$c2$summary);  print(out[[k]]$mix_c$c5$summary);  print(out[[k]]$mix_c$c10$summary)
  }

  out
}
