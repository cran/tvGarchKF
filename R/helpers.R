#' @importFrom utils modifyList
#' @importFrom stats approx
# ===================================================
# ===================== Helpers =====================
# ===================================================
.normalize_tv_names <- function(df) {
  stopifnot(is.data.frame(df))
  cn <- names(df)
  if (is.null(cn)) stop("'tvparameter' should have columns names.")

  num_idx <- which(vapply(df, is.numeric, logical(1)))
  if (length(num_idx) < 4L) {
    stop("Need at least 4 numeric columns in 'tvparameter' (t, c, alpha, beta).")
  }
  lc <- tolower(gsub("\\s+", "", cn))

  find_by_pat <- function(pats, exclude = integer(0)) {
    idx <- setdiff(seq_along(lc), exclude)
    for (p in pats) {
      hit <- grep(p, lc[idx], perl = TRUE)
      if (length(hit) == 1L) return(idx[hit])
      if (length(hit) > 1L)  return(idx[hit[1L]])
    }
    NA_integer_
  }

  t_pats <- c("\u005e t$", "\u005et1$", "\u005e time(_?index)?$", "\u005e time$", "\u005e idx$", "\u005e index$", "\u005e fecha$", "\u005e tiempo$")
  o_pats <- c("\u005e c?$", "\u005e w1?$", "\u005e om(ega)?$")
  a_pats <- c("\u005e alpha?$", "\u005e a1?$", "\u005e al(pha)?$")
  b_pats <- c("\u005e beta?$", "\u005e b1?$", "\u005e be(ta)?$")

  used <- integer(0)

  it <- find_by_pat(t_pats)
  if (is.na(it)) {
    mono <- sapply(num_idx, function(j) {
      x <- as.numeric(df[[j]]); if (length(x) > 1) mean(diff(x) > 0, na.rm = TRUE) else 0
    })
    it <- num_idx[which.max(mono)]
    if (all(!is.finite(mono)) || mono[which.max(mono)] < 0.6) {
      warning("It's not detect column of time wrong; using '", cn[it], "'.")
    }
  }
  used <- c(used, it)

  io <- find_by_pat(o_pats, used)
  ia <- find_by_pat(a_pats, c(used, io))
  ib <- find_by_pat(b_pats, c(used, io, ia))

  rem <- setdiff(num_idx, c(it, io, ia, ib))
  pick_next <- function() { if (length(rem)) { j <- rem[1]; rem <<- rem[-1]; j } else NA_integer_ }
  if (is.na(io)) io <- pick_next()
  if (is.na(ia)) ia <- pick_next()
  if (is.na(ib)) ib <- pick_next()

  if (any(is.na(c(it, io, ia, ib)))) {
    stop("Could not automatically identify t/c/alpha/beta.",
         "Rename the columns or make sure you have \u2265 4 numeric.")
  }

  out <- data.frame(
    t      = as.numeric(df[[it]]),
    c  = as.numeric(df[[io]]),
    alpha  = as.numeric(df[[ia]]),
    beta   = as.numeric(df[[ib]]),
    check.names = FALSE
  )
  attr(out, "mapping") <- c(t = cn[it], c = cn[io], alpha = cn[ia], beta = cn[ib])
  out
}

.merge_opts <- function(user, defaults) {
  if (is.null(user)) return(defaults)
  modifyList(defaults, user, keep.null = TRUE)
}
.get_arg <- function(dots, name, default) {
  if (!is.null(dots[[name]])) dots[[name]] else default
}
.pred_ss_vec <- function(fit, x) {
  pr <- try(predict(fit, x), silent = TRUE)
  if (inherits(pr, "try-error") || is.null(pr)) return(rep(NA_real_, length(x)))
  if (is.list(pr) && !is.null(pr$y)) as.numeric(pr$y) else as.numeric(pr)
}

.stationarity_metrics <- function(alpha, beta, tol = 1e-3) {
  ab <- alpha + beta
  max_ab <- suppressWarnings(max(ab, na.rm = TRUE))
  viol_idx <- which(ab >= (1 - tol))
  n_viol <- length(viol_idx)
  prop_viol <- n_viol / length(ab)
  list(
    ok            = (n_viol == 0),
    max_alpha_beta= max_ab,
    n_viol        = n_viol,
    prop_viol     = prop_viol,
    first_viol    = if (n_viol) viol_idx[1] else NA_integer_,
    tol           = tol
  )
}

.norm_names_if <- function(tvparameter) {
  if (exists(".normalize_tv_names", mode = "function")) {
    tryCatch(.normalize_tv_names(tvparameter), error = function(e) tvparameter)
  } else tvparameter
}

.merge_opts <- function(user, defaults) {
  if (is.null(user)) return(defaults)
  modifyList(defaults, user, keep.null = TRUE)
}

get_df_grid <- function(df_grid, name, default) {
  if (is.list(df_grid)) {
    if (!is.null(df_grid[[name]])) return(df_grid[[name]])
    return(default)
  }
  df_grid
}


eval_model <- function(tvparameter, Fittvparameter, na.rm = TRUE) {
  stopifnot(is.data.frame(tvparameter), is.data.frame(Fittvparameter))

  num_tv  <- which(vapply(tvparameter, is.numeric, logical(1)))
  num_fit <- which(vapply(Fittvparameter, is.numeric, logical(1)))
  if (length(num_tv) == 0L || length(num_fit) == 0L) {
    stop("Both data frames must contain numeric columns.")
  }

  k <- min(length(num_tv), length(num_fit))
  true <- as.matrix(tvparameter[ , num_tv[seq_len(k)], drop = FALSE])
  pred <- as.matrix(Fittvparameter[ , num_fit[seq_len(k)], drop = FALSE])

  n <- min(nrow(true), nrow(pred))
  true <- true[seq_len(n), , drop = FALSE]
  pred <- pred[seq_len(n), , drop = FALSE]

  res <- true - pred

  rmse_col <- sqrt(colMeans(res^2, na.rm = na.rm))
  mae_col  <- colMeans(abs(res),   na.rm = na.rm)

  rmse_macro <- mean(rmse_col, na.rm = na.rm)
  mae_macro  <- mean(mae_col,  na.rm = na.rm)
  rmse_micro <- sqrt(mean(res^2, na.rm = na.rm))
  mae_micro  <- mean(abs(res),  na.rm = na.rm)

  param_names <- c("c","alpha","beta")[seq_len(k)]

  list(
    per_param = data.frame(
      param = param_names,
      RMSE  = as.numeric(rmse_col),
      MAE   = as.numeric(mae_col)
    ),
    overall = c(
      RMSE_macro = rmse_macro,
      MAE_macro  = mae_macro,
      RMSE_micro = rmse_micro,
      MAE_micro  = mae_micro
    )
  )
}
print_eval_model <- function(x, digits = 4) {
  cat("Per Parameter\n")
  cat("-------------\n")
  print(x$per_param, row.names = FALSE, digits = digits)

  cat("\nOverall\n")
  cat("-------\n")
  overall_df <- as.data.frame(t(x$overall))
  ord <- c("RMSE_macro","MAE_macro","RMSE_micro","MAE_micro")
  overall_df <- overall_df[, intersect(ord, names(overall_df)), drop = FALSE]
  print(overall_df, row.names = FALSE, digits = digits)
}


#' @export
print.auto.tvGarchKF <- function(x, ...) {
  cat("Stationarity check (alpha + beta < 1 - tol) by method\n")
  cat("----------------------------------------------------\n")
  df <- x$summary
  cr <- x$criteria
  df$ok <- ifelse(df$ok, "OK", "FAIL")
  df$prop_viol <- sprintf("%.1f%%", 100 * df$prop_viol)
  print(df, row.names = FALSE)
  if (any(df$ok == "FAIL")) {
    worst <- df[which.max(df$max_alpha_beta), , drop = FALSE]
    cat(sprintf("\nWorst max(alpha+beta)=%.4f in method: %s (tol=%.1e)\n",
                worst$max_alpha_beta, worst$method, x$tol))
  } else {
    cat(sprintf("\nOK: stationarity satisfied for all methods (tol=%.1e)\n", x$tol))
  }
  invisible(x)
}

#' @method summary auto.tvGarchKF
#' @export
summary.auto.tvGarchKF <- function(object, ...) {

  cat("\n")
  cat("Automatic Non-parametric Parameter Selection for tv-GARCH(1,1)\n")
  cat("==============================================================\n\n")

  cat("Input information\n")
  cat("-----------------\n")
  cat("Input class:", paste(object$class_input, collapse = ", "), "\n")
  cat("Number of observations:", object$n, "\n")
  cat("Tolerance:", object$tol, "\n\n")

  gam_ref <- object$gam_reference

  cat("GAM reference parameters\n")
  cat("------------------------\n")
  print(gam_ref, row.names = FALSE)

  cat("\n")
  cat("Interpretation: GAM provides a reference level of smoothness and complexity.\n")
  cat("Recommended intervals are built using the selected values as upper bounds,\n")
  cat("except for smoothing-spline spar, where the upper bound is 2 times the selected value.\n\n")


  # ======================================================
  # Smoothing spline intervals
  # ======================================================

  ss <- object$selected_parameters$smoothing_spline

  spar_lower <- 0.30
  spar_upper <- pmin(1.3, 2 * ss$spar)

  ss_out <- data.frame(
    parameter = ss$parameter,
    `edf-spar` = round(ss$spar, 4),
    `recommended spar` = paste0(
      "[",
      spar_lower,
      ", ",
      round(spar_upper, 4),
      "]"
    ),
    `nknots-selected`  = ss$nknots,
    `recommended nknots` = paste0(
      "[3, ",
      ss$nknots,
      "]"
    ),
    check.names = FALSE
  )


  cat("Smoothing Spline\n")
  cat("----------------------------------------------------\n")
  print(ss_out, row.names = FALSE)

  cat("\n")

  # ======================================================
  # Local polynomial intervals
  # ======================================================

  lp <- object$selected_parameters$local_polynomial

  lp_out <- data.frame(
    parameter = lp$parameter,
    `edf-span` = round(lp$span, 4),
    `recommended span` = paste0(
      "[0.50, ",
      round(2 * lp$span, 4),
      "]"
    ),
    `degree-selected` = lp$degree,
    `recommended degree` = paste0(
      "[1, ",
      lp$degree,
      "]"
    ),
    check.names = FALSE
  )

  cat("Local Polynomial\n")
  cat("----------------------------------------------------\n")
  print(lp_out, row.names = FALSE)

  cat("\n")

  # ======================================================
  # B-spline intervals
  # ======================================================

  bs <- object$selected_parameters$b_spline

  bs_out <- data.frame(
    parameter = bs$parameter,
    `df-selected` = bs$df,
    `recommended df` = paste0(
      "[3, ",
      bs$df,
      "]"
    ),
    `degree-selected` = bs$degree,
    `recommended degree` = paste0(
      "[1, ",
      bs$degree,
      "]"
    ),
    check.names = FALSE
  )

  cat("B-Spline\n")
  cat("----------------------------------------------------\n")
  print(bs_out, row.names = FALSE)

  cat("\n")
  cat("Note:\n")
  cat("For smoothing splines, the upper bound for spar is min(1.3, 2 * spar_selected).\n")
  cat("This allows smoother curves when the EDF-based matching returns a low spar value.\n")

  invisible(object)
}


##########################
#### PLOT HELPERS ####
##########################

plot_tvGarchKalmanFit <- function(fit,
                                  tv_curve,
                                  type = c("linear", "linear", "quadratic"),
                                  trig = "cos",
                                  arg = "2*pi*u",
                                  col.points = "black",
                                  col.fit = "black",
                                  col.stationary = "black",
                                  lwd.fit = 3,
                                  lwd.stationary = 2,
                                  lty.fit = 1,
                                  lty.stationary = 2,
                                  pch = 20,
                                  cex.points = 3,
                                  cex.axis = 3,
                                  cex.main = 3,
                                  mfrow = c(1, 3),
                                  mar = c(5, 6, 4, 2),margen=0.01) {

  coef_tv <- tv_curve@coef
  n       <- length(fit@sigma2)
  t       <- coef_tv$t
  c       <- coef_tv$c
  alpha   <- coef_tv$alpha
  beta    <- coef_tv$beta
  t_plot = (t-min(t))/(max(t)-min(t))

  u <- t / n

  G <- tv_curve@stationary_coef

  fit <- as.numeric(fit@coef)
  idx <- 1
  build_curve <- function(model_type, u) {

    model_type <- tolower(model_type)
    if (model_type == "constant") {
      y <- fit[idx] + 0 * u
      idx <<- idx + 1

    } else if (model_type == "linear") {

      y <- fit[idx] + fit[idx + 1] * u
      idx <<- idx + 2

    } else if (model_type == "quadratic") {

      y <- fit[idx] +
        fit[idx + 1] * u +
        fit[idx + 2] * u^2

      idx <<- idx + 3

    } else if (model_type == "cubic") {

      y <- fit[idx] +
        fit[idx + 1] * u +
        fit[idx + 2] * u^2 +
        fit[idx + 3] * u^3

      idx <<- idx + 4

    } else if (model_type == "trigonometric") {

      z <- eval(parse(text = arg), envir = list(u = u, pi = pi))

      if (trig == "cos") {
        y <- fit[idx] + fit[idx + 1] * cos(z)
      } else {
        y <- fit[idx] + fit[idx + 1] * sin(z)
      }

      idx <<- idx + 2

    } else if (model_type %in% c("nonlinear", "nolineal")) {

      y <- fit[idx] +
        fit[idx + 1] * u +
        fit[idx + 2] * sqrt(u)

      idx <<- idx + 3

    } else if (model_type %in% c("coslog", "logcos")) {

      u_safe <- pmax(u, .Machine$double.eps)
      y <- fit[idx] +
        fit[idx + 1] * cos(3 * (1 - log(u_safe)))
      idx <<- idx + 2

    } else {
      stop("Unknown model type: ", model_type)
    }
    y
  }

  hat_c     <- build_curve(type[1], u)
  hat_alpha <- build_curve(type[2], u)
  hat_beta  <- build_curve(type[3], u)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow = mfrow, mar = mar)

  plot_one <- function(y_true, y_hat, g, main_expr,margen) {
    ylim <- c(
      min(c(y_true, y_hat,g), na.rm = TRUE) - margen,
      max(c(y_true, y_hat, g), na.rm = TRUE) + margen
    )
    plot(t_plot, y_true,
         xlim = c(0, max(t_plot)), #c(0,n)
         ylim = ylim,
         bty = "n",
         main = main_expr,
         las = 1,
         ylab = "",
         xlab = "",
         pch = pch,
         cex = cex.points,
         cex.axis = cex.axis,
         cex.main = cex.main,
         col = col.points)
    lines(t_plot, y_hat,
          col = col.fit,
          lwd = lwd.fit,
          lty = lty.fit)

    abline(h = g,
           col = col.stationary,
           lwd = lwd.stationary,
           lty = lty.stationary)
  }

  plot_one(c, hat_c, G[1], expression(c(u)), margen)
  mtext(side=1, text = "u", line=3, cex=1.5)
  plot_one(alpha, hat_alpha, G[2], expression(alpha(u)), margen)
  mtext(side=1, text = "u", line=3, cex=1.5)
  plot_one(beta,  hat_beta,  G[3], expression(beta(u)), margen)
  mtext(side=1, text = "u", line=3, cex=1.5)

  invisible(data.frame(
    t = t,
    u = u,
    c_hat = hat_c,
    alpha_hat = hat_alpha,
    beta_hat = hat_beta
  ))
}


plot_tvGarchNonParamFit <- function(object,
                                    tv_curve = NULL,
                                    method = NULL,
                                    col.points = "black",
                                    col.fit = "black",
                                    col.stationary = "black",
                                    lwd.fit = 3,
                                    lwd.stationary = 2,
                                    lty.fit = 1,
                                    lty.stationary = 2,
                                    pch = 20,
                                    cex.points = 3,
                                    cex.axis = 3,
                                    cex.main = 3,
                                    mfrow = c(1, 3),
                                    mar = c(5, 6, 4, 2), margen=0.01) {

  ## --------------------------------------------------
  ## 1. Extraer curvas heuristicas por bloques
  ## --------------------------------------------------
  if (is.null(tv_curve)) {
    if (!is.null(object@tvparameter)) {
      tv_curve <- object@tvparameter
    } else if (!is.null(object@TVparameter)) {
      tv_curve <- object@TVparameter
    } else {
      stop("Please provide 'tv_curve' or include it inside 'object'.")
    }
  }

  coef_tv <- tv_curve@coef

  n     <- length(object@data)
  t     <- coef_tv$t
  c     <- coef_tv$c
  alpha <- coef_tv$alpha
  beta  <- coef_tv$beta
  u     <- t/n
  t_plot <- (t-min(t))/(max(t)-min(t))

  G <- tv_curve@stationary_coef

  ## --------------------------------------------------
  ## 2. Extraer valores ajustados no parametricos
  ## --------------------------------------------------
  if (!is.null(method)) {
    method <- match.arg(
      method,
      choices = c("Smoothing Spline", "Local Polynomial", "B-Spline")
    )
  }

  get_fitted <- function(object, method = NULL) {

    ## Caso 1: object$fitted ya contiene omega, alpha, beta
    if (!is.null(object@fitted.values) &&
        all(c("c", "alpha", "beta") %in% names(object@fitted.values))) {
      return(object@fitted.values[, c("c", "alpha", "beta")])
    }

    ## Caso 2: object tiene lista por metodos
    if (!is.null(method) && !is.null(object[[method]]@fitted)) {
      return(object[[method]]@fitted[, c("c", "alpha", "beta")])
    }

    ## Caso 3: nombres alternativos
    if (!is.null(object$fit) &&
        all(c("c", "alpha", "beta") %in% names(object$fit))) {
      return(object$fit[, c("c", "alpha", "beta")])
    }

    stop("Could not find fitted values with columns: c, alpha, beta.")
  }

  fitted_values <- get_fitted(object, method)

  hat_omega <- fitted_values$c
  hat_alpha <- fitted_values$alpha
  hat_beta  <- fitted_values$beta

  ## Si la longitud no coincide, interpolar sobre la grilla de bloques
  if (length(hat_omega) != length(t)) {
    u_fit <- seq(0, 1, length.out = length(hat_omega))

    hat_omega <- approx(u_fit, hat_omega, xout = u, rule = 2)$y
    hat_alpha <- approx(u_fit, hat_alpha, xout = u, rule = 2)$y
    hat_beta  <- approx(u_fit, hat_beta,  xout = u, rule = 2)$y
  }

  ## --------------------------------------------------
  ## 3. Graficar
  ## --------------------------------------------------
  old_mfrow <- par("mfrow")
  old_mar   <- par("mar")
  on.exit(par(mfrow = old_mfrow, mar = old_mar), add = TRUE)
  par(mfrow = mfrow, mar = mar)

  plot_one <- function(y_true, y_hat, g, main_expr,margen) {
    ylim <- c(
      min(c(y_true,y_hat, g), na.rm = TRUE) - margen,
      max(c(y_true,y_hat, g), na.rm = TRUE) + margen
    )

    plot(t_plot, y_true,
         xlim = c(0, max(t_plot)),
         ylim = ylim,
         bty = "n",
         main = main_expr,
         las = 1,
         ylab = "",
         xlab = "",
         pch = pch,
         cex = cex.points,
         cex.axis = cex.axis,
         cex.main = cex.main,
         col = col.points)

      lines(t_plot, y_hat,
          col = col.fit,
          lwd = lwd.fit,
          lty = lty.fit)

    abline(h = g,
           col = col.stationary,
           lwd = lwd.stationary,
           lty = lty.stationary)
  }

  plot_one(c, hat_omega, G[1], expression(c(u)), margen)
  mtext(side=1, text = "u", line=3, cex=1.5)
  plot_one(alpha, hat_alpha, G[2], expression(alpha(u)), margen)
  mtext(side=1, text = "u", line=3, cex=1.5)
  plot_one(beta, hat_beta, G[3], expression(beta(u)), margen)
  mtext(side=1, text = "u", line=3, cex=1.5)

  invisible(data.frame(
    t = t,
    u = u,
    c = c,
    alpha = alpha,
    beta = beta,
    omega_hat = hat_omega,
    alpha_hat = hat_alpha,
    beta_hat = hat_beta
  ))
}

#########################
#### PREDICT HELPERS ####
#########################
tvGarchKalmanPredict <- function(object, n.ahead = 10, trace.log = FALSE)
{
  series0  <- object@data
  nsample0 <- length(series0)

  series <- c(series0, rep(NA, n.ahead))
  nsample <- length(series)

  ## ===============================
  ## 1. Construir c(u), alpha(u), beta(u)
  ## ===============================

  if (inherits(object, "tvGarchKF")) {
    params <- tvCoef(
      x = as.numeric(object@coef),
      data = series,        # use data, no series
      c = object@pars$c,
      alpha = object@pars$alpha,
      beta = object@pars$beta,
      type = object@type,
      exponentes = object@exponentes,
      trig = object@trig,
      arg = object@arg
    )

    c <- as.numeric(params$c)
    a <- as.numeric(params$a)
    b <- as.numeric(params$b)

  } else if (inherits(object, "tvGarchKFNonParam_Fit")) {
    grid_in <- seq_len(nsample0) / nsample0

    fit_omega <- object@Param_model$c
    fit_alpha <- object@Param_model$alpha
    fit_beta  <- object@Param_model$beta

    if (identical(object@method, "Local Polynomial")) {

      c_in <- predict(fit_omega, newdata = data.frame(u = grid_in), type = "response")
      a_in <- predict(fit_alpha, newdata = data.frame(u = grid_in), type = "response")
      b_in <- predict(fit_beta,  newdata = data.frame(u = grid_in), type = "response")

    } else if (identical(object@method, "Smoothing Spline")) {

      c_in <- predict(fit_omega, grid_in)$y
      a_in <- predict(fit_alpha, grid_in)$y
      b_in <- predict(fit_beta,  grid_in)$y

    } else if (identical(object@method, "B-Spline")) {

      c_in <- predict(fit_omega, newdata = data.frame(u = grid_in), type = "response")
      a_in <- predict(fit_alpha, newdata = data.frame(u = grid_in), type = "response")
      b_in <- predict(fit_beta,  newdata = data.frame(u = grid_in), type = "response")

    } else {
      stop("Parametric method not recognized.")
    }

    c <- c(as.numeric(c_in), rep(tail(as.numeric(c_in), 1), n.ahead))
    a <- c(as.numeric(a_in), rep(tail(as.numeric(a_in), 1), n.ahead))
    b <- c(as.numeric(b_in), rep(tail(as.numeric(b_in), 1), n.ahead))

  } else {
    stop("The object must inherit from 'tvGarchKF' or 'tvGarchKFNonParam'.")
  }

  ## ===============================
  ## 2. Correccion solo para prediccion
  ## ===============================

  c <- pmax(c, .Machine$double.eps)
  a <- pmax(a, 0)
  b <- pmax(b, 0)

  ab <- a + b
  idx <- ab >= 0.999

  if (any(idx, na.rm = TRUE)) {
    #warning("Se corrigieron valores con alpha(u) + beta(u) >= 1.")
    a[idx] <- 0.999 * a[idx] / ab[idx]
    b[idx] <- 0.999 * b[idx] / ab[idx]
  }

  if (trace.log) {
    for (i in seq_len(nsample)) {
      cat("i:", i,
          "| c:", c[i],
          "| alpha:", a[i],
          "| beta:", b[i],
          "| alpha+beta:", a[i] + b[i], "\n")
    }
  }

  ## ===============================
  ## 3. Filtro de Kalman
  ## ===============================

  y <- numeric(nsample)

  for (i in seq_len(nsample)) {
    y[i] <- ifelse(
      is.na(series[i]),
      100000000,
      series[i]^2 - c[i] / (1 - a[i] - b[i])
    )
  }

  P   <- rep(0, nsample)
  X   <- rep(0, nsample)
  F.m <- rep(0, nsample)
  K   <- rep(0, nsample)

  denom <- 1 - (a[1] + b[1])^2
  denom <- max(denom, .Machine$double.eps)

  P[1] <- a[1]^2 / denom

  resultado <- .C(
    "tvGarch1ab2",
    as.integer(nsample), as.double(a), as.double(b),
    as.double(c), as.double(P), as.double(X),
    as.double(F.m), as.double(K), as.double(y)
  )

  P   <- resultado[[5]]
  X   <- resultado[[6]]
  F.m <- resultado[[7]]
  K   <- resultado[[8]]

  sigma2 <- as.numeric(X + c / (1 - a - b))
  sigma2 <- pmax(sigma2, .Machine$double.eps)

  list(
    sigma = sqrt(sigma2),
    sigma2 = sigma2,
    P = P,
    X = X,
    Fm = F.m,
    c = c,
    alpha = a,
    beta = b
  )
}

#' @export
print.tvGarchKF_predict <- function(x, digits = 6, max.print = 10, ...) {

  n.ahead <- x$n.ahead
  n <- x$nsample

  sigma_h <- tail(x$sigma_forecast, n.ahead)

  if (!is.null(x$interval)) {

    result <- data.frame(
      `h-step` = seq_len(n.ahead),
      sigma = sigma_h,
      lower = x$interval$lower,
      upper = x$interval$upper,
      check.names = FALSE
    )

  } else {

    result <- data.frame(
      `h-step` = seq_len(n.ahead),
      sigma = sigma_h,
      check.names = FALSE
    )
  }

  cat(
    "\nForecast of volatility sigma (t+",
    n,
    " ..., t+",
    n + n.ahead,
    "):\n",
    sep = ""
  )

  if (!is.null(x$conf)) {
    cat(
      "Confidence level: ",
      round(100 * x$conf, 1),
      "%\n\n",
      sep = ""
    )
  }

  if (nrow(result) > max.print) {

    n_head <- floor(max.print / 2)
    n_tail <- max.print - n_head

    result_print <- rbind(
      head(result, n_head),
      tail(result, n_tail)
    )

    result_print <- round(result_print, digits = digits)

    print(
      result_print[seq_len(n_head), ],
      row.names = FALSE
    )

    cat("      .        .        .        .\n")
    cat("      .        .        .        .\n")
    cat("      .        .        .        .\n")

    print(
      result_print[(n_head + 1):nrow(result_print), ],
      row.names = FALSE
    )

    cat(
      "\n[",
      nrow(result) - max.print,
      " forecast horizons omitted]\n",
      sep = ""
    )

  } else {

    print(
      round(result, digits = digits),
      row.names = FALSE
    )
  }
  invisible(x)
}

#' @method plot tvGarchKF_predict
#' @export
plot.tvGarchKF_predict <- function(x,
                                   main = "Prediction with confidence intervals",
                                   xlab = "Time",
                                   ylab = expression(hat(sigma)[t+h]),
                                   col.sigma = "black",
                                   col.ic = "gray60",
                                   lwd.sigma = 3,
                                   lwd.ic = 2,
                                   lty.sigma = 1,
                                   lty.ic = 2,
                                   legend = TRUE,
                                   legend.pos = "topleft",
                                   legend.cex = 1.2,
                                   cex.axis = 2,
                                   cex.lab = 2,
                                   cex.main = 2,
                                   mar = c(5, 6, 4, 2),
                                   ...)
{
  n <- x$nsample
  n.ahead <- x$n.ahead

  time_future <- n + seq_len(n.ahead)
  sigma_h <- tail(x$sigma_forecast, n.ahead)

  if (!is.null(x$interval)) {
    lower <- x$interval$lower
    upper <- x$interval$upper
  } else {
    lower <- rep(NA_real_, n.ahead)
    upper <- rep(NA_real_, n.ahead)
  }

  vals <- c(sigma_h, lower, upper)
  vals <- vals[is.finite(vals)]

  if (length(vals) == 0)
    stop("No finite values available for plotting.")

  ylim <- range(vals)

  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)

  par(mar = mar, xpd = TRUE)

  plot(
    time_future,
    sigma_h,
    type = "l",
    col = col.sigma,
    lwd = lwd.sigma,
    lty = lty.sigma,
    ylim = ylim,
    xlab = "",
    ylab = "",
    # xaxt = "n",
    # yaxt = "n",
    main = main,
    bty = "n",
    las = 1,
    cex.axis = cex.axis,
    cex.lab = cex.lab,
    cex.main = cex.main,
    ...
  )

  mtext(side = 1, text = xlab, line = 3, cex = 2)
  mtext(side = 2, text = ylab, line = 3.5, cex = 2)

  if (!all(is.na(lower)) && !all(is.na(upper))) {
    lines(time_future, lower, col = col.ic, lwd = lwd.ic, lty = lty.ic)
    lines(time_future, upper, col = col.ic, lwd = lwd.ic, lty = lty.ic)
  }

  if (isTRUE(legend)) {

    if (identical(legend.pos, "right")) {
      usr <- par("usr")
      legend.x <- usr[2] + 0.05 * diff(usr[1:2])
      legend.y <- usr[4]

      legend(
        x = legend.x,
        y = legend.y,
        legend = c(
          expression(hat(sigma)[t+h]),
          paste0(round(100 * x$conf), "% CI")
        ),
        col = c(col.sigma, col.ic),
        lwd = c(lwd.sigma, lwd.ic),
        lty = c(lty.sigma, lty.ic),
        bty = "n",
        cex = legend.cex
      )
    } else {
      legend(
        legend.pos,
        legend = c(
          expression(hat(sigma)[t+h]),
          paste0(round(100 * x$conf), "% CI")
        ),
        col = c(col.sigma, col.ic),
        lwd = c(lwd.sigma, lwd.ic),
        lty = c(lty.sigma, lty.ic),
        bty = "n",
        cex = legend.cex
      )
    }
  }

  invisible(x)
}
