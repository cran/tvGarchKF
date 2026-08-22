#' @importFrom mgcv gam
#' @importFrom stats loess
#' @importFrom utils capture.output
#' @importFrom stats AIC
#' @importFrom stats gaussian
select_best_method_tv_stationary <- function(object,
                                             nknots_candidates = c(4, 5, 6, 8, 10),
                                             methods = c("REML"),
                                             metric = c("mse", "mae"),
                                             max_knots = NULL,
                                             tol = 1e-3,
                                             bs = "ps",
                                             lambda_r = 1,
                                             max_slope = Inf,
                                             max_curvature = Inf) {

  metric <- match.arg(metric)

  tvparameter <- as.data.frame(object@coef)

  x     <- tvparameter$t
  c <- tvparameter$c
  alpha <- tvparameter$alpha
  beta  <- tvparameter$beta

  data_vec <- object@data
  if (is.null(data_vec)) stop("object@data is missing.")

  n <- length(data_vec)
  u0 <- x / n
  pred_grid <- seq(min(u0), max(u0), length.out = n)

  uniq_x <- length(unique(u0))
  max_eff <- if (is.null(max_knots)) uniq_x else min(max_knots, uniq_x)

  k_candidates <- unique(pmin(nknots_candidates, max_eff))
  k_candidates <- k_candidates[k_candidates >= 3]

  if (!length(k_candidates)) {
    stop("No valid k candidates available.")
  }

  gam_method <- methods[1]

  roughness_fun <- function(u, yhat) {
    ord <- order(u)
    u <- u[ord]
    yhat <- yhat[ord]

    du <- diff(u)
    if (any(du <= 0)) return(Inf)

    dy <- diff(yhat) / du
    d2 <- diff(dy) / diff(u[-1])

    mean(d2^2, na.rm = TRUE)
  }

  smoothness_fun <- function(u, yhat) {
    ord <- order(u)
    u <- u[ord]
    yhat <- yhat[ord]

    du <- diff(u)
    if (any(du <= 0)) {
      return(list(ok = FALSE, max_slope = Inf, max_curvature = Inf))
    }

    dy <- diff(yhat) / du
    d2 <- diff(dy) / diff(u[-1])

    ms <- max(abs(dy), na.rm = TRUE)
    mc <- max(abs(d2), na.rm = TRUE)

    list(
      ok = ms <= max_slope && mc <= max_curvature,
      max_slope = ms,
      max_curvature = mc
    )
  }

  score_fun <- function(y, yhat, u) {

    fit_error <- if (metric == "mse") {
      mean((y - yhat)^2, na.rm = TRUE)
    } else {
      mean(abs(y - yhat), na.rm = TRUE)
    }

    rough <- roughness_fun(u, yhat)

    fit_error + lambda_r * rough
  }

  fit_gam_one <- function(y, k) {
    dat <- data.frame(u = u0, y = y)

    try(
      mgcv::gam(
        y ~ s(u, bs = bs, k = k),
        data = dat,
        method = gam_method
      ),
      silent = TRUE
    )
  }

  extract_info <- function(fit) {
    list(
      sp = if (!is.null(fit$sp)) as.numeric(fit$sp[1]) else NA_real_,
      edf = sum(fit$edf, na.rm = TRUE),
      fitted = as.numeric(fitted(fit))
    )
  }

  # =========================================================
  # Omega: selección independiente con suavidad
  # =========================================================

  best_omega_score <- Inf
  best_omega <- NULL

  for (k_omega in k_candidates) {

    fit_omega <- fit_gam_one(c, k_omega)
    if (inherits(fit_omega, "try-error") || is.null(fit_omega)) next

    info_omega <- extract_info(fit_omega)

    omega_pred <- as.numeric(
      predict(fit_omega, newdata = data.frame(u = pred_grid))
    )

    if (length(omega_pred) != length(pred_grid) || anyNA(omega_pred)) next

    sm_omega <- smoothness_fun(pred_grid, omega_pred)
    if (!isTRUE(sm_omega$ok)) next

    err_omega <- score_fun(c, info_omega$fitted, u0)

    if (is.finite(err_omega) && err_omega < best_omega_score) {

      best_omega_score <- err_omega

      best_omega <- list(
        method = gam_method,
        k = k_omega,
        nknots = k_omega,
        sp = info_omega$sp,
        spar = info_omega$sp,
        df = info_omega$edf,
        edf = info_omega$edf,
        fitted = info_omega$fitted,
        predicted = omega_pred,
        smoothness = sm_omega,
        score = err_omega,
        fit = fit_omega
      )
    }
  }

  if (is.null(best_omega)) {
    stop("No valid smooth GAM fit was found for c")
  }

  # =========================================================
  # Alpha + Beta: selección conjunta
  # =========================================================

  best_score <- Inf
  best_alpha <- NULL
  best_beta  <- NULL
  best_stat  <- NULL
  best_sm_alpha <- NULL
  best_sm_beta  <- NULL

  for (k_alpha in k_candidates) {

    fit_alpha <- fit_gam_one(alpha, k_alpha)
    if (inherits(fit_alpha, "try-error") || is.null(fit_alpha)) next

    info_alpha <- extract_info(fit_alpha)
    alpha_fit <- info_alpha$fitted

    alpha_pred <- as.numeric(
      predict(fit_alpha, newdata = data.frame(u = pred_grid))
    )

    if (length(alpha_pred) != length(pred_grid) || anyNA(alpha_pred)) next

    sm_alpha <- smoothness_fun(pred_grid, alpha_pred)
    if (!isTRUE(sm_alpha$ok)) next

    for (k_beta in k_candidates) {

      fit_beta <- fit_gam_one(beta, k_beta)
      if (inherits(fit_beta, "try-error") || is.null(fit_beta)) next

      info_beta <- extract_info(fit_beta)
      beta_fit <- info_beta$fitted

      beta_pred <- as.numeric(
        predict(fit_beta, newdata = data.frame(u = pred_grid))
      )

      if (length(beta_pred) != length(pred_grid) || anyNA(beta_pred)) next

      sm_beta <- smoothness_fun(pred_grid, beta_pred)
      if (!isTRUE(sm_beta$ok)) next

      stat <- .stationarity_metrics(alpha_pred, beta_pred, tol = tol)
      if (!isTRUE(stat$ok)) {
        message(
          "Rejected: k_alpha=", k_alpha,
          ", k_beta=", k_beta,
          ", max(alpha+beta)=",
          round(stat$max_alpha_beta, 6)
        )
        next
      }

      score_alpha <- score_fun(alpha, alpha_fit, u0)
      score_beta  <- score_fun(beta, beta_fit, u0)

      joint_score <- mean(c(score_alpha, score_beta), na.rm = TRUE)

      if (is.finite(joint_score) && joint_score < best_score) {

        best_score <- joint_score

        best_alpha <- list(
          method = gam_method,
          k = k_alpha,
          nknots = k_alpha,
          sp = info_alpha$sp,
          spar = info_alpha$sp,
          df = info_alpha$edf,
          edf = info_alpha$edf,
          fitted = alpha_fit,
          predicted = alpha_pred,
          smoothness = sm_alpha,
          score = score_alpha,
          fit = fit_alpha
        )

        best_beta <- list(
          method = gam_method,
          k = k_beta,
          nknots = k_beta,
          sp = info_beta$sp,
          spar = info_beta$sp,
          df = info_beta$edf,
          edf = info_beta$edf,
          fitted = beta_fit,
          predicted = beta_pred,
          smoothness = sm_beta,
          score = score_beta,
          fit = fit_beta
        )

        best_stat <- stat
        best_sm_alpha <- sm_alpha
        best_sm_beta  <- sm_beta
      }
    }
  }

  if (is.null(best_alpha) || is.null(best_beta)) {
    stop("No stationary and smooth alpha-beta GAM combination was found.")
  }

  summary <- data.frame(
    parameter = c("c", "alpha", "beta"),
    method = c(best_omega$method, best_alpha$method, best_beta$method),
    nknots = c(best_omega$nknots, best_alpha$nknots, best_beta$nknots),
    k = c(best_omega$k, best_alpha$k, best_beta$k),
    spar = c(best_omega$spar, best_alpha$spar, best_beta$spar),
    sp = c(best_omega$sp, best_alpha$sp, best_beta$sp),
    df = c(best_omega$df, best_alpha$df, best_beta$df),
    score = c(best_omega$score, best_alpha$score, best_beta$score),
    max_slope = c(
      best_omega$smoothness$max_slope,
      best_alpha$smoothness$max_slope,
      best_beta$smoothness$max_slope
    ),
    max_curvature = c(
      best_omega$smoothness$max_curvature,
      best_alpha$smoothness$max_curvature,
      best_beta$smoothness$max_curvature
    )
  )

  list(
    summary = summary,
    stationarity = best_stat,

    smoothness = list(
      c = best_omega$smoothness,
      alpha = best_sm_alpha,
      beta  = best_sm_beta
    ),

    fits = list(
      c = best_omega,
      alpha = best_alpha,
      beta  = best_beta
    ),

    fitted = data.frame(
      t = x,
      c = best_omega$fitted,
      alpha = best_alpha$fitted,
      beta  = best_beta$fitted,
      alpha_beta = best_alpha$fitted + best_beta$fitted
    ),

    predicted = data.frame(
      t = pred_grid,
      c = best_omega$predicted,
      alpha = best_alpha$predicted,
      beta  = best_beta$predicted,
      alpha_beta = best_alpha$predicted + best_beta$predicted
    )
  )
}


# métricas de estacionariedad
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


select_best_k_gam <- function(x, y,
                              k_grid = 3:8,
                              bs = "ps",
                              method = "REML") {

  out <- lapply(k_grid, function(K) {

    fit <- try(
      mgcv::gam(
        y ~ s(x, bs = bs, k = K),
        method = method
      ),
      silent = TRUE
    )

    if (inherits(fit, "try-error")) return(NULL)

    data.frame(
      k = K,
      sp = fit$sp[1],
      edf = sum(fit$edf),
      AIC = AIC(fit),
      GCV = fit$gcv.ubre,
      REML_score = fit$gcv.ubre
    )
  })

  tab <- do.call(rbind, out)

  tab[which.min(tab$AIC), ]
}

select_best_k_by_edf <- function(x, y,
                                 k_grid = 3:8,
                                 bs = "ps",
                                 method = "REML",
                                 edf_min = 2,
                                 edf_max = sqrt(length(y))) {

  res <- lapply(k_grid, function(K) {

    fit <- try(
      mgcv::gam(
        y ~ s(x, bs = bs, k = K),
        method = method
      ),
      silent = TRUE
    )

    if (inherits(fit, "try-error")) return(NULL)

    data.frame(
      k = K,
      sp = as.numeric(fit$sp[1]),
      edf = sum(fit$edf),
      edf_target = (edf_min + edf_max) / 2,
      fit = I(list(fit))
    )
  })

  tab <- do.call(rbind, res)

  tab$edf_distance <- abs(tab$edf - tab$edf_target)
  tab$valid_edf <- tab$edf >= edf_min & tab$edf <= edf_max

  if (any(tab$valid_edf)) {
    best <- tab[tab$valid_edf, ]
    best <- best[which.min(best$edf_distance), ]
  } else {
    best <- tab[which.min(tab$edf_distance), ]
  }

  list(
    best_k = best$k,
    best_sp = best$sp,
    best_edf = best$edf,
    edf_table = tab[, c("k", "sp", "edf", "edf_distance", "valid_edf")],
    fit = best$fit[[1]]
  )
}


match_ss_spar_by_edf <- function(x, y,
                                 edf_target,
                                 nknots,
                                 spar_grid = seq(0.3, 1.3, by = 0.01),
                                 method = NULL,
                                 control.spar = list(lower = 0.3, upper = 1.3)) {

  stopifnot(length(x) == length(y))

  res <- lapply(spar_grid, function(sp) {

    fit <- try(
      npreg::ss(
        x = x,
        y = y,
        spar = sp,
        nknots = nknots
      ),
      silent = TRUE
    )

    if (inherits(fit, "try-error") || is.null(fit)) {
      return(NULL)
    }

    info <- .pluck_ss_info(fit)

    if (is.null(info$df) || !is.finite(info$df)) {
      return(NULL)
    }

    data.frame(
      spar = sp,
      nknots = nknots,
      edf_ss = info$df,
      edf_target = edf_target,
      edf_distance = abs(info$df - edf_target),
      fit = I(list(fit))
    )
  })

  tab <- do.call(rbind, res)

  if (is.null(tab) || nrow(tab) == 0) {
    stop("No valid npreg::ss fit was obtained for the supplied spar_grid.")
  }

  best <- tab[which.min(tab$edf_distance), ]

  list(
    best_spar = best$spar,
    best_nknots = best$nknots,
    best_edf_ss = best$edf_ss,
    edf_target = best$edf_target,
    edf_distance = best$edf_distance,
    fit = best$fit[[1]],
    grid = tab[, c("spar", "nknots", "edf_ss", "edf_target", "edf_distance")]
  )
}

.pluck_ss_info <- function(fit) {

  out <- list(
    yhat = NULL,
    spar = NA_real_,
    lambda = NA_real_,
    df = NA_real_,
    knots = NULL
  )

  if (is.null(fit)) return(out)

  # Valores ajustados
  if (!is.null(fit$y)) {
    out$yhat <- as.numeric(fit$y)
  } else if (!is.null(fit$fitted.values)) {
    out$yhat <- as.numeric(fit$fitted.values)
  } else if (!is.null(fit$fitted)) {
    out$yhat <- as.numeric(fit$fitted)
  }

  # spar
  if (!is.null(fit$spar)) {
    out$spar <- as.numeric(fit$spar)
  }

  # lambda
  if (!is.null(fit$lambda)) {
    out$lambda <- as.numeric(fit$lambda)
  }

  # df / edf
  if (!is.null(fit$df)) {
    out$df <- as.numeric(fit$df)
  } else if (!is.null(fit$edf)) {
    out$df <- as.numeric(fit$edf)
  } else {
    sm <- try(summary(fit), silent = TRUE)

    if (!inherits(sm, "try-error")) {
      # Intentar extraer df desde summary
      txt <- capture.output(sm)
      df_line <- grep("Residual standard error", txt, value = TRUE)

      if (length(df_line) > 0) {
        tmp <- regmatches(df_line, regexpr("[0-9]+\\.?[0-9]* degrees of freedom", df_line))
        if (length(tmp) > 0) {
          out$df <- as.numeric(sub(" degrees of freedom", "", tmp))
        }
      }
    }
  }

  # knots
  if (!is.null(fit$knots)) {
    out$knots <- fit$knots
  }

  out
}

match_ss_spar_by_edf_tv <- function(tvparameter,
                                    best_k_omega,
                                    best_k_alpha,
                                    best_k_beta,
                                    data_length,
                                    spar_grid = seq(0.3, 1.3, by = 0.01)) {

  # tvparameter <- .normalize_tv_names(as.data.frame(tvparameter))

  u <- tvparameter$t / data_length

  best_list <- list(
    c = best_k_omega,
    alpha = best_k_alpha,
    beta  = best_k_beta
  )

  res <- lapply(names(best_list), function(par) {

    best_gam <- best_list[[par]]

    fit_ss <- match_ss_spar_by_edf(
      x = u,
      y = tvparameter[[par]],
      edf_target = best_gam$edf,
      nknots = best_gam$k,
      spar_grid = spar_grid
    )

    data.frame(
      parameter = par,
      gam_k = best_gam$k,
      gam_edf = best_gam$edf,
      ss_spar = fit_ss$best_spar,
      ss_nknots = fit_ss$best_nknots,
      ss_edf = fit_ss$best_edf_ss,
      edf_distance = fit_ss$edf_distance
    )
  })

  do.call(rbind, res)
}

match_loess_span_by_edf <- function(x, y,
                                    edf_target,
                                    degree = 1,
                                    span_grid = seq(0.50, 2.50, by = 0.05)) {

  stopifnot(length(x) == length(y))

  res <- lapply(span_grid, function(spn) {

    fit <- try(
      suppressWarnings(
        loess(
          y ~ x,
          span = spn,
          degree = degree,
          control = loess.control(
            surface = "direct",
            trace.hat = "exact"
          )
        )
      ),
      silent = TRUE
    )

    if (inherits(fit, "try-error") || is.null(fit)) return(NULL)

    edf_loess <- fit$trace.hat

    if (is.null(edf_loess) || !is.finite(edf_loess)) return(NULL)

    data.frame(
      span = spn,
      degree = degree,
      edf_loess = edf_loess,
      edf_target = edf_target,
      edf_distance = abs(edf_loess - edf_target),
      fit = I(list(fit))
    )
  })

  tab <- do.call(rbind, res)

  if (is.null(tab) || nrow(tab) == 0) {
    stop("No valid loess fit was obtained.")
  }

  best <- tab[which.min(tab$edf_distance), ]

  list(
    best_span = best$span,
    best_degree = best$degree,
    best_edf_loess = best$edf_loess,
    edf_target = best$edf_target,
    edf_distance = best$edf_distance,
    fit = best$fit[[1]],
    grid = tab[, c("span", "degree", "edf_loess", "edf_target", "edf_distance")]
  )
}


match_loess_span_by_edf_tv <- function(tvparameter,
                                       best_k_omega,
                                       best_k_alpha,
                                       best_k_beta,
                                       data_length,
                                       degree = c(1, 2, 2),
                                       span_grid = seq(0.50, 2.50, by = 0.05)) {

  tvparameter <- .normalize_tv_names(as.data.frame(tvparameter))

  u <- tvparameter$t / data_length

  if (length(degree) == 1) degree <- rep(degree, 3)

  best_list <- list(
    c = best_k_omega,
    alpha = best_k_alpha,
    beta  = best_k_beta
  )

  deg_list <- c(
    c = degree[1],
    alpha = degree[2],
    beta  = degree[3]
  )

  res <- lapply(names(best_list), function(par) {

    best_gam <- best_list[[par]]

    fit_loess <- match_loess_span_by_edf(
      x = u,
      y = tvparameter[[par]],
      edf_target = best_gam$edf,
      degree = deg_list[par],
      span_grid = span_grid
    )

    data.frame(
      parameter = par,
      gam_k = best_gam$k,
      gam_edf = best_gam$edf,
      loess_span = fit_loess$best_span,
      loess_degree = fit_loess$best_degree,
      loess_edf = fit_loess$best_edf_loess,
      edf_distance = fit_loess$edf_distance
    )
  })

  do.call(rbind, res)
}

##### Bs-Spline

match_bs_by_gam_edf <- function(x, y,
                                edf_target,
                                degree_grid = 1:3,
                                df_grid = 3:8,
                                family = gaussian()) {

  stopifnot(length(x) == length(y))

  res <- list()

  for (deg in degree_grid) {

    df_ok <- df_grid[df_grid > deg & df_grid <= length(unique(x))]

    for (df in df_ok) {

      fit <- try(
        suppressWarnings(
          glm(
            y ~ splines::bs(
              x,
              degree = deg,
              df = df,
              Boundary.knots = range(x, finite = TRUE)
            ),
            data = data.frame(x = x, y = y),
            family = family
          )
        ),
        silent = TRUE
      )

      if (inherits(fit, "try-error") || is.null(fit)) next

      yhat <- as.numeric(fitted(fit))

      # aproximación simple del EDF: número de coeficientes efectivos
      edf_bs <- fit$rank

      res[[length(res) + 1L]] <- data.frame(
        degree = deg,
        df = df,
        edf_bs = edf_bs,
        edf_target = edf_target,
        edf_distance = abs(edf_bs - edf_target),
        fit = I(list(fit))
      )
    }
  }

  tab <- do.call(rbind, res)

  if (is.null(tab) || nrow(tab) == 0) {
    stop("No valid B-spline fit was obtained.")
  }

  best <- tab[which.min(tab$edf_distance), ]

  list(
    best_degree = best$degree,
    best_df = best$df,
    best_edf_bs = best$edf_bs,
    edf_target = best$edf_target,
    edf_distance = best$edf_distance,
    fit = best$fit[[1]],
    grid = tab[, c("degree", "df", "edf_bs", "edf_target", "edf_distance")]
  )
}

match_bs_by_gam_edf_tv <- function(tvparameter,
                                   best_k_omega,
                                   best_k_alpha,
                                   best_k_beta,
                                   data_length,
                                   degree_grid = 1:3,
                                   df_grid = 3:8) {

  tvparameter <- .normalize_tv_names(as.data.frame(tvparameter))

  u <- tvparameter$t / data_length

  best_list <- list(
    c = best_k_omega,
    alpha = best_k_alpha,
    beta  = best_k_beta
  )

  res <- lapply(names(best_list), function(par) {

    best_gam <- best_list[[par]]

    fit_bs <- match_bs_by_gam_edf(
      x = u,
      y = tvparameter[[par]],
      edf_target = best_gam$edf,
      degree_grid = degree_grid,
      df_grid = df_grid
    )

    data.frame(
      parameter = par,
      gam_k = best_gam$k,
      gam_edf = best_gam$edf,
      bs_degree = fit_bs$best_degree,
      bs_df = fit_bs$best_df,
      bs_edf = fit_bs$best_edf_bs,
      edf_distance = fit_bs$edf_distance
    )
  })

  do.call(rbind, res)
}
