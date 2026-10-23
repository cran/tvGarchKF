### Smoothing Spline ###
run_smoothing_and_plotting <- function(data, tvparameter, spar = NULL, nknots, methods = NULL,
                                       control.spar = list(lower = 0.5, upper = 1.2), max_knots = NULL, plot = TRUE)
  {
  t  <- tvparameter$t
  w1 <- tvparameter$c
  a1 <- tvparameter$alpha
  b1 <- tvparameter$beta
  n  <- length(data)
  u0 <- t/n

  if (is.null(max_knots)) {
    max_knots <- length(unique(u0))
  }
  if (length(nknots) == 1) {
    nknots <- rep(nknots, 3)
  }
  nknots <- pmin(nknots, max_knots)
  if (length(nknots) != 3) {
    stop("nknots must have length 1 or 3.")
  }
  if (!is.null(spar)) {
    if (length(spar) == 1) spar <- rep(spar, 3)
    if (length(spar) != 3) stop("spar must be NULL, length 1, or length 3.")
  }
  if (!is.null(methods)) {
    if (length(methods) == 1) methods <- rep(methods, 3)
    if (length(methods) != 3) stop("methods must be NULL, length 1, or length 3.")
  }

  estacionario1 <- garchFit(~ garch(1,1),
                            data = data, include.mean = FALSE,
                            cond.dist = "norm", trace = FALSE)

  G1 <- estacionario1@fit$matcoef

  omega.ss <- fit_ss_one(
    x = u0,
    y = w1,
    spar_i = if (!is.null(spar)) spar[1] else NULL,
    method_i = if (!is.null(methods)) methods[1] else NULL,
    nknots_i = nknots[1],
    control.spar = control.spar
  )

  alpha.ss <- fit_ss_one(
    x = u0,
    y = a1,
    spar_i = if (!is.null(spar)) spar[2] else NULL,
    method_i = if (!is.null(methods)) methods[2] else NULL,
    nknots_i = nknots[2],
    control.spar = control.spar
  )

  beta.ss <- fit_ss_one(
    x = u0,
    y = b1,
    spar_i = if (!is.null(spar)) spar[3] else NULL,
    method_i = if (!is.null(methods)) methods[3] else NULL,
    nknots_i = nknots[3],
    control.spar = control.spar
  )

  curva.smooth <- data.frame(c = omega.ss$y,
                             alpha = alpha.ss$y,
                             beta = beta.ss$y)

  new_t <- (1:n)/n
  pred.omega.smooth <- predict(omega.ss, new_t)$y
  pred.alpha.smooth <- predict(alpha.ss, new_t)$y
  pred.beta.smooth <- predict(beta.ss, new_t)$y

  Condition_final <- summary(pred.alpha.smooth + pred.beta.smooth)
  param_tv.smooth <- data.frame(
    t = new_t,
    c = pred.omega.smooth,
    alpha = pred.alpha.smooth,
    beta  = pred.beta.smooth)

  if (plot) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow = c(1, 3), mar = c(5, 6, 4, 2))

    t_plot = (t-min(t))/(max(t)-min(t))

    plot(t_plot, w1, xlim = c(0, max(t_plot)), ylim = c(min(w1) * 0.5, max(w1) * 1.5),
         bty = "n", main = expression(~c(u)), las = 1, lwd = 3,
         ylab = "", cex.axis = 1.2, xlab = "", cex.main = 2, pch = 20)
    lines(t_plot, omega.ss$y, lwd = 3, col = "deeppink", lty = 1)
    abline(h = G1[1, 1], lwd = 2, lty = 2)
    mtext(side=1, text = "u", line=3, cex=1.5)

    plot(t_plot, a1, xlim = c(0, max(t_plot)), ylim = c(min(a1) * 0.5, max(a1) * 1.5),
         bty = "n", main = expression(~alpha(u)), las = 1, lwd = 3,
         ylab = "", cex.axis = 1.2, xlab = "", cex.main = 2, pch = 20)
    lines(t_plot, alpha.ss$y, lwd = 3, col = "deeppink", lty = 1)
    abline(h = G1[2, 1], lwd = 2, lty = 2)
    mtext(side=1, text = "u", line=3, cex=1.5)

    plot(t_plot, b1, xlim = c(0, max(t_plot)), ylim = c(min(b1) * 0.5, max(b1) * 1.5),
         bty = "n", main = expression(~beta(u)), las = 1, lwd = 3,
         ylab = "", cex.axis = 1.2, xlab = "", cex.main = 2, pch = 20)
    lines(t_plot, beta.ss$y, lwd = 3, col = "deeppink", lty = 1)
    abline(h = G1[3, 1], lwd = 2, lty = 2)
    mtext(side=1, text = "u", line=3, cex=1.5)
  }

  result <- list(
    "Check condition" = Condition_final,
    tvparameter = tvparameter,
    fitted.values   = curva.smooth,
    ParamStationary = G1,
    HattvParameters = param_tv.smooth,
    data = data,
    smoothing_used = list(
      spar = spar,
      methods = methods,
      nknots = nknots
    ),
    Param_model = list(
      c = omega.ss,
      alpha = alpha.ss,
      beta  = beta.ss
    )
  )
  class(result) <- c("tvgarchKF", class(result))
  return(result)
}

fit_ss_one <- function(x, y, spar_i = NULL, method_i = NULL, nknots_i, control.spar) {
  args <- list(
    x = x,
    y = y,
    nknots = nknots_i
  )
  if (!is.null(spar_i)) {
    args$spar <- spar_i
  } else {
    if (!is.null(method_i)) args$method <- method_i
    args$control.spar <- control.spar
  }
  do.call(npreg::ss, args)
}

