#' @importFrom stats loess
#' @importFrom stats loess.control
### Local Polynomial
run_polynomial_and_local_fit <- function(data, tvparameter, plot = TRUE,
                                         span, degree) {
  n <- length(data)
  estacionario1 <- garchFit(~ garch(1,1),
                            data = data, include.mean = FALSE,
                            cond.dist = "norm", trace = FALSE)
  G1 <- estacionario1@fit$matcoef

  t  <- tvparameter$t
  w1 <- tvparameter$c
  a1 <- tvparameter$alpha
  b1 <- tvparameter$beta

  tvparameter$u <- t/n

  fit_omega <- loess(
    c ~ u,
    data = tvparameter,
    span = span[1],
    degree = degree[1],
    control = loess.control(surface = "direct")
  )
  fit_alpha <- loess(
    alpha ~ u,
    data = tvparameter,
    span = span[2],
    degree = degree[2],
    control = loess.control(surface = "direct")
  )
  fit_beta <- loess(
    beta ~ u,
    data = tvparameter,
    span = span[3],
    degree = degree[3],
    control = loess.control(surface = "direct")
  )

  grid_t <- seq(min(tvparameter$u), max(tvparameter$u), length.out = n)

  pred.omega.lp <- try(
    suppressWarnings(
      as.numeric(predict(fit_omega, newdata = data.frame(u = grid_t)))
    ),
    silent = TRUE
  )
  pred.alpha.lp <- try(
    suppressWarnings(
      as.numeric(predict(fit_alpha, newdata = data.frame(u = grid_t)))
    ),
    silent = TRUE
  )
  pred.beta.lp <- try(
    suppressWarnings(
      as.numeric(predict(fit_beta, newdata = data.frame(u = grid_t)))
    ),
    silent = TRUE
  )

  alpha_beta_sum <- pred.alpha.lp + pred.beta.lp
  max_ab   <- suppressWarnings(max(alpha_beta_sum, na.rm = TRUE))
  Condition <- summary(alpha_beta_sum)
  if (all((alpha_beta_sum) == 1)) {
    message("The stationarity condition is not satisfied: try changing h_frac/degree_grid.")
  }

  curva.smooth <- data.frame(c = fitted(fit_omega),
                             alpha = fitted(fit_alpha),
                             beta = fitted(fit_beta))

  param_tv.loc <- data.frame(
    t     = grid_t,
    c = pred.omega.lp,
    alpha = pred.alpha.lp,
    beta  = pred.beta.lp
  )

  local_used = list(
    span = span,
    degree = degree
  )

  if (isTRUE(plot)) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow = c(1, 3), mar = c(5, 6, 4, 2))

    t_plot = (t-min(t))/(max(t)-min(t))

    plot(t_plot, w1, xlim = c(0, max(t_plot)), ylim = c(min(w1)*0.5, max(w1)*1.5),
         bty = "n", main = expression(~c(u)), las = 1, lwd = 3,
         ylab = "", cex.axis = 1.2, xlab = "", cex.main = 2, pch = 20)
    lines(t_plot, fit_omega$fitted, lwd = 3, col = "orange")
    abline(h = G1[1, 1], lwd = 2, lty = 2)
    mtext(side=1, text = "u", line=3, cex=1.5)

    plot(t_plot, a1, xlim = c(0, max(t_plot)), ylim = c(min(a1)*0.5, max(a1)*1.5),
         bty = "n", main = expression(~alpha(u)), las = 1, lwd = 3,
         ylab = "", cex.axis = 1.2, xlab = "", cex.main = 2, pch = 20)
    lines(t_plot, fit_alpha$fitted, lwd = 3, col = "orange")
    abline(h = G1[2, 1], lwd = 2, lty = 2)
    mtext(side=1, text = "u", line=3, cex=1.5)

    plot(t_plot, b1, xlim = c(0, max(t_plot)), ylim = c(min(b1)*0.5, max(b1)*1.5),
         bty = "n", main = expression(~beta(u)), las = 1, lwd = 3,
         ylab = "", cex.axis = 1.2, xlab = "", cex.main = 2, pch = 20)
    lines(t_plot, fit_beta$fitted, lwd = 3, col = "orange")
    abline(h = G1[3, 1], lwd = 2, lty = 2)
    mtext(side=1, text = "u", line=3, cex=1.5)
  }

  out <- list(
    "Check condition" = Condition,
    tvparameter       = tvparameter,
    fitted.values     = curva.smooth,
    HattvParameters   = param_tv.loc,
    ParamStationary   = G1,
    data              = data,
    Param_model       = list(c = fit_omega, alpha = fit_alpha, beta = fit_beta),
    local_used        = list(span = span, degree = degree)
  )
  class(out) <- c("tvgarchKF")
  return(out)
}
