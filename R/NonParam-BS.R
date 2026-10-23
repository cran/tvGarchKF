#' @importFrom splines bs
#' @importFrom stats glm
##### B-Spline model #####
fit_b_spline_models <- function(data, tvparameter, plot = TRUE,
                                degree, df) {
  n <- length(data)
  estacionario1 <- garchFit(~ garch(1,1),
                            data = data, include.mean = FALSE,
                            cond.dist = "norm", trace = FALSE)
  G1 <- estacionario1@fit$matcoef

  t  <- tvparameter$t
  c <- tvparameter$c
  alpha <- tvparameter$alpha
  beta <- tvparameter$beta

  degree <- degree
  optimal_df <- df
  tvparameter$u <- t/n

  fit_omega.bs <- glm(
    c ~ splines::bs(u, degree = degree[1], df = optimal_df[1], Boundary.knots = c(0, 1)),
    data = tvparameter
  )
  fit_alpha.bs <- glm(
    alpha ~ splines::bs(u, degree = degree[2], df = optimal_df[2], Boundary.knots = c(0, 1)),
    data = tvparameter
  )
  fit_beta.bs <- glm(
    beta ~ splines::bs(u, degree = degree[3], df = optimal_df[3], Boundary.knots = c(0, 1)),
    data = tvparameter
  )
  grid_u <- seq_len(n) / n

  pred.omega.bs <- predict(fit_omega.bs, newdata = data.frame(u = grid_u), type="response")
  pred.alpha.bs <- predict(fit_alpha.bs, newdata = data.frame(u = grid_u), type="response")
  pred.beta.bs  <- predict(fit_beta.bs,  newdata = data.frame(u = grid_u), type="response")

  curva.bs <- data.frame(
    c = fit_omega.bs$fitted.values,
    alpha = fit_alpha.bs$fitted.values,
    beta  = fit_beta.bs$fitted.values
  )

  Condition <- summary(pred.alpha.bs + pred.beta.bs)
  if (all((pred.alpha.bs + pred.beta.bs) == 1)) {
    message("The stationarity condition is not satisfied: try changing the degrees/df grids")
  }
  param_tv.bs <- data.frame(t = seq_len(n),c = pred.omega.bs, alpha = pred.alpha.bs, beta = pred.beta.bs)

  if (isTRUE(plot)) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    x <- tvparameter$t
    par(mfrow = c(1, 3), mar = c(5, 6, 4, 2))

    t_plot = (t-min(t))/(max(t)-min(t))

    plot(t_plot, tvparameter$c, xlim = c(0, max(t_plot)),
         ylim = c(min(tvparameter$c)*0.5, max(tvparameter$c)*1.5),
         bty = "n", main = expression(~c(u)), las = 1, lwd = 3, ylab = "",
         cex.axis = 1.2, xlab = "", cex.main = 2.0, pch = 20)
    lines(t_plot, fit_omega.bs$fitted.values, lwd = 3, col = "blue")
    abline(h = G1[1, 1], lwd = 2, lty = 2)
    mtext(side=1, text = "u", line=3, cex=1.5)

    plot(t_plot, tvparameter$alpha, xlim = c(0, max(t_plot)),
         ylim = c(min(tvparameter$alpha)*0.5, max(tvparameter$alpha)*1.5),
         bty = "n", main = expression(~alpha(u)), las = 1, lwd = 3, ylab = "",
         cex.axis = 1.2, xlab = "", cex.main = 2.0, pch = 20)
    lines(t_plot, fit_alpha.bs$fitted.values, lwd = 3, col = "blue")
    abline(h = G1[2, 1], lwd = 2, lty = 2)
    mtext(side=1, text = "u", line=3, cex=1.5)

    plot(t_plot, tvparameter$beta, xlim = c(0, max(t_plot)),
         ylim = c(min(tvparameter$beta)*0.5, max(tvparameter$beta)*1.5),
         bty = "n", main = expression(~beta(u)), las = 1, lwd = 3, ylab = "",
         cex.axis = 1.2, xlab = "", cex.main = 2.0, pch = 20)
    lines(t_plot, fit_beta.bs$fitted.values, lwd = 3, col = "blue")
    abline(h = G1[3, 1], lwd = 2, lty = 2)
    mtext(side=1, text = "u", line=3, cex=1.5)
  }

  out <- list(
    method            = "B-Spline",
    "Check condition" = Condition,
    tvparameter       = tvparameter,
    fitted.values     = curva.bs,
    ParamStationary   = G1,
    HattvParameters   = param_tv.bs,
    data              = data,
    Param_model       = list(c = fit_omega.bs, alpha = fit_alpha.bs, beta = fit_beta.bs),
    best_degree       = degree,
    best_df           = optimal_df
  )
  class(out) <- c("tvgarchKF", class(out))
  return(out)
}
