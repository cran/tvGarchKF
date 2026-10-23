#' @title Estimate the local structure of tv-GARCH(1,1) parameters.
#' @description This function performs an exploratory analysis of the dynamic structure of the time-varying GARCH(1,1) parameters.
#' Specifically, the observation domain \eqn{\{1, \ldots,T\}} is partitioned into \eqn{M} overlapping blocks, each of length \eqn{N}, with a
#' constant shift of size \eqn{S} between consecutive blocks. These quantities satisfy \eqn{T = S(M-1)+N}. The midpoint
#' of the \eqn{j}-th block, for \eqn{j=1,\ldots,M}, is defined as \eqn{t_j=S(j-1)+N/2}. For each block, a stationary GARCH(1,1) model
#' is fitted using the obervations contained in that block. The resulting sequence of local estimates provides an
#' empirical trajectory of the time-varying parameters \eqn{(c(u))}, \eqn{\alpha(u)} and \eqn{\beta(u)}. These trajectories can
#' can be used as an exploratory tool for identifying suitable functional forms or as initial estimates for subsequent parametric
#' or non-parametric estimation of the tv-GARCH parameter curves.
#' @param data A univariate time series, where could be an object of class "ts" or "zoo".
#' @param S A positive integer specifying the step size or shift between consecutive overlapping blocks.
#' @param N A positive integer specifying the number of observations contained in each block.
#' @param plot A logical value indicating whether the local estimates of the time-varying parameter curves should be displayed. The default is TRUE.
#' @param title An optional character string specifying a title associated with the fitted object. The default is NULL.
#' @param description An optional character string providing a user-defined description. The default is NULL.
#' @param cex.axis A numeric expansion factor controlling the size of the axis annotations. The default is 2.
#' @param cex.main A numeric expansion factor controlling the size of the main title. The default is 2.5.
#' @param cex.points A numeric expansion factor controlling the size of the points representing the local estimates of the time-varying parameters. The default is 2.
#' @param mar A numeric vector specifying the plot margins in the order bottom, left, top and right. The default is c(5,6,4,2).
#' @param mfrow A numeric vector of length two specifying the graphical layout. The default, c(1,3), arranges the local estimates of \eqn{c(u)}, \eqn{\alpha(u)} and \eqn{\beta(u)} in a single row.
#' @return Data frame who contains c, alpha, beta of GARCH(1,1) model and midpoint each block.
#' @importFrom fGarch garchFit
#' @importFrom graphics abline
#' @importFrom graphics par
#' @importFrom graphics lines
#' @importFrom timeSeries description
#' @export
#' @examples
#' ipsa<-diff(log(indipsa))*100
#' S = 100
#' N = 800
#' tv <- tvParameter(ipsa,S,N)

tvParameter<-function(data,
                      S,
                      N,
                      plot = TRUE,
                      title = NULL,
                      description = NULL,
                      cex.axis = 2,
                      cex.main = 2.5,
                      cex.points = 2,
                      mar = c(5,6,4,2),
                      mfrow = c(1,3)){

  CALL = as.call(match.call())
  omega <- vector("numeric")
  alpha <- vector("numeric")
  beta <- vector("numeric")

  y = data
  T. = length(y)
  S  = S
  N  = N
  # T. = 0:T./T.
  M  = trunc((T.-N)/S+1)
  estacionario <- garchFit(~garch(1,1), data=y, include.mean=F, trace = FALSE)
  G <- estacionario@fit$coef

  for(j in 1:M){
    model= garchFit(~garch(1,1),data=y[(S*(j-1)+1):(S*(j - 1)+N)], include.mean=F, trace = FALSE)
    omega[j]=model@fit$matcoef[1,1]
    alpha[j]=model@fit$matcoef[2,1]
    beta[j]=model@fit$matcoef[3,1]
  }

  t = S * (1:M - 1) + N/2
  t_plot = (t-min(t))/(max(t)-min(t))
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  if (plot == TRUE) {
    par(mfrow = mfrow, mar = mar)
    plot(t_plot, omega, xlim=c(0,max(t_plot)), ylim =c(min(omega)* 0.5, max(omega)*1.5), bty = "n",
         main = expression(~c(u)), las = 1, lwd = 3,
         ylab = "", cex.axis = cex.axis, xlab = "", cex.main = cex.main, cex=cex.points, pch = 20, cex.lab=3)
    abline(h=G[1], lwd=2, lty=2)
    mtext(side=1, text = "u", line=3, cex=1.5)

    plot(t_plot, alpha, xlim=c(0,max(t_plot)), ylim = c(min(alpha)* 0.5, max(alpha)*1.5), bty = "n",
         main = expression(~alpha(u)), las = 1, lwd = 3,
         ylab = "", cex.axis = cex.axis, xlab = "", cex.main = cex.main, cex=cex.points, pch = 20, cex.lab=3)
    abline(h=G[2], lwd=2, lty=2)
    mtext(side=1, text = "u", line=3, cex=1.5)

    plot(t_plot, beta, xlim=c(0,max(t_plot)), ylim = c(min(beta)*0.5, max(beta)*1.5), bty = "n",
         main =expression(~beta(u)), las = 1, lwd = 3, ylab = "", cex.axis = cex.axis, xlab = "",
         cex.main = cex.main, cex=cex.points, pch = 20, cex.lab=3)#, mgp = c(6, 1, 0))
    abline(h=G[3], lwd=2, lty=2)
    mtext(side=1, text = "u", line=3, cex=1.5)
  }

  if(is.null(title)) title = "TV-GARCH Modelling applied Kalman Filter"
  if(is.null(description)) description = description()
  tvparameter <- data.frame(t=t, c=omega, alpha=alpha, beta=beta)

  new("tvGarchKF", title = title, call = CALL, coef = tvparameter, data = data, description = description,
      stationary_coef = as.numeric(G),var.coef = "NOT APPLIED", residuals = "NOT APPLIED", sigma2 = "NOT APPLIED", loglik = "NOT APPLIED",
      aic = "NOT APPLIED", method  = "NOT APPLIED", pars = "NOT APPLIED", type = "NOT APPLIED",
      exponentes = "NOT APPLIED", trig = "NOT APPLIED", arg = "NOT APPLIED")
}
