#' @title Plot method for \code{tvGarchKF}, \code{tvGarchKFNonParam} and \code{tvGarchKFNonParam_Fit} objects.
#' @rdname plot-methods
#' @aliases plot,tvGarchKF,missing-method  plot,tvGarchKFNonParam_Fit,missing-method  plot,tvGarchKFNonParam,missing-method
#' @description
#' This method provides graphical summaries for fitted tv-GARCH models. By default, it displays the estimated conditional volatility. Alternatively,
#' setting \code{which = "tvparameters"} produces a three-panel plot of the estimated time-varying GARCH parameters \eqn{c(u)}, \eqn{\alpha(u)} and \eqn{\beta(u)}.
#' The function also allows the user to superimpose fitted parametric curves and stationary GARCH benchmarks, and to customize the main graphical elements.
#' @param x An object class \code{tvGarchKF} or \code{tvGarchKFNonParam_Fit} or \code{tvGarchKFNonParam}.
#' @param y No use.
#' @param which A character string indicating the type of plot. By default "volatility", the function displays the estimated conditional volatility \eqn{\widehat{\sigma}_t}. Setting which = "tvparameters" produces a three-panel plot of the estimated time-varying GARCH parameters \eqn{c(u)}, \eqn{\alpha(u)}, and \eqn{\beta(u)}.
#' @param x.index An optional vector defining the horizontal axis. If NULL, the observations are plotted against their sequential index. This argument is useful for displaying the series against calendar dates or another user-defined time scale.
#' @param tvParameter An optional object containing the block estimates of the time-varying GARCH parameters. This argument is used only when which = "tvparameters", in which case the fitted parametric curves are superimposed on the corresponding block estimates.
#' @param type A character vector specifying the functional form assigned to each time-varying parameter. Available specifications include "polynomial", "NoLineal" and "trigonometric".
#' @param trig A character vector indicating the trigonometric function, "cos" or "sin", used when type = "trigonometric".
#' @param arg A character vector specifying the argument of the trigonometric function. This argument is only used when type = "trigonometric".
#' @param col.points Color used for the points representing the estimated time-varying parameters. Default is "black".
#' @param col.fit Color used for the fitted time-varying parameter curves. If NULL, colors are assigned internally.
#' @param col.stationary Color used for the stationary GARCH reference line. Default is "black".
#' @param lwd.y Line width used for the estimated conditional volatility. Default is 3.
#' @param lwd.fit Line width used for the fitted time-varying parameter curves. Default is 3.
#' @param lwd.stationary Line width used for the stationary GARCH reference line. Default is 2.
#' @param lty.fit Line type used for the fitted time-varying parameter curves. Default is 1.
#' @param lty.stationary Line type used for the stationary GARCH reference line. Default is 2.
#' @param pch Plotting symbol used for the points representing the estimated time-varying parameters. Default is 20.
#' @param cex.points Expansion factor controlling the size of the plotted points. Default is 3.
#' @param cex.axis Expansion factor controlling the size of the axis annotations. Default is 2.
#' @param cex.lab Expansion factor controlling the size of the axis labels Default is 2.
#' @param cex.main Expansion factor controlling the size of the main title Default is 2.
#' @param cex.legend Expansion factor controlling the size of the legend text. Default is 1.2.
#' @param mfrow Numeric vector of length two specifying the graphical layout. The default, \code{c(1,3)}, arranges the plots of \eqn{c(u)}, \eqn{\alpha(u)} and \eqn{\beta(u)} in a single row.
#' @param mar Numeric vector specifying the plot margins in the order bottom, left, top and right. Default is \code{c(5,5,4,2)}.
#' @param xlab Label for the horizontal axis. Default is "Time".
#' @param ylab Label for the vertical axis. Default is "Log-Returns".
#' @param main Optional main title for the plot. If NULL, no user-specified title is added.
#' @param legend Logical. If TRUE, a legend is added to the plot. Default is TRUE.
#' @param legend.pos Position of the legend. Default is "bottomleft".
#' @param legend.label Optional character vector specifying the legend labels. If NULL, default labels are generated internally.
#' @param bty Character specifying the box type. Default is "n", so no box is drawn around the plot.
#' @param las Numeric value controlling the orientation of axis labels. Default is 1, corresponding to horizontal labels.
#' @param ann Logical. If False, the default plot annotation is suppressed.
#' @param ... Additional graphical arguments passed to the underlying base \code{R} plotting functions.
#' @importFrom zoo index
#' @importFrom graphics title
#' @importFrom graphics mtext
#' @importFrom stats is.ts
#' @export
#' @examples
#' \donttest{
#' library("tvGarchKF")
#' n = 2000
#' seed = 321
#' type = c("polynomial","trigonometric","polynomial")
#' model<-tvGarch_Sim(n=n,c=c(0.18,0.10),alpha=c(0.45,0.4),beta=c(0,1.05,-1.05,0.1),
#'                       type=type,trig="cos",arg="2*pi*u",seed=seed)
#' tv_curve<-tvParameter(model[,1],S=120,N=600,plot=TRUE,cex.axis=3,cex.main=3,cex.points=3,
#'                         mar=c(5,6,4,2),mfrow=c(1,3))
#' # Parametric
#' c0     <- c(0.20,0.08)
#' alpha0 <- c(0.40,0.45)
#' beta0  <- c(0.00,0.95,-0.95,0.1)
#' fit<-tvGarchKalmanFit(data=model[,1],c=c0,alpha=alpha0,beta=beta0,type=type,trig="cos",arg="2*pi*u")
#' plot(fit, which = "volatility", mar = c(5, 4, 4, 2))
#' plot(fit,which = "tvparameters",tvParameter = tv_curve,
#'         type = c("linear", "trigonometric", "cubic"),trig = "cos",arg = "2*pi*u")
#'
#' # Smoothing Spline
#' FitNonParam_ss<-tvGarchNonParamFit(data=model[,1],tv_curve,smoothing=list(spar=c(0.5,0.68,0.8),
#'                                     nknots=c(3,3,3)),plot=TRUE,method="Smoothing Spline")
#' plot(FitNonParam_ss,which = "tvparameters",tvParameter = tv_curve)
#'
#' TVparameter_ss<-predict(FitNonParam_ss)
#' FitNonParamKF_ss<-tvGarchNonParamKF(model[, 1], TVparameter_ss)
#' plot(FitNonParamKF_ss)
#'
#' # Local Polynomial
#' FitNonParam_lp<-tvGarchNonParamFit(data=model[,1],tvparameter=tv_curve,
#'                                     localPoly=list(span=c(0.70,2.5,1.4),degree=c(1,2,2)),
#'                                     plot=TRUE,method="Local Polynomial")
#' plot(FitNonParam_lp,which = "tvparameters",tvParameter = tv_curve)
#'
#' TVparameter_lp<-predict(FitNonParam_lp)
#' FitNonParamKF_lp<-tvGarchNonParamKF(model[, 1], TVparameter_lp)
#' plot(FitNonParamKF_lp)
#'
#' # B-Spline
#' FitNonParam_bs<-tvGarchNonParamFit(data=model[,1],tvparameter=tv_curve,
#'                                      bspline=list(degree=c(2,2,2),df=c(3,4,2)),
#'                                      plot=TRUE,method="B-Spline")
#' plot(FitNonParam_bs,which = "tvparameters",tvParameter = tv_curve)
#'
#' TVparameter_bs<-predict(FitNonParam_bs)
#' FitNonParamKF_bs<-tvGarchNonParamKF(model[,1],TVparameter_bs)
#' plot(FitNonParamKF_bs)
#' }

setMethod(f = "plot", signature(x = "tvGarchKF"), definition =
            function(x, y,
                     which = c("volatility", "tvparameters"),
                     x.index = NULL,
                     tvParameter = NULL,
                     type = NULL,
                     trig = NULL,
                     arg = NULL,

                     ## colores
                     col.points = "black",
                     col.fit = NULL,
                     col.stationary = "black",

                     ## líneas
                     lwd.y = 3,
                     lwd.fit = 3,
                     lwd.stationary = 2,
                     lty.fit = 1,
                     lty.stationary = 2,

                     ## símbolos y texto
                     pch = 20,
                     cex.points = 3,
                     cex.axis = 2,
                     cex.lab = 2,
                     cex.main = 2,
                     cex.legend = 1.2,

                     ## disposición gráfica
                     mfrow = c(1, 3),
                     mar = c(5, 5, 4, 2),

                     ## etiquetas
                     xlab = "Time",
                     ylab = "Log-Returns",
                     main = NULL,

                     ## leyenda
                     legend = TRUE,
                     legend.pos = "bottomleft",
                     legend.label = NULL,

                     ## opciones base R
                     bty = "n",
                     las = 1,
                     ann = FALSE,
                     ...)
            {
              which <- match.arg(which)
              if (which == "tvparameters") {

                if (is.null(tvParameter)) {
                  stop("For which = 'tvparameters', the argument 'tvParameter' must be provided.")
                }
                if (is.null(col.fit)) {
                  col.fit <- switch(
                    x@method,
                    "Smoothing Spline"  = "deeppink",
                    "Local Polynomial" = "orange",
                    "B-Spline"         = "blue",
                    "Parametric" = "black"
                  )
                }
                if (is.null(type)) type <- x@type
                if (is.null(trig)) trig <- x@trig
                if (is.null(arg))  arg  <- x@arg

                return(
                  plot_tvGarchKalmanFit(
                    fit = x,
                    tv_curve = tvParameter,
                    type = type,
                    trig = trig,
                    arg = arg,
                    lwd.fit = lwd.fit,
                    pch = pch,
                    cex.points = cex.lab,
                    cex.axis = cex.axis,
                    cex.main = cex.main,
                    mar = mar,
                    ...
                  )
                )

              }
              if (which == "volatility") {
                if (is.null(col.fit)) {
                  col.fit <- "red"
                }
                y <- x@data
                sigma <- x@sigma2
                if (is.null(x.index)) {

                  if (is.ts(y) || inherits(y, "zoo") || inherits(y, "xts")) {
                    x.index <- index(y)
                  } else {
                    x.index <- seq_along(y)
                  }

                } else {

                  if (length(x.index) != length(y)) {
                    stop("'x.index' must have the same length as the observed series.")
                  }
                }
                oldpar <- par(no.readonly = TRUE)
                on.exit(par(oldpar))
                par(mar = mar)

                plot(x.index, y,
                     type = "l",
                     bty = bty,
                     main = "",
                     las = las,
                     lwd = lwd.y,
                     pch = pch,
                     cex.axis = cex.axis,
                     cex.lab = cex.lab,
                     ann = ann,
                     ...)

                mtext(side = 1, text = xlab, line = 3, cex = cex.lab)
                mtext(side = 2, text = ylab, line = 3, cex = cex.lab)
                lines(x.index, -sqrt(sigma),
                      col = col.fit,
                      lwd = lwd.fit,
                      lty = lty.fit)
                legend.label <- "tv GARCH parametric"
                title(main = main, cex.main = cex.main)
                if (legend) {
                  legend(legend.pos,
                         legend = legend.label,
                         col = col.fit,
                         lty = lty.fit,
                         lwd = lwd.fit,
                         cex = cex.legend,
                         bty = "n")
                }
              }
            }
)

#' @rdname plot-methods
#' @export
setMethod(f = "plot", signature(x = "tvGarchKFNonParam_Fit"), definition =
            function(
                  x,y,
                  which = c("volatility", "tvparameters"),
                  x.index = NULL,
                  tvParameter = NULL,
                  type = NULL,
                  trig = NULL,
                  arg = NULL,

                  ## colores
                  col.points = "black",
                  col.fit = NULL,
                  col.stationary = "black",

                  ## líneas
                  lwd.y = 3,
                  lwd.fit = 3,
                  lwd.stationary = 2,
                  lty.fit = 1,
                  lty.stationary = 2,

                  ## símbolos y texto
                  pch = 20,
                  cex.points = 3,
                  cex.axis = 2,
                  cex.lab = 2,
                  cex.main = 2,
                  cex.legend = 1.2,

                  ## disposición gráfica
                  mfrow = c(1, 3),
                  mar = c(5, 5, 4, 2),

                  ## etiquetas
                  xlab = "Time",
                  ylab = "Log-Returns",
                  main = NULL,

                  ## leyenda
                  legend = TRUE,
                  legend.pos = "bottomleft",
                  legend.label = NULL,

                  ## opciones base R
                  bty = "n",
                  las = 1,
                  ann = FALSE,
                  ...)
          {
              which <- match.arg(which)
              if (which == "tvparameters") {

                if (is.null(tvParameter)) {
                  stop("For which = 'tvparameters', the argument 'tvParameter' must be provided.")
                }
                if (is.null(col.fit)) {
                  col.fit <- switch(
                    x@method,
                    "Smoothing Spline"  = "deeppink",
                    "Local Polynomial" = "orange",
                    "B-Spline"         = "blue",
                    "Parametric" = "black"
                  )}
                return(
                  plot_tvGarchNonParamFit(
                    object = x,
                    tv_curve = tvParameter,
                    method = x@method,
                    col.points = col.points,
                    col.fit = col.fit,
                    col.stationary = col.stationary,
                    lwd.fit = lwd.fit,
                    lwd.stationary = lwd.stationary,
                    lty.stationary = lty.stationary,
                    pch = pch,
                    cex.points = cex.points,
                    cex.axis = cex.axis,
                    cex.main = cex.main,
                    mar = mar,
                    ...
                  ))
              }
              if (which == "volatility") {
                if (is.null(col.fit)) {
                  col.fit <- "red"
                }
                y <- x@data
                sigma <- x@fitted.values
                if (is.null(x.index)) {

                  if (is.ts(y) || inherits(y, "zoo") || inherits(y, "xts")) {
                    x.index <- index(y)
                  } else {
                    x.index <- seq_along(y)
                  }

                } else {

                  if (length(x.index) != length(y)) {
                    stop("'x.index' must have the same length as the observed series.")
                  }
                }
                oldpar <- par(no.readonly = TRUE)
                on.exit(par(oldpar))
                par(mar = mar)

                plot(x.index, y,
                     type = "l",
                     bty = bty,
                     main = "",
                     las = las,
                     lwd = lwd.y,
                     pch = pch,
                     cex.axis = cex.axis,
                     cex.lab = cex.lab,
                     ann = ann,
                     ...)

                mtext(side = 1, text = xlab, line = 3, cex = cex.lab)
                mtext(side = 2, text = ylab, line = 3, cex = cex.lab)
                lines(x.index, -sqrt(sigma),
                      col = col.fit,
                      lwd = lwd.fit,
                      lty = lty.fit)
                legend.label <- x@method
                title(main = main, cex.main = cex.main)
                if (legend) {
                  legend(legend.pos,
                         legend = legend.label,
                         col = col.fit,
                         lty = lty.fit,
                         lwd = lwd.fit,
                         cex = cex.legend,
                         bty = "n")
                }
              }
}
)
#' @rdname plot-methods
#' @export
setMethod(f = "plot", signature(x = "tvGarchKFNonParam"), definition =
            function(x, y,
                     which = c("volatility", "tvparameters"),
                     x.index = NULL,
                     tvParameter = NULL,
                     type = NULL,
                     trig = NULL,
                     arg = NULL,

                     ## colores
                     col.points = "black",
                     col.fit = NULL,
                     col.stationary = "black",

                     ## líneas
                     lwd.y = 3,
                     lwd.fit = 3,
                     lwd.stationary = 2,
                     lty.fit = 1,
                     lty.stationary = 2,

                     ## símbolos y texto
                     pch = 20,
                     cex.points = 3,
                     cex.axis = 2,
                     cex.lab = 2,
                     cex.main = 2,
                     cex.legend = 1.2,

                     ## disposición gráfica
                     mfrow = c(1, 3),
                     mar = c(5, 5, 4, 2),

                     ## etiquetas
                     xlab = "Time",
                     ylab = "Log-Returns",
                     main = NULL,

                     ## leyenda
                     legend = TRUE,
                     legend.pos = "bottomleft",
                     legend.label = NULL,

                     ## opciones base R
                     bty = "n",
                     las = 1,
                     ann = FALSE,
                     ...)
          {
              which <- match.arg(which)
              if (which == "tvparameters") {

                if (is.null(tvParameter)) {
                  stop("For which = 'tvparameters', the argument 'tvParameter' must be provided.")
                }
                if (is.null(col.fit)) {
                  col.fit <- switch(
                    x@method,
                    "Smoothing Spline"  = "deeppink",
                    "Local Polynomial" = "orange",
                    "B-Spline"         = "blue",
                    "Parametric" = "black"
                  )}
                return(
                  plot_tvGarchNonParamFit(
                    object = x,
                    tv_curve = tvParameter,
                    method = x@method,
                    col.points = col.points,
                    col.fit = col.fit,
                    col.stationary = col.stationary,
                    lwd.fit = lwd.fit,
                    lwd.stationary = lwd.stationary,
                    lty.stationary = lty.stationary,
                    pch = pch,
                    cex.points = cex.points,
                    cex.axis = cex.axis,
                    cex.main = cex.main,
                    mar = mar,
                    ...
                  ))
              }
              if (which == "volatility") {
                if (is.null(col.fit)) {
                  col.fit <- "red"
                }
                y <- x@data
                sigma <- x@fitted.values
                if (is.null(x.index)) {

                  if (is.ts(y) || inherits(y, "zoo") || inherits(y, "xts")) {
                    x.index <- index(y)
                  } else {
                    x.index <- seq_along(y)
                  }

                } else {

                  if (length(x.index) != length(y)) {
                    stop("'x.index' must have the same length as the observed series.")
                  }
                }
                oldpar <- par(no.readonly = TRUE)
                on.exit(par(oldpar))
                par(mar = mar)

                plot(x.index, y,
                     type = "l",
                     bty = bty,
                     main = "",
                     las = las,
                     lwd = lwd.y,
                     pch = pch,
                     cex.axis = cex.axis,
                     cex.lab = cex.lab,
                     ann = ann,
                     ...)

                mtext(side = 1, text = xlab, line = 3, cex = cex.lab)
                mtext(side = 2, text = ylab, line = 3, cex = cex.lab)
                lines(x.index, -sqrt(sigma),
                      col = col.fit,
                      lwd = lwd.fit,
                      lty = lty.fit)
                legend.label <- x@method
                title(main = main, cex.main = cex.main)
                if (legend) {
                  legend(legend.pos,
                         legend = legend.label,
                         col = col.fit,
                         lty = lty.fit,
                         lwd = lwd.fit,
                         cex = cex.legend,
                         bty = "n")
                  }}})
