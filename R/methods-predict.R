#' @title Forecast conditional volatility from a tv-GARCH model.
#' @rdname predict-methods
#' @aliases predict,tvGarchKF-method predict,tvGarchKFNonParam-method predict,tvGarchKFNonParam_Fit-method
#' @description
#' This method computes out-of-sample forecasts of the conditional volatility from a fitted tv-GARCH model using the Kalman Filter prediction recursions.
#' Approximate confidence intervals for the forecasted volatility are constructed from the state prediction error variance.
#' @param object An object class \code{tvGarchKF} or \code{tvGarchKFNonParam} or \code{tvGarchKFNonParam_Fit}.
#' @param n.ahead A positive integer specifying the number of steps ahead for which the condition volatility is forecast. The default is 10.
#' @param conf Confidence level used to construct prediction intervals for the forecasted volatility. The default is 0.95.
#' @param crit_val Optional critical value used to construct the confidence intervals. If NULL, the critical value is automatically obtained from the standard normal distribution according to the confidence level specified by \code{conf}.
#' @param interval Logical. If TRUE, confidence intervals for the forecasted conditional volatility are computed and returned. The default is TRUE.
#' @param Plot Logical. If TRUE, the forecasted conditional volatility and the corresponding confidence intervals are displayed. The default is TRUE.
#' @param ... Additional graphical arguments passed to the plot method for objects of class "tvGarchKF_predict".
#' @return An object of class "tvGarchKF_predict" containing the forecast of the conditional volatility and the quantities required to asses their uncertainty.
#' The returned object includes:
#' \itemize{
#'    \item \code{n.ahead}: Forecast horizon.
#'    \item \code{nsample}: Number of observations in the fitted sample.
#'    \item \code{sigma_forecast}: Forecasted conditional volatility (\eqn{\widehat{\sigma}_{t+h}}, for \eqn{h = 1,\ldots, \texttt{n.ahead}}).
#'    \item \code{se_sigma}: Standard errors of the forecasted conditional volatility.
#'    \item \code{X_forecast}: Predicted state values.
#'    \item \code{P_forecast}: State prediction error variances.
#'    \item \code{method}: Estimation method used in the fitted tv-GARCH model.
#'    \item \code{conf}: Confidence level used to construct intervals.
#'    \item \code{crit_val}: Critical values used to construct the confidence intervals.
#'    \item \code{interval}: A data frame containing the forecast horizon, forecasted volatility, standard error, and lower and upper confidence limits. This components is included only when \code{interval = TRUE}.
#' }
#' If \code{Plot = TRUE}, the forecasted volatility and its confidence intervals are also displayed.
#' @importFrom utils tail
#' @importFrom stats qnorm
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
#' predict(fit,n.ahead=30,Plot=TRUE,legend.cex=1.5,legend.pos="topleft")
#' # Smoothing Spline
#' FitNonParam_ss<-tvGarchNonParamFit(data=model[,1],tv_curve,smoothing=list(spar=c(0.5,0.68,0.8),
#'                                     nknots=c(3,3,3)),plot=TRUE,method="Smoothing Spline")
#' TVparameter_ss<-predict(FitNonParam_ss)
#' FitNonParamKF_ss<-tvGarchNonParamKF(model[, 1], TVparameter_ss)
#' predict(FitNonParamKF_ss, n.ahead = 20, Plot=TRUE)
#' # Local Polynomial
#' FitNonParam_lp<-tvGarchNonParamFit(data=model[,1],tvparameter=tv_curve,
#'                                     localPoly=list(span=c(0.70,2.5,1.4),degree=c(1,2,2)),
#'                                     plot=TRUE,method="Local Polynomial")
#' TVparameter_lp<-predict(FitNonParam_lp)
#' FitNonParamKF_lp<-tvGarchNonParamKF(model[, 1], TVparameter_lp)
#' predict(FitNonParamKF_lp, n.ahead = 20, Plot=TRUE)
#' # B-Spline
#' FitNonParam_bs<-tvGarchNonParamFit(data=model[,1],tvparameter=tv_curve,
#'                                      bspline=list(degree=c(2,2,2),df=c(3,4,2)),
#'                                      plot=TRUE,method="B-Spline")
#' TVparameter_bs<-predict(FitNonParam_bs)
#' FitNonParamKF_bs<-tvGarchNonParamKF(model[,1],TVparameter_bs)
#' predict(FitNonParamKF_bs,n.ahead=20,Plot=TRUE)
#' }

setMethod(f = "predict", signature = (object = "tvGarchKF"),
          definition = function(object, n.ahead = 10,
                                conf = 0.95,
                                crit_val = NULL,
                                interval = TRUE,
                                Plot = TRUE,
                                ...){
            if (!inherits(object, "tvGarchKF")) {
              stop("Object must inherit from class 'tvGarchKF'.")
            }
            nsample <- length(object@data)
            rec <- tvGarchKalmanPredict(
              object = object,
              n.ahead = n.ahead,
              trace.log = FALSE
            )
            sigma_h <- tail(rec$sigma, n.ahead)
            X_h     <- tail(rec$X, n.ahead)
            P_h     <- tail(rec$P, n.ahead)

            if (is.null(crit_val)) {
              alpha <- 1 - conf
              crit_val <- qnorm(1 - alpha / 2)
            }

            z <- abs(crit_val[1])
            se_sigma <- sqrt(pmax(P_h, 0)) / (2 * sigma_h)

            lower <- sigma_h * exp(-z * se_sigma / pmax(sigma_h, .Machine$double.eps))
            upper <- sigma_h + z * se_sigma

            lower <- pmax(lower, 0)

            out <- list(
              n.ahead = n.ahead,
              nsample = nsample,
              sigma_forecast = sigma_h,
              se_sigma = se_sigma,
              X_forecast = X_h,
              P_forecast = P_h,
              method = object@method,
              conf = conf,
              crit_val = c(-z, z)
            )
            if (isTRUE(interval)) {
              out$interval <- data.frame(
                h = seq_len(n.ahead),
                sigma = sigma_h,
                se = se_sigma,
                lower = lower,
                upper = upper
              )
            }
            class(out) <- c("tvGarchKF_predict", "list")
            if (isTRUE(Plot)){
              plot(out, ...)
            }
            return(out)
          }
            )

#' @rdname predict-methods
#' @export

setMethod(f = "predict", signature(object = "tvGarchKFNonParam"),
          definition = function(object, n.ahead = 10, ...)
{
    is_nonpar <- !is.null(object@method) && object@method %in% c("Smoothing Spline", "Local Polynomial", "B-Spline")

    if (is_nonpar) {
      if (identical(object@method, "Local Polynomial")) {
        fitted <- object@HattvParameters
        out <- list(
          fitted      = fitted,
          method      = object@method,
          Param_model = object@Param_model
        )
        return(out)
      }

      if (identical(object@method, "Smoothing Spline")) {
        fitted = object@HattvParameters
        out <- list(
          fitted = fitted,
          method = object@method,
          Param_model = object@Param_model
        )
        return(out)
      }

      if (identical(object@method, "B-Spline")) {
        fitted<-object@HattvParameters
        out<- list(
          fitted = fitted,
          method = object@method,
          Param_model = object@Param_model
        )
        return(out)
      }
    }
}
)


#' @rdname predict-methods
#' @export

setMethod(f = "predict", signature = (object = "tvGarchKFNonParam_Fit"),
          definition = function(object, n.ahead = 10,
                                conf = 0.95,
                                crit_val = NULL,
                                interval = TRUE,
                                Plot = TRUE,
                                ...){
            if (!inherits(object, "tvGarchKFNonParam_Fit")) {
              stop("Object must inherit from class 'tvGarchKFNonParam_Fit'.")
            }
            nsample <- length(object@data)
            rec <- tvGarchKalmanPredict(
              object = object,
              n.ahead = n.ahead,
              trace.log = FALSE
            )
            sigma_h <- tail(rec$sigma, n.ahead)
            X_h     <- tail(rec$X, n.ahead)
            P_h     <- tail(rec$P, n.ahead)

            if (is.null(crit_val)) {
              alpha <- 1 - conf
              crit_val <- qnorm(1 - alpha / 2)
            }

            z <- abs(crit_val[1])
            se_sigma <- sqrt(pmax(P_h, 0)) / (2 * sigma_h)

            lower <- sigma_h * exp(-z * se_sigma / pmax(sigma_h, .Machine$double.eps))
            upper <- sigma_h + z * se_sigma

            lower <- pmax(lower, 0)

            out <- list(
              n.ahead = n.ahead,
              nsample = nsample,
              sigma_forecast = sigma_h,
              se_sigma = se_sigma,
              X_forecast = X_h,
              P_forecast = P_h,
              method = object@method,
              conf = conf,
              crit_val = c(-z, z)
            )
            if (isTRUE(interval)) {
              out$interval <- data.frame(
                h = seq_len(n.ahead),
                sigma = sigma_h,
                se = se_sigma,
                lower = lower,
                upper = upper
              )
            }
            class(out) <- c("tvGarchKF_predict", "list")
            if (isTRUE(Plot)){
              plot(out, ...)
            }
            return(out)
          })


