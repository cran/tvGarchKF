#' @title EStimate tv-GARCH(1,1) parameter curves using non-parametric smoothing methods.
#' @description Given a set of initial local estimates of the time-varying parameters, this function applies smoothing spline, local polynomial, or B-spline
#' to obtain smooth estimates tv-GARCH(1,1) parameter curves (\eqn{c(u)}, \eqn{\alpha(u)}, \eqn{\beta(u)}). The tuning parameters can be specified directly
#' by the user or selected automatically using \code{auto.tvGarchKF()}.
#' @param data A univariate time series, which can be an object of class "ts", "zoo" or "xts".
#' @param tvparameter An object containing the initial estimates of the time-varying GARCH parameters \eqn{c(u)}, \eqn{\alpha(u)} and \eqn{\beta(u)}, typically obtained from \code{tvParameter()} or extracted from an object of class "tvGarchKF".
#' @param smoothing A list specifying the tuning parameters for the smoothing spline method. The component spar = \eqn{(spar_c,spar_\alpha,spar_\beta)} controls the amount of smoothing for \eqn{c(u)}, \eqn{\alpha(u)} and \eqn{\beta(u)}, respectively, while nknots = \eqn{(nknots_c,nknots_\alpha,nknots_\beta)} specifies the corresponding number of knots. The optional components methods, control.spar and max_nknots provide additional control over the smoothing procedure.
#' @param bspline A list specifying the tuning parameters for the B-spline method. The component degree specifies the degree of the B-spline basis, whereas df specifies the degrees of freedom used to estimate the parameter curves \eqn{c(u)}, \eqn{\alpha(u)} and \eqn{\beta(u)}. Curve-specific values can be supplied for each time-varying parameter.
#' @param localPoly A list specifying the tuning parameters for the local polynomial method. The component span controls the proportion of observations used in each local fit, whereas degree specifies the polynomial degree. Curve-specific values can be supplied for \eqn{c(u)}, \eqn{\alpha(u)} and \eqn{\beta(u)}.
#' @param method A character string specifying the non-parametric smoothing method. Available options are "Smoothing Spline", "Local Polynomial", and "B-Spline". A single method is applied to all three curves \eqn{c(u)}, \eqn{\alpha(u)}, and \eqn{\beta(u)}, although curve-specific tuning parameters can be specified.
#' @param plot Logical. If TRUE, the estimated time-varying parameter curves are displayed. The default is TRUE.
#' @param title An optional character string specifying a title associated with the fitted object. The default is NULL.
#' @param description An optional character string providing a user-defined description of the fitted model. The default is NULL.
#' @return An object of class "tvGarchKFNonParam" containing the non-parametric estimates of the tv-GARCH parameter curves. The returned object includes:
#' \itemize{
#'    \item \code{tvparameter}: The initial estimates of the time-varying GARCH parameters used as input.
#'    \item \code{fitted.values}: A data frame containing the fitted non-parametric curves \eqn{c(u)}, \eqn{\alpha(u)} and \eqn{\beta(u)} at the observed time points.
#'    \item \code{HattvParameters}: The estimated non-parametric curves evaluated over the complete time domain.
#'    \item \code{method}: The selected non-parametric smoothing method.
#'    \item \code{C, Alpha, Beta}: Fitted objects associated with the non-parametric estimation of (\eqn{c(u)}, \eqn{\alpha(u)}, \eqn{\beta(u)}), respectively.
#'    \item \code{Method-specific tuning parameters}: The tuning parameters associated with the selected non-parametric method, such as spar, nknots, span, degree or df.
#'    \item \code{Check_condition}: Diagnostic information for the local stationarity condition \eqn{(\alpha(u) + \beta(u) < 1)}, including the maximum value of \eqn{(\alpha(u)+\beta(u))} over the observed time domain and an indicator of whether the condition is satisfied.
#' }
#' The resulting object can be used with \code{predict()} to obtain the estimated time-varying parameter curves required as input to \code{tvGarchNonParamKF()}.
#' @importFrom timeSeries description
#' @importFrom methods new
#' @importFrom locfit lp
#' @importFrom splines bs
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
#' fit_par<-auto.tvGarchKF(object=tv_curve,data=model[,1],k_grid=3:8,bs_gam="ps",gam_method="REML",
#'                           spar_grid=seq(0.3,1.3,by=0.01),span_grid=seq(0.50,2.50,by=0.05),
#'                           loess_degree=c(1,2,2), bs_degree_grid=c(2,2,2),tol=1e-3)
#' # Smoothing Spline
#' FitNonParam_ss<-tvGarchNonParamFit(data=model[,1],tv_curve,smoothing=list(spar=c(0.5,0.68,0.8),
#'                                     nknots=c(3,3,3)),plot=TRUE,method="Smoothing Spline")
#' # Local Polynomial
#' FitNonParam_lp<-tvGarchNonParamFit(data=model[,1],tvparameter=tv_curve,
#'                                     localPoly=list(span=c(0.70,2.5,1.4),degree=c(1,2,2)),
#'                                     plot=TRUE,method="Local Polynomial")
#' # B-Spline
#' FitNonParam_bs<-tvGarchNonParamFit(data=model[,1],tvparameter=tv_curve,
#'                                      bspline=list(degree=c(2,2,2),df=c(3,4,2)),
#'                                      plot=TRUE,method="B-Spline")
#' }

tvGarchNonParamFit <- function(
    data,
    tvparameter,
    smoothing = list(
      spar = c(0.5,0.5,0.5),
      nknots = c(12L, 12L, 12L),
      methods = NULL,
      control.spar = list(lower = 0.5, upper = 1.2),
      max_knots = NULL),
    bspline = list(
      degree = 1:3,
      df = c(3,3,4)),
    localPoly = list(
      span = seq(0.10, 0.45, by = 0.05),
      degree = 1:2),
    plot = TRUE,
    method = c("Smoothing Spline", "Local Polynomial", "B-Spline"),
    title = NULL,
    description = NULL
    ){
  data <- as.numeric(data)
  tvparameter <- .normalize_tv_names(tvparameter@coef)
  method <- match.arg(method)

  smoothing_def <- list(
    spar=c(0.5,0.5,0.5),
    nknots = c(12L, 12L, 12L),
    methods = NULL,
    control.spar = list(lower = 0.5, upper = 1.2),
    max_knots = NULL)
  bspline_def <- list(
    degree = 2:3,
    df = c(3,3,4))
  local_def <- list(
    span = seq(0.10, 0.45, by = 0.05),
    degree = 1:2)

  smoothing <- .merge_opts(smoothing, smoothing_def)
  bspline   <- .merge_opts(bspline, bspline_def)
  localPoly <- .merge_opts(localPoly, local_def)

  res <- switch(
    method,
    "Smoothing Spline" = run_smoothing_and_plotting(
      data = data,
      tvparameter = tvparameter,
      spar= smoothing$spar,
      nknots = smoothing$nknots,
      methods = smoothing$methods,
      control.spar = smoothing$control.spar,
      max_knots = smoothing$max_knots,
      plot = plot
    ),
    "Local Polynomial" = run_polynomial_and_local_fit(
      data = data,
      tvparameter = tvparameter,
      span = localPoly$span,
      degree = localPoly$degree,
      plot = plot
    ),
    "B-Spline" = fit_b_spline_models(
      data = data,
      tvparameter = tvparameter,
      degree = bspline$degree,
      df = bspline$df,
      plot = plot
    )
  )

  CALL = as.call(match.call())
  if(is.null(title)) title = "TV-GARCH Modelling Non Parametric applied Kalman Filter"
  if(is.null(description)) description = description()

  switch(method,
         "Smoothing Spline" = new("tvGarchKFNonParam", title = title, call = CALL, Check_condition = res$"Check condition",
                                  tvparameter = res$tvparameter, fitted.values = res$fitted.values,
                                  ParamStationary = res$ParamStationary, HattvParameters = res$HattvParameters, data = data,
                                  Omega = res$Param_model$c, Alpha = res$Param_model$alpha,Beta =  res$Param_model$beta,
                                  smoothing_used = res$smoothing_used, local_used = "NOT APPLIED", best_df = "NOT APPLIED",
                                  best_degree = "NOT APPLIED", Param_model = res$Param_model,
                                  method = method, description = description),
         "Local Polynomial" = new("tvGarchKFNonParam", title = title, call = CALL, Check_condition = res$"Check condition",
                                  tvparameter = res$tvparameter, fitted.values = res$fitted.values,
                                  ParamStationary = res$ParamStationary, HattvParameters = res$HattvParameters, data = data,
                                  Omega = res$Param_model$c, Alpha = res$Param_model$alpha,Beta =  res$Param_model$beta,
                                  local_used = res$local_used, smoothing_used = "NOT APPLIED", best_degree = "NOT APPLIED",
                                  Param_model = res$Param_model, best_df = "NOT APPLIED",
                                  method = method, description = description),
         "B-Spline" = new("tvGarchKFNonParam", title = title, call = CALL, Check_condition = res$"Check condition",
                          tvparameter = res$tvparameter, fitted.values = res$fitted.values,
                          ParamStationary = res$ParamStationary, HattvParameters = res$HattvParameters, data = data,
                          Omega = res$Param_model$c, Alpha = res$Param_model$alpha, Beta =  res$Param_model$beta,
                          best_degree = res$best_degree, best_df = res$best_df,
                          local_used = "NOT APPLIED", smoothing_used = "NOT APPLIED", Param_model = res$Param_model,
                          method = method, description = description)
  )

}


