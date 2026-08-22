#' @title Automatic selection of tuning parameters for non-parametric tv-GARCH estimation
#' @description This function provides an automatic procedure for selecting tuning parameters for non-parametric
#' estimation of the tv-GARCH parameter curves. Generalized additive models (GAMs) are used as reference smoothers to estimate the effective
#' degrees of freedom (EDF) associated with each time-varying parameter curve. The resulting EDF values are then used
#' as a common complexity criterion to calibrate comparable tuning parameters for smoothing spline, local polynomial and B-spline estimators.
#' @param object An object class "tvGarchKF" containing the initial estimates of the time-varying GARCH parameters obtained from the local estimation procedure.
#' @param data A numeric vector containing the observed time series used to fit the tv-GARCH model.
#' @param k_grid An integer vector specifying the candidate basis dimensions (k) evaluated when fitting the GAM reference models. The default is 3:8.
#' @param bs_gam A character string indicating the spline basis used in the GAM reference models. The default "ps" corresponding to penalized B-splines.
#' @param gam_method A character string specifying the criterion used to estimate the smoothing parameters of the GAM reference models. The default is "REML". Other available choices include "ML" and "GCV.Cp".
#' @param spar_grid A numeric vector containing candidate smoothing parameters used to select the smoothing spline whose effective degrees of freedom most closely match those of the corresponding GAM reference models. The default is seq(0.3, 1.3, by = 0.01).
#' @param span_grid A numeric vector containing candidate span values used to select the local polynomial smoother whose effective degrees of freedom most closely match those of the corresponding GAM reference models. The default is seq(0.5, 2.5, by = 0.05).
#' @param loess_degree An integer vector specifying the polynomial degree used for the local polynomial approximation of each time-varying parameter curve. The elements correspond, respectively, to \eqn{c(u)}, \eqn{\alpha(u)} and \eqn{\beta(u)}. The default is c(1,2,2).
#' @param bs_degree_grid An integer vector specifying the degree of the B-spline basis used for each time-varying parameter curve. The elements correspond, respectively, to \eqn{c(u)}, \eqn{\alpha(u)} and \eqn{\beta(u)}. The default is c(2,2,2).
#' @param tol A positive numeric value specifying the tolerance used in the EDF-matching procedure. The default is 0.001.
#' @details
#' The automatic selection procedure first fits a GAM reference model to each initial time-varying parameter curve and extracts its effective degrees of freedom.
#' These EDF values provide a common measure of curve complexity. For each non-parametric estimator, the function then searches over the corresponding tuning-parameter
#' grid and selects the specification whose effective degrees of freedom most closely match those of the GAM reference model.
#'
#' Consequently, the smoothing spline, local polynomial and B-spline estimators are calibrated to have approximately comparable levels of smoothness,
#' despite being bases on different smoothing mechanisms.
#' @importFrom locfit lp
#' @importFrom stats time
#' @import npreg
#' @return An object of class "auto.tvGarchKF" containing the fitted GAM reference models and
#' the tuning parameters automatically selected for smoothing spline, local polynomial, and B-spline estimators.
#'
#' The returned object includes the estimated effective degrees of freedom
#' (edf), smoothing parameters (spar), local polynomial spans (span), polynomial and spline
#' degrees (degree), basis dimensions (k), numbers of knots (nknots), and measures of
#' agreement between the effective degrees of freedom of each non-parametric estimator
#' and those of the corresponding GAM reference model.
#'
#' The selected tuning parameters are organised
#' by estimation method and time-varying parameter curve, c(u), \eqn{\alpha(u)} and \eqn{\beta(u)},
#' and can subsequently be supplied to \code{tvGarchNonParamFit()} for non-parametric estimation of the time-varying GARCH parameter curves.
#' @export
#' @examples
#' library("tvGarchKF")
#' n = 2000
#' seed = 321
#' type = c( "polynomial","trigonometric","polynomial")
#' model<-tvGarch_Sim(n=n,c=c(0.18,0.10),alpha=c(0.45,0.4),beta=c(0,1.05,-1.05,0.1),
#'                       type=type,trig="cos",arg="2*pi*u",seed=seed)
#' tv_curve<-tvParameter(model[,1],S=120,N=600,plot=TRUE,cex.axis=3,cex.main=3,cex.points=3,
#'                         mar=c(5,6,4,2),mfrow=c(1,3))
#' fit_par<-auto.tvGarchKF(object=tv_curve,data=model[,1],k_grid=3:8,bs_gam="ps",gam_method="REML",
#'                           spar_grid=seq(0.3,1.3,by=0.01),span_grid=seq(0.50,2.50,by=0.05),
#'                           loess_degree=c(1,2,2), bs_degree_grid=c(2,2,2),tol=1e-3)

auto.tvGarchKF <- function(object, data, k_grid = 3:8, bs_gam = "ps",
                           gam_method = "REML", spar_grid = seq(0.3, 1.3, by = 0.01),
                           span_grid = seq(0.50, 2.50, by = 0.05), loess_degree = c(1, 2, 2),
                           bs_degree_grid = c(2, 2, 2), tol = 1e-3) {

  if (!inherits(object, "tvGarchKF")) {
    stop("Object must be of class 'tvGarchKF'.")
  }

  tvparameter <- .normalize_tv_names(as.data.frame(object@coef))

  n <- length(data)
  t <- tvparameter$t
  u <- t/n

  # =========================================================
  # 1) GAM como puente: obtener k, sp, edf por curva
  # =========================================================

  best_k_omega <- select_best_k_gam(
    x = u,
    y = tvparameter$c,
    k_grid = k_grid,
    bs = bs_gam,
    method = gam_method
  )

  best_k_alpha <- select_best_k_gam(
    x = u,
    y = tvparameter$alpha,
    k_grid = k_grid,
    bs = bs_gam,
    method = gam_method
  )

  best_k_beta <- select_best_k_gam(
    x = u,
    y = tvparameter$beta,
    k_grid = k_grid,
    bs = bs_gam,
    method = gam_method
  )

  selected_gam <- data.frame(
    parameter = c("c", "alpha", "beta"),
    k = c(best_k_omega$k, best_k_alpha$k, best_k_beta$k),
    sp = c(best_k_omega$sp, best_k_alpha$sp, best_k_beta$sp),
    edf = c(best_k_omega$edf, best_k_alpha$edf, best_k_beta$edf)
  )

  # =========================================================
  # 2) Match GAM -> npreg::ss mediante EDF
  # =========================================================

  selected_smoothing <- match_ss_spar_by_edf_tv(
    tvparameter = tvparameter,
    best_k_omega = best_k_omega,
    best_k_alpha = best_k_alpha,
    best_k_beta  = best_k_beta,
    data_length = n,
    spar_grid = spar_grid
  )

  selected_smoothing <- data.frame(
    parameter = selected_smoothing$parameter,
    spar = selected_smoothing$ss_spar,
    nknots = selected_smoothing$ss_nknots,
    edf_target = selected_smoothing$gam_edf,
    edf_ss = selected_smoothing$ss_edf,
    edf_distance = selected_smoothing$edf_distance
  )

  # =========================================================
  # 3) Match GAM -> loess mediante EDF
  # =========================================================

  selected_local <- match_loess_span_by_edf_tv(
    tvparameter = tvparameter,
    best_k_omega = best_k_omega,
    best_k_alpha = best_k_alpha,
    best_k_beta  = best_k_beta,
    data_length = n,
    degree = loess_degree,
    span_grid = span_grid
  )

  selected_local <- data.frame(
    parameter = selected_local$parameter,
    span = selected_local$loess_span,
    degree = selected_local$loess_degree,
    edf_target = selected_local$gam_edf,
    edf_loess = selected_local$loess_edf,
    edf_distance = selected_local$edf_distance
  )

  # =========================================================
  # 4) B-spline: usar k como df objetivo práctico
  # =========================================================

  selected_bspline <- data.frame(
    parameter = c("c", "alpha", "beta"),
    degree = bs_degree_grid,
    df = c(best_k_omega$k, best_k_alpha$k, best_k_beta$k),
    edf_target = c(best_k_omega$edf, best_k_alpha$edf, best_k_beta$edf)
  )

  # =========================================================
  # 5) Salida organizada
  # =========================================================

  out <- list(
    class_input = class(object),
    n = n,
    tol = tol,
    grid_time = t,

    gam_reference = selected_gam,

    selected_parameters = list(
      smoothing_spline = selected_smoothing,
      local_polynomial = selected_local,
      b_spline = selected_bspline
    ),

    details = list(
      smoothing_spline = list(
        spar = selected_smoothing$spar,
        nknots = selected_smoothing$nknots,
        table = selected_smoothing
      ),

      local_polynomial = list(
        span = selected_local$span,
        degree = selected_local$degree,
        table = selected_local
      ),

      b_spline = list(
        degree = selected_bspline$degree,
        df = selected_bspline$df,
        table = selected_bspline
      )
    )
  )

  class(out) <- c("auto.tvGarchKF", "list")
  out
}

