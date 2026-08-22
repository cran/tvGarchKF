#' @title Estimate a tv-GARCH(1,1) by using the Kalman Filter.
#' @description This function estimates a time-varying GARCH(1,1) (tv-GARCH(1,1)) model using the Kalman Filter (KF). The model parameters are allowed
#' to evolve smoothly over rescaled time \eqn{(u=t/T)} according to user-specified functional forms. The time-varying parameters \eqn{c(t)},
#' \eqn{\alpha(t)} and \eqn{\beta(t)} determine the dynamics of the condition variance through
#' \eqn{c(u_t) +\alpha(u_t) r^2_{t-1} +\beta(u_t)\sigma_{t-1}, u_t = \frac{t}{T}}. The functional forms of the time-varying
#' parameters are specified by the user and the corresponding coefficients are estimated using the KF-based procedure.
#' @param data A univariate time series, which can be an object of class "ts" or "zoo".
#' @param c A numeric vector containing the initial values of the coefficients defining the function \eqn{c(u)}. These values are used to initialize the estimation procedure.
#' @param alpha A numeric vector containing the initial values of the coefficients defining the function \eqn{\alpha(u)}. These values are used to initialize the estimation procedure.
#' @param beta A numeric vector containing the initial values of the coefficients defining the function \eqn{\beta(u)}. These values are used to initialize the estimation procedure.
#' @param type A character vector specifying the functional form assigned to each time-varying parameter. Available options are "polynomial", "NoLineal" and "trigonometric". The three elements correspond to \eqn{(c(u)),\alpha(u),\beta(u)}.
#' @param exponentes A numeric vector specifying the exponents used in the nonlinear specification when type = "NoLineal".
#' @param trig A character value specifying the trigonometric function used when type = "trigonometric". Available options are "cos" and "sin".
#' @param arg A character string specifying the argument of the trigonometric function used when type = "trigonometric".
#' @param predict A non-negative integer specifying the number of steps ahead for volatility forecasting. The default is 0, in which case no out-of-sample forecasts are computed during model fitting.
#' @param trace.log Logical. If TRUE, additional information about the estimation procedure is printed. The default is FALSE.
#' @param title An optional character string specifying a title associated with the fitted object. The default is NULL.
#' @param description An optional character string providing a user-defined description of the fitted model. The default is NULL.
#' @return A list of class "tvGarchKF".
#' @importFrom stats optim
#' @importFrom stats optimHess
#' @importFrom numDeriv hessian
#' @importFrom numDeriv jacobian
#' @importFrom timeSeries description
#' @importFrom methods new
#' @useDynLib tvGarchKF, .registration = TRUE
#' @export
#' @details
#' The time-varying parameters can be represented using polynomial, nonlinear
#' or trigonometric functional forms. For the polynomial specification,
#' \deqn{c(t) = c_0 + c_1u + c_2 u^2 + \ldots + c_p u^p,}
#' \deqn{\alpha(t) = a_0 + a_1u + a_2 u^2 + \ldots + a_p u^p,}
#' \deqn{\beta(t) = b_0 + b_1u + b_2 u^2 + \ldots + b_p u^p,}
#' where \eqn{u=t/T}, with \eqn{t=1, 2, \ldots, T}.
#' For the nonlinear specification,
#' \deqn{c(t) = c_0 + \sum_{j=1}^k c_j u^{q_{c,j}},}
#' \deqn{\alpha(t) = a_0 + \sum_{j=1}^k a_j u^{q_{\alpha,j}},}
#' \deqn{\beta(t) = b_0 + \sum_{j=1}^k b_j u^{q_{\beta,j}},}
#' where the exponents \eqn{q_{c,j}}, \eqn{q_{\alpha,j}} and \eqn{q_{\beta,j}} are specified through the exponentes argument
#' For the trigonometric specification,
#' \deqn{c(t) = c_0 + c_1 g(u),}
#' \deqn{\alpha(t) = a_0 + a_1 g(u),}
#' \deqn{\beta(t) = b_0 + b_1 g(u),}
#' where (\eqn{g(u)}) denotes the trigonometric function specified through trig, with its argument defined by arg.
#' @examples
#' \donttest{
#' ipsa<-diff(log(indipsa))
#' c <- c(0.05,0.05)
#' alpha <- c(0.05,0.05)
#' beta <- c(0.05,0.05)
#' type_fit <- c("trigonometric","trigonometric","trigonometric")
#' fit<-tvGarchKalmanFit(ipsa,c=c,alpha=alpha,beta=beta,type=type_fit,trig="cos",arg="3*(1-log(u))")
#' }

tvGarchKalmanFit <- function(data,
                             c,
                             alpha,
                             beta,
                             type=c("polynomial","NoLineal","trigonometric"),
                             exponentes = NULL,
                             trig = NULL,
                             arg = NULL,
                             predict = 0,
                             trace.log = FALSE,
                             title = NULL,
                             description = NULL){
  pars = c(c, alpha, beta)
  checkInput(type, c, alpha, beta, length(exponentes))
  CALL = as.call(match.call())

  if(type[1]!="polynomial" | type[2]!="polynomial" | type[3]!="polynomial"){
    aux <- optim(par = pars, fn = tvGarchKalmanLoglike, data = data, c = c, alpha = alpha, beta = beta, type = type, exponentes = exponentes, trig = trig, arg = arg, method = c("BFGS"))
  }
  else{
    aux <- optim(par = pars, fn = tvGarchKalmanLoglike, data = data, c = c, alpha = alpha, beta = beta, type = type, exponentes = exponentes, trig = trig, arg = arg, method = c("CG"))
  }
  par <- round(aux$par,5)

  # if(plot == TRUE){
  #   params_parameter <- tvCoef(x=par, data, c, alpha, beta, type, exponentes, trig, arg, S, N, flag = FALSE)
  #   tvParameter(data = data, S=S, N=N, plot=plot, fit = par, c = params_parameter$c, a = params_parameter$a, b = params_parameter$b)
  # }
  x <- aux$par
  aux. <- tvGarchKalmanPrint(x = x, data = data, c = c, alpha = alpha,
                             beta = beta, type = type, exponentes = exponentes,
                             trig = trig, arg = arg, predict = predict,
                             trace.log = F)

  loglik <- aux.$loglik
  fitted.values <- aux.$sigma

  if (is.null(fitted.values) || all(is.na(fitted.values))) {
    stop("El modelo no pudo ser ajustado: la verosimilitud fue penalizada.")
  }

  aic <- -2*loglik + 2*length(x)/length(data)

  h <- hessian(func = tvGarchKalman, x = x, data = data, c = c, alpha = alpha, beta = beta, type = type, exponentes = exponentes, trig = trig, arg = arg, predict = predict)
  H <- solve(h)/length(data)

  table <- data.frame(rbind(x))
  lista <- c(round(x,5))
  names(lista) <- c(
    paste0("c", seq_len(length(c)) - 1),
    paste0("a", seq_len(length(alpha)) - 1),
    paste0("b", seq_len(length(beta)) - 1))
  colnames(table) <- names(lista)

  HH <- optimHess(par = x, fn = tvGarchKalman, data = data, c = c, alpha = alpha, beta = beta, type = type, exponentes = exponentes, trig = trig, arg = arg, predict = predict)

  vcov_mat <- tryCatch(
    solve(HH),
    error = function(e) matrix(NA, nrow = length(x), ncol = length(x))
  )

  rownames(vcov_mat) <- colnames(vcov_mat) <- names(lista)

  residuals <- as.matrix(data/sqrt(fitted.values))
  colnames(residuals) <- "innovations"

  if(is.null(title)) title = "TV-GARCH Modelling applied Kalman Filter"
  if(is.null(description)) description = description() # {timeSeries}

  pars = list(c=c, alpha=alpha, beta=beta)


  new("tvGarchKF", title = title, call = CALL, coef = lista, var.coef = vcov_mat,
      sigma2 = fitted.values, residuals = residuals, data = data, loglik = loglik,
      aic = aic, method = "Parametric", description = description, pars = pars,
      type = type, exponentes=exponentes, trig = trig, arg = arg, stationary_coef = "NOT APPLIED")

}

