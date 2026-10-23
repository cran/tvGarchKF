#' Clase tvGarchKF
#'
#' Define principal S4 class of parametrics methods.
#'
#' @slot call User call.
#' @slot coef Estimation of coefficients \eqn{c}, \eqn{\alpha} and \eqn{beta}.
#' @slot var.coef Matrix of covariance.
#' @slot residuals Residuals
#' @slot sigma2 Conditional variance.
#' @slot data A univariate time series, where could be an object of class "ts" or "zoo".
#' @slot loglik Log-likelihood.
#' @slot aic AIC value.
#' @slot title Title.
#' @slot description Description.
#' @slot method A character string specifying the parametric estimation approach used to the tv-parameter curves.
#' @slot type A charactertype=c("polynomial","NoLineal","trigonometric"), which indicates what type of function will be assigned to each tv-parameters.
#' @slot pars A vector with initial coefficients (\eqn{c}, \eqn{\alpha} and \eqn{beta}).
#' @slot exponentes A vector for the exponents for the "NoLinear" case.
#' @slot trig A character parameter indicating the trigonometric function to be used, "cos" or "sin".
#' @slot arg This is the value of the argument to be used in the trigonometric function.
#' @slot stationary_coef A vector with stationary coefficients (\eqn{c}, \eqn{\alpha} and \eqn{beta}).
#' @importFrom utils head
#' @importFrom methods slot
#' @importFrom methods slotNames
#' @exportClass tvGarchKF
setClass(
  "tvGarchKF",
  slots = list(
    title = "character",
    call = "ANY",
    coef = "ANY",
    var.coef = "ANY",
    residuals = "ANY",
    sigma2 = "ANY",
    data = "ANY",
    loglik = "ANY",
    aic = "ANY",
    method = "character",
    description = "character",
    pars = "ANY",
    type = "ANY",
    exponentes = "ANY",
    trig = "ANY",
    arg = "ANY",
    stationary_coef = "ANY"
  )
)

setMethod("show", "tvGarchKF", function(object) {
  for (s in slotNames(object)) {
    val <- slot(object, s)
    if(!identical(val,"NOT APPLIED"))
    {
      cat(s,":\n")
      if (is.numeric(val) && !is.null(names(val))) {
      print(val)
      } else if (is.vector(val) && length(val) < 10) {
        cat(" ", paste(val, collapse = " "), "\n")
        cat("\n")
      } else if (is.vector(val)) {
        cat(" ", paste(head(val, 15), collapse = " "), "... (", length(val), "elementos )\n")
        cat("\n")
      } else if (is.matrix(val) || is.data.frame(val)) {
        cat(" [", nrow(val), "x", ncol(val), "]\n", sep = "")
        print(utils::head(val, 15))
        cat("\n")
      } else {
        print(val)
        cat("\n")
      }
    }
  }
})



