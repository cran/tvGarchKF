#' Clase tvGarchKFNonParam
#'
#' Define principal S4 class of non-parametrics methods.
#'
#' @slot call User call.
#' @slot Check_condition Stationarity condition
#' @slot fitted.values A data frame containing the estimated non-parametric curves \eqn{c(u)}, \eqn{\alpha(u)} and \eqn{\beta(u)}.
#' @slot tvparameter The original time-varying parameter structure used as input.
#' @slot data A univariate time series, where could be an object of class "ts" or \"zoo".
#' @slot ParamStationary Estimation of stationary coefficients.
#' @slot Omega Estimation of \eqn{c} parameter.
#' @slot Alpha Estimation of \eqn{\alpha} parameter.
#' @slot Beta Estimation of \eqn{\beta} parameter.
#' @slot title Title.
#' @slot HattvParameters Predicted values of the non-parametric parameter curves over the complete time domain.
#' @slot Param_model This variable it's the stationary model.
#' @slot best_degree The degree of the spline basis functions.
#' @slot best_df The degrees of freedom used to estimate the tv-parameter curves \eqn{c(u)}, \eqn{\alpha(u)} and \eqn{\beta(u)}.
#' @slot local_used A list containing local polynomial tuning parameters.
#' @slot smoothing_used A list containing the smoothing spline tuning parameters.
#' @slot method The selected non-parametric estimation method.
#' @slot description Description.
#' @importFrom utils head
#' @importFrom methods slot
#' @importFrom methods slotNames
#' @exportClass tvGarchKFNonParam

setClass(
  "tvGarchKFNonParam",
  slots = list(
    title = "character",
    call = "ANY",
    Check_condition = "ANY",
    tvparameter = "ANY",
    fitted.values = "data.frame",
    ParamStationary = "matrix",
    HattvParameters = "data.frame",
    data = "numeric",
    Omega = "ANY",
    Alpha = "ANY",
    Beta = "ANY",
    best_degree = "ANY",
    best_df = "ANY",
    smoothing_used = "ANY",
    local_used = "ANY",
    Param_model = "ANY",
    method = "character",
    description = "character"
  )
)

setMethod("show", "tvGarchKFNonParam", function(object) {
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
