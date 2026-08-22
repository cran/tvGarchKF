#' Clase tvGarchKFNonParam
#'
#' Define S4 class of non-parametrics methods applied Filter Kalman.
#'
#' @slot title Title.
#' @slot call User call.
#' @slot fitted.values This value it's the conditional variance.
#' @slot tvparameter The original time-varying parameter structure used as input.
#' @slot data A univariate time series, where could be an object of class "ts" or "zoo".
#' @slot Omega Estimation of \eqn{c} parameter.
#' @slot Alpha Estimation of \eqn{\alpha} parameter.
#' @slot Beta Estimation of \eqn{\beta} parameter.
#' @slot Param_model This variable it's the stationary model.
#' @slot method The selected non-parametric estimation method.
#' @slot description Description.
#' @importFrom utils head
#' @importFrom methods slot
#' @importFrom methods slotNames
#' @exportClass tvGarchKFNonParam_Fit
setClass(
  "tvGarchKFNonParam_Fit",
  slots = list(
    title = "character",
    call = "ANY",
    tvparameter = "ANY",
    fitted.values = "numeric",
    Omega = "ANY",
    Alpha = "ANY",
    Beta = "ANY",
    Param_model = "ANY",
    data = "numeric",
    method = "character",
    description = "character"
  )
)

setMethod("show", "tvGarchKFNonParam_Fit", function(object) {
  for (s in slotNames(object)) {
  val <- slot(object, s)
    cat(s,":\n")
    if (is.numeric(val) && !is.null(names(val))) {
      print(val)
    } else if (is.vector(val) && length(val) < 10) {
      cat(" ", paste(val, collapse = " "), "\n")
      cat("\n")
    } else if (is.vector(val)) {
      cat(" ", paste(head(val, 6), collapse = " "), "... (", length(val), "elementos )\n")
      cat("\n")
    } else if (is.matrix(val) || is.data.frame(val)) {
      cat(" [", nrow(val), "x", ncol(val), "]\n", sep = "")
      print(utils::head(val, 15))
      cat("\n")
    } else {
      print(val)
      cat("\n")
    }
}})
