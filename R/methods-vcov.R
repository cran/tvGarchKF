#' @title Vcov method.
#' @param object Object class \code{tvGarchKF}
#' @param ... Others.
#' @exportMethod vcov

setMethod(f = "vcov", signature(object = "tvGarchKF"), definition = function(object, ...)
{
  return(object@var.coef)
}
)
