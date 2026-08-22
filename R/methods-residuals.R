#' @title Residuals method.
#' @param object Object class \code{tvGarchKF}
#' @param as.zoo Boolean flag.
#' @param ... Others.
#' @importFrom zoo zoo
#' @exportMethod residuals

setMethod(f = "residuals", signature(object = "tvGarchKF"), definition = function(object, as.zoo = TRUE, ...)
{
  if (as.zoo == TRUE) {
    object@residuals <- zoo(object@residuals, order.by = index(object@data))
  }
  return(object@residuals)
}
)
