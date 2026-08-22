#' @title Fitted method.
#' @rdname fitted-methods
#' @aliases fitted,tvGarchKF-method fitted,tvGarchKFNonParam-method fitted,tvGarchKFNonParam_Fit-method
#' @param object Object class \code{tvGarchKF} or \code{tvGarchKFNonParam} or \code{tvGarchKFNonParam_Fit}.
#' @export

setMethod(f = "fitted", signature(object = "tvGarchKF"), definition = function(object)
{
  results <- object@sigma2
  return(results)
})
#' @rdname fitted-methods
#' @export
setMethod(f = "fitted", signature(object = "tvGarchKFNonParam"), definition = function(object)
{
  is_nonpar <- !is.null(object@method) &&
    object@method %in% c("Smoothing Spline", "Local Polynomial", "B-Spline")
  if (is_nonpar) {
    if (identical(object@method, "Local Polynomial")) {
      out<- object@fitted.values
      return(out)
    }
    if (identical(object@method, "Smoothing Spline")) {
      out<- object@fitted.values
      return(out)
    }

    if (identical(object@method, "B-Spline")) {
      out<- object@fitted.values
      return(out)
    }
  }
})
#' @rdname fitted-methods
#' @export
setMethod(f = "fitted", signature(object = "tvGarchKFNonParam_Fit"), definition = function(object)
{
  results <- object@fitted.values
  return(results)
})


