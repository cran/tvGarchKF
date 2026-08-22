#' Number of Observations
#'
#' Generic function to extract the number of observations from an object.
#' @param object Object class \code{tvGarchKF}.
#' @param ... Others.
#' @return Number of observations.
#' @importFrom stats nobs
#' @exportMethod nobs

setMethod(f = "nobs", signature(object = "tvGarchKF"),
          function(object, ...) {
            length(object@sigma2)
          }
)
