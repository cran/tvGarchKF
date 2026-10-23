#' @title Summary method.
#' @rdname summary-methods
#' @aliases summary,tvGarchKF-method summary,tvGarchKFNonParam-method
#' @param object Object class \code{tvGarchKF} or \code{tvGarchKFNonParam}.
#' @param ... Others.
#' @importFrom stats pnorm
#' @exportMethod summary

setMethod(f = "summary", signature(object = "tvGarchKF"), definition = function(object, ...)
{
  cat("\n")
  if (!is.null(object@description)) cat("Date:", object@description, "\n")
  cat("\n")
  cat("Time-Varying GARCH(1,1) via Kalman Filter\n")
  cat("========================================\n\n")

  if (!is.null(object@type) && !is.null(object@pars)) {

    cat("Time-varying parameter functions:\n\n")

    cat(
      "c(u)     =",
      tv_function_expression(
        type = object@type[1],
        par_name = "c",
        npar = length(object@pars$c),
        trig = object@trig,
        arg = object@arg,
        exponentes = object@exponentes
      ),
      "\n"
    )

    cat(
      "alpha(u) =",
      tv_function_expression(
        type = object@type[2],
        par_name = "a",
        npar = length(object@pars$alpha),
        trig = object@trig,
        arg = object@arg,
        exponentes = object@exponentes
      ),
      "\n"
    )

    cat(
      "beta(u)  =",
      tv_function_expression(
        type = object@type[3],
        par_name = "b",
        npar = length(object@pars$beta),
        trig = object@trig,
        arg = object@arg,
        exponentes = object@exponentes
      ),
      "\n\n"
    )
  }

  cat("========================================\n\n")
  cat("Coefficients:\n\n")

  aux1 <- object@coef
  if (is.null(aux1) || length(aux1) == 0) {
    cat("There is no 'coef' in the object; omitting the coefficients table.\n")
  } else {
    # Nombres solo si hay 6 parámetros
    if (is.null(names(aux1)) && length(aux1) == 6L) {
      names(aux1) <- c("c0", "a0", "a1", "a2", "b0", "b1")
    }

    # Errores estándar
    se <- if (!is.null(object@var.coef) && length(diag(object@var.coef)) == length(aux1)) {
      sqrt(diag(object@var.coef))
    } else {
      rep(NA_real_, length(aux1))
    }

    zval <- aux1 / se
    pval <- 2 * (1 - stats::pnorm(abs(zval)))

    Tab <- data.frame(
      Estimate   = as.numeric(aux1),
      `Std. Error` = as.numeric(se),
      `z-value`  = as.numeric(zval),
      `Pr(>|z|)` = as.numeric(pval),
      row.names  = names(aux1)
    )
    print(Tab)
    cat("\n")
  }

  if (!is.null(object@loglik)) cat("Log-Likelihood:", round(object@loglik, 4), "\n")
  if (!is.null(object@aic))    cat("AIC:", round(object@aic, 4), "\n")
  # if (!is.null(object@bic))    cat("BIC:", round(object@bic, 4), "\n")
  cat("\n")

  invisible(object)

  invisible(object)
}
)

#' @rdname summary-methods
#' @export
setMethod(f = "summary", signature(object = "tvGarchKFNonParam"), definition = function(object, ...)
{

  is_nonpar <- !is.null(object@method) &&
    object@method %in% c("Smoothing Spline", "Local Polynomial", "B-Spline")

  if (is_nonpar) {
    cat("\n")
    if (!is.null(object@description)) cat("Date:", object@description, "\n")
    cat("Time-Varying GARCH(1,1) via Non-parametric Method\n")
    cat("=================================================\n")
    cat("Nonparametric Method:", object@method, "\n\n")

    if (identical(object@method, "Local Polynomial")) {

      max_val <- unname(object@Check_condition["Max."])

      cat("Stationarity condition (alpha(u) + beta(u) < 1)\n")
      cat("=================================================\n")
      cat("Maximum value:", round(max_val, 4), "\n")
      cat("Condition satisfied?", max_val < 1, "\n\n")

      cat("Local Polynomial selected parameters\n")
      cat("====================================\n")

      local_tbl <- data.frame(
        parameter = c("c", "alpha", "beta"),
        span = object@local_used$span,
        degree = object@local_used$degree
      )

      print(local_tbl, row.names = FALSE)

      cat("Performance metrics: RMSE and MAE \n")
      cat("=================================================\n")
      cat("Nonparametric Method:", object@method, "\n\n")

      res <- eval_model(object@tvparameter[,2:4], object@fitted.values, na.rm = TRUE)
      print_eval_model(res)
    }

    if (identical(object@method, "Smoothing Spline")) {

      max_val <- unname(object@Check_condition["Max."])

      cat("Stationarity condition (alpha(u) + beta(u) < 1)\n")
      cat("=================================================\n")
      cat("Maximum value:", round(max_val, 4), "\n")
      cat("Condition satisfied?", max_val < 1, "\n\n")

      cat("Smoothing Spline selected parameters\n")
      cat("====================================\n")

      smoothing_tbl <- data.frame(
        parameter = c("c", "alpha", "beta"),
        spar = object@smoothing_used$spar,
        method = if (is.null(object@smoothing_used$methods)) {
          rep(NA_character_, 3)
        } else {
          object@smoothing_used$methods
        },
        nknots = object@smoothing_used$nknots
      )

      print(smoothing_tbl, row.names = FALSE)

      cat("Performance metrics: RMSE and MAE \n")
      cat("=================================================\n")
      cat("Nonparametric Method:", object@method, "\n\n")

      res <- eval_model(object@tvparameter[,2:4], object@fitted.values, na.rm = TRUE)
      print_eval_model(res)
    }

    if (identical(object@method, "B-Spline")) {

      max_val <- unname(object@Check_condition["Max."])

      cat("Stationarity condition (alpha(u) + beta(u) < 1)\n")
      cat("=================================================\n")
      cat("Maximum value:", round(max_val, 4), "\n")
      cat("Condition satisfied?", max_val < 1, "\n\n")

      cat("B-Spline selected parameters\n")
      cat("============================\n")

        bs_tbl <- data.frame(
          parameter = c("c", "alpha", "beta"),
          degree = object@best_degree,
          df = object@best_df
        )
        print(bs_tbl, row.names = FALSE)
      # }

      cat("Performance metrics: RMSE and MAE \n")
      cat("=================================================\n")
      cat("Nonparametric Method:", object@method, "\n\n")

      res <- eval_model(object@tvparameter[,2:4], object@fitted.values, na.rm = TRUE)
      print_eval_model(res)
    }

    return(invisible(object))
  }
}
)

