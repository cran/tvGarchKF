polynomial <- function(n, coeficientes=numeric()){
  aux1 <- length(coeficientes)
  u <- (1:n)/n

  eta <- matrix(NA, ncol = aux1, nrow = n)
  for(j in 1:(aux1)){
    eta[,j] <- coeficientes[j]*u^(j-1)
  }
  result <- numeric(n)
  for (i in 1:n) {
    result[i] <- sum(eta[i, ], na.rm = TRUE)
  }
  return(result)
}

NoLineal <- function(n, coeficientes=numeric(), exponentes=numeric()){
  aux1 <- length(coeficientes)
  u <- (1:n)/n

  eta <- matrix(NA, ncol = aux1, nrow = n)
  for(j in 1:(aux1)){
    eta[,j] <- coeficientes[j]*u^(exponentes[j])
  }
  result <- numeric(n)
  for (i in 1:n) {
    result[i] <- sum(eta[i, ], na.rm = TRUE)
  }
  return(result)
}

trigonometric <- function(n, coeficientes = numeric(), trig = c("sin","cos"), arg = "u"){

  aux1 <- length(coeficientes)
  u <- (1:n)/n

  eta <- matrix(NA, ncol = aux1, nrow = n)
  eta[,1] <- coeficientes[1]
  arg <- parse(text = arg)

  for(j in 2:(aux1)){
    if (trig == "cos"){
      eta[,j] <- coeficientes[j]*cos(eval(arg, envir = list(u = u)))
    }else{
      eta[,j] <- coeficientes[j]*sin(eval(arg, envir = list(u = u)))
    }
  }

  result <- numeric(n)
  for (i in 1:n) {
    result[i] <- sum(eta[i, ], na.rm = TRUE)
  }
  return(result)
}

polynomial_param <- function(n, coeficientes=numeric(), S, N){
  aux1 <- length(coeficientes)
  T. = n
  S  = S
  N  = N
  M  = trunc((T.-N)/S+1)
  t = S * (1:M - 1) + N/2
  u. = t/T.

  eta <- matrix(NA, ncol = aux1, nrow = length(u.))
  for(j in 1:(aux1)){
    eta[,j] <- coeficientes[j]*u.^(j-1)
  }
  result <- numeric(length(u.))
  for (i in 1:length(u.)) {
    result[i] <- sum(eta[i, ], na.rm = TRUE)
  }
  return(result)
}

NoLineal_param <- function(n, coeficientes=numeric(), exponentes=numeric(), S, N){
  aux1 <- length(coeficientes)
  T. = n
  S  = S
  N  = N
  M  = trunc((T.-N)/S+1)
  t = S * (1:M - 1) + N/2
  u. = t/T.

  eta <- matrix(NA, ncol = aux1, nrow = length(u.))
  for(j in 1:(aux1)){
    eta[,j] <- coeficientes[j]*u.^(exponentes[j])
  }
  result <- numeric(length(u.))
  for (i in 1:length(u.)) {
    result[i] <- sum(eta[i, ], na.rm = TRUE)
  }
  return(result)
}

trigonometric_param <- function(n, coeficientes = numeric(), trig = c("sin","cos"), arg = "u", S, N){
  aux1 <- length(coeficientes)
  T. = n
  S  = S
  N  = N
  M  = trunc((T.-N)/S+1)
  t = S * (1:M - 1) + N/2
  u. = t/T.

  eta <- matrix(NA, ncol = aux1, nrow = length(u.))
  eta[,1] <- coeficientes[1]
  arg <- parse(text = arg)

  for(j in 2:(aux1)){
    if (trig == "cos"){
      eta[,j] <- coeficientes[j]*cos(eval(arg, envir = list(u = u.)))
    }else{
      eta[,j] <- coeficientes[j]*sin(eval(arg, envir = list(u = u.)))
    }
  }

  result <- numeric(length(u.))
  for (i in 1:length(u.)) {
    result[i] <- sum(eta[i, ], na.rm = TRUE)
  }
  return(result)
}

checkInput <- function(type, c, a, b, exponents) {
  if(length(type) != 3) stop("invalid type length")
  if(type[1] != "polynomial" && type[1] != "NoLineal" && type[1] != "trigonometric") stop("Invalid type for c.")
  if(type[2] != "polynomial" && type[2] != "NoLineal" && type[2] != "trigonometric") stop("Invalid type for a.")
  if(type[3] != "polynomial" && type[3] != "NoLineal" && type[3] != "trigonometric") stop("Invalid type for b.")
  sup <- c()
  #if(exponents == 0) exponents <- 1
  if(type[1] == "NoLineal") sup <- c(sup,c)
  if(type[2] == "NoLineal") sup <- c(sup,a)
  if(type[3] == "NoLineal") sup <- c(sup,b)
  if(length(sup) != exponents) stop("Not match length of exponents and length of coefficients No Lineal functions.")
}

tv_function_expression <- function(type, par_name, npar,
                                   trig = NULL, arg = NULL, exponentes = NULL) {

  type <- as.character(type)

  if (type == "polynomial") {

    terms <- paste0(par_name, 0:(npar - 1))

    if (npar >= 2) {
      for (j in 2:npar) {
        power <- j - 1
        if (power == 1) {
          terms[j] <- paste0(terms[j], "u")
        } else {
          terms[j] <- paste0(terms[j], "u^", power)
        }
      }
    }

    return(paste(terms, collapse = " + "))
  }

  if (type == "trigonometric") {

    if (is.null(trig)) trig <- "cos"
    if (is.null(arg)) arg <- "u"

    return(
      paste0(
        par_name, "0 + ",
        par_name, "1", trig, "(", arg, ")"
      )
    )
  }

  if (type == "NoLineal") {

    if (is.null(exponentes)) {
      exponentes <- 0:(npar - 1)
    }

    terms <- character(npar)

    for (j in seq_len(npar)) {
      power <- exponentes[j]

      if (power == 0) {
        terms[j] <- paste0(par_name, j - 1)
      } else if (power == 1) {
        terms[j] <- paste0(par_name, j - 1, "u")
      } else if (power == 0.5) {
        terms[j] <- paste0(par_name, j - 1, "sqrt(u)")
      } else {
        terms[j] <- paste0(par_name, j - 1, "u^", power)
      }
    }

    return(paste(terms, collapse = " + "))
  }
}
