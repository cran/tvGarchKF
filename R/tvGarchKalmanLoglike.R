#' @title Models tv-Garch Filter Kalman LogLikehood.
#' @description It is the function to use in the process to fit coefficients.
#' @param x A vector of coefficents to fit.
#' @param data A univariate time series, where could be an object of class "ts" or "zoo".
#' @param c A vector for the coefficients for the function \eqn{c(u)} in the model.
#' @param alpha A vector for the coefficients for the function \eqn{\alpha(u)} in the model.
#' @param beta A vector for the coefficients for the function \eqn{\beta(u)} in the model.
#' @param nsample A value of data length.
#' @param type A character type=c("polynomial","NoLineal","trigonometric"), which indicates what type of function will be assigned to each tv-parameters.
#' @param exponentes A vector for the exponents for the "NoLinear" case.
#' @param trig A character parameter indicating the trigonometric function to be used, "cos" or "sin".
#' @param arg This is the value of the argument to be used in the trigonometric function.
#' @param predict A numeric value that indicates the h-steps of volatility prediction.
#' @param trace.log Variable to print names of coefficients.
#' @return Value of loglike in model.
#' @importFrom stats na.exclude
#' @export
tvGarchKalmanLoglike <- function(x, data, c, alpha, beta, nsample=length(data), type=c("polynomial","NoLineal","trigonometric"), exponentes = NULL, trig = NULL, arg = NULL, predict, trace.log=FALSE){

  pred<-rep(NA,predict)
  data<-c(data,pred)

  nsample <- length(data)
  u <- (1:nsample)/nsample
  params <- tvCoef(x=x, data, c, alpha, beta, type, exponentes, trig, arg)
  c <- params$c
  a <- params$a
  b <- params$b

  for(i in 1:nsample){
    if(any(is.nan(a[i])) || any(is.nan(b[i])) || any(is.nan(c[i]))) {
      stop("Se produjo un error al estimar los parametros")
    }
  }

  for(i in 1:nsample){
    if(trace.log) {
      cat("a[i]:",a[i],"|b[i]:",b[i],"|c:",c,"\n")
    }
  }

  for(i in 1:nsample){
    if ((c[i]<=0) & (a[i]<=0) & (b[i]<=0) & ((b[i]+a[i])>=1) ) {
      return(100000000.)
    }
  }

  for (i in 1:nsample) {
    if ((b[i]+a[i])>=1){
      return(100000000.)
    }
  }

  y <- NULL

  for(i in 1:nsample){
    y[i] <- ifelse (is.na(data[i]), 100000000., data[i]^2-c[i]/(1-a[i]-b[i]))
  }

  P <-rep(0, nsample) # Varianza del error de estimaci?n
  X <-rep(0, nsample)   # Representa el estado (alfa)
  F.m <-rep(0, nsample)
  K <-rep(0, nsample)

  # Inicializo P_1
  P[1]<-(a[1]^2)/(1-(a[1]+b[1])^2)

  if (length(c) > 1){
    resultado <- .C("tvGarch1ab2", as.integer(nsample), as.double(a), as.double(b), as.double(c), as.double(P), as.double(X), as.double(F.m), as.double(K), as.double(y))
  }else{
    resultado <- .C("tvGarch1ab", as.integer(nsample), as.double(a), as.double(b), as.double(c), as.double(P), as.double(X), as.double(F.m), as.double(K), as.double(y))
  }

  X <- resultado[[6]]
  F.m <- resultado[[7]]
  K <- resultado[[8]]
  ####
  e.sig2=X+(c/(1-a-b))


  if( min(e.sig2) <=0 | min(F.m) <=0 ){loglik=100000000}
  else{
    loglik <- sum(na.exclude((sqrt(F.m)*data^2) / e.sig2 + log(e.sig2) - 0.5*log(F.m) ))
    attr(loglik,"sigma")<-e.sig2
    loglik}
}
