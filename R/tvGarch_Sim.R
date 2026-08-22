#' @title Generating Simulations using a tv-Garch Model
#' @description Simulate from a tv-Garch(1,1) model, using specific functions for the tv-parameters.
#' @param n This is the value of the length of the series to be simulated.
#' @param c A vector for the coefficients for the function \eqn{c(u)} in the model.
#' @param alpha A vector for the coefficients for the function \eqn{\alpha(u)} in the model.
#' @param beta A vector for the coefficients for the function \eqn{\beta(u)} in the model.
#' @param type A character type=c("polynomial","NoLineal","trigonometric"), which indicates what type of function will be assigned to each tv-parameters.
#' @param exponentes A vector for the exponents for the "NoLinear" case.
#' @param trig A character parameter indicating the trigonometric function to be used, "cos" or "sin".
#' @param arg This is the value of the argument to be used in the trigonometric function.
#' @param seed This value is the seed to simulate series.
#' @return An object of class "zoo" with two components: the first component represents returns "x", while the second component denotes conditional variance "sigma".
#' @importFrom stats na.exclude
#' @importFrom stats rnorm
#' @importFrom withr with_seed
#' @export
#' @examples
#' ## Simulate from a tv-GARCH(1,1) model lineal:
#' alpha_sim <- c(0.2, 0.2)
#' beta_sim <- c(0.45, 0.5, -0.85)
#' type_sim <- c("polynomial","polynomial","polynomial")
#' Sim1 <- tvGarch_Sim(n = 6000, c = 0.1, alpha = alpha_sim, beta = beta_sim, type = type_sim)
#' plot(Sim1[,1], type="l", main="Simulated tvGARCH(1, 1) process",
#'     ylim=c(-max(Sim1[,2]), max(Sim1[,2])))
#' lines(Sim1[,2], type="l", col="red")
#' legend("topright",legend=c("tvGARCH(1,1)",expression(sigma(u))),
#'       col=c("black","red"),lty=1,bty="n",lwd=1)
#' ## Simulate from a tv-GARCH(1,1) model non linear:
#' alpha_sim2 <- c(0.75, 0.08)
#' beta_sim2 <- c(0.05, 0.03, 0.06)
#' type_sim2 <- c("polynomial","polynomial","NoLineal")
#' expo <- c(0, 1, 1/2)
#' Sim2<-tvGarch_Sim(n=6000,c=0.05,alpha=alpha_sim2,beta=beta_sim2,type=type_sim2,exponentes=expo)
#' plot(Sim2[,1], type="l", main="Simulated tvGARCH(1, 1) process",
#'      ylim=c(-max(Sim2[,2]),max(Sim2[,2])))
#' lines(Sim2[,2], type="l", col="red")
#' legend("topright",legend=c("tvGARCH(1,1)",expression(sigma(u))),
#'        col=c("black","red"),lty=1,bty="n",lwd=1)
tvGarch_Sim <- function(n, c, alpha, beta, type = c("polynomial", "NoLineal","trigonometric"),
                              exponentes = NULL, trig = NULL, arg = NULL, seed = NULL){

  if(!is.null(seed)){
    e <- withr::with_seed(seed,{rnorm(n)})
  }else{
    e <- rnorm(n)
    }
  x<-rep(0, n)
  sigma<-rep(0,n)
  u<-(1:n)/n

  aux <- list(c=c, alpha=alpha, beta=beta)
  aux3 <- list()

  for(i in 1:length(type)){
    if(type[i] == "polynomial"){
      aux2 <- polynomial(n, coeficientes = unlist(aux[i]))
    }
    if(type[i] == "NoLineal"){
      aux2 <- NoLineal(n, coeficientes = unlist(aux[i]), exponentes = exponentes)
    }
    if(type[i] == "trigonometric"){
      aux2 <- trigonometric(n, coeficientes = unlist(aux[i]), trig = trig, arg = arg)
    }
    aux2 <- list(aux2)
    aux3 <- c(aux3,aux2)
  }

  c <- unlist(aux3[1])
  a <- unlist(aux3[2])
  b <- unlist(aux3[3])

  for(i in 2:n){
    sigma[i]<- c[i]+a[i]*x[i-1]**2+b[i]*sigma[i-1]
    x[i]<-e[i]*sqrt(sigma[i])
  }
  Z<-cbind(x,sigma)
  return(Z)
}
