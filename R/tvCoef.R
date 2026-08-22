tvCoef <- function(x, data, c, alpha, beta, type, exponentes, trig, arg, S, N, flag = TRUE){
  nsample <- length(data)
  aux1 <- length(c)
  aux2 <- length(alpha)
  aux3 <- length(beta)
  if (aux1<=1){
    c = x[1]
  }else{
    c = c()
    for(i in 1:aux1){
      c[i] = x[i]
    }
  }
  if (aux2<=1){
    alpha = x[1+aux1]
  }else{
    alpha = c()
    for(i in 1:aux2){
      alpha[i] = x[i+aux1]
    }
  }
  if (aux3<=1){
    beta = x[1+aux1+aux2]
  }else{
    beta = c()
    for(i in 1:aux3){
      beta[i] = x[i+aux2+aux1]
    }
  }
  auxA <- list(c=c, alpha=alpha, beta=beta)
  auxB <- list()
  #auxC <- list()
  if(flag==FALSE){
    for(i in 1:length(type)){
      if(type[i] == "polynomial"){
        auxC <- polynomial_param(nsample, coeficientes = unlist(auxA[i]), S, N)
      }
      if(type[i] == "NoLineal"){
        auxC <- NoLineal_param(nsample, coeficientes = unlist(auxA[i]), exponentes = exponentes, S, N)
      }
      if(type[i] == "trigonometric"){
        auxC <- trigonometric_param(nsample, coeficientes = unlist(auxA[i]), trig = trig, arg = arg, S, N)
      }
      auxB <- c(auxB,list(auxC))
    }
  }else{
    for(i in 1:length(type)){
      if(type[i] == "polynomial"){
        auxC <- polynomial(nsample, coeficientes = unlist(auxA[i]))
      }
      if(type[i] == "NoLineal"){
        exponent <- exponentes[1:length(unlist(auxA[i]))]
        auxC <- NoLineal(nsample, coeficientes = unlist(auxA[i]), exponentes = exponent)
        exponentes <- exponentes[(length(exponent)+1):length(exponentes)]
      }
      if(type[i] == "trigonometric"){
        auxC <- trigonometric(nsample, coeficientes = unlist(auxA[i]), trig = trig, arg = arg)
      }
      auxB <- c(auxB,list(auxC))
    }
  }

  c <- unlist(auxB[1])
  a <- unlist(auxB[2])
  b <- unlist(auxB[3])
  return(list(a=a,b=b,c=c))
}
