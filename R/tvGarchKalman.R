tvGarchKalman <- function(x, data, c, alpha, beta, type, exponentes = NULL, trig = NULL, arg = NULL, predict = NULL){
  nsample=length(data)
  x <- x
  tvGarchKalmanLoglike(x=x, data = data, c = c, alpha = alpha, beta = beta, type = type, exponentes = exponentes, trig = trig, arg = arg, predict)
}
