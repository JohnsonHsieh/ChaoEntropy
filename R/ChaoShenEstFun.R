ChaoShenEstFun <-
function(x){ 
  #   entropy(x, method = "CS")
  x <- x[x > 0]
  n <- sum(x)
  f <- function(i) { sum(x == i) }
  f1 = f(1); f2 = f(2)
  Chat <- Candf0Fun(f1, f2, n)[1]
  adjpi <- Chat * x / n
  est <- - sum(adjpi * log(adjpi) / (1 - (1 - adjpi)^n))
  return(est)
}
