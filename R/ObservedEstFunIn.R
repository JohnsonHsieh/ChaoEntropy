ObservedEstFunIn <-
function(y, t) {
  y <- y[y > 0]
  U <- sum(y)
  pi <- (y / t) / (U / t)
  est <- -sum(pi * log(pi))
  return(est)
}
