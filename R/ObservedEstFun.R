ObservedEstFun <-
function(x) { 
  #   entropy(x, method = "ML")
  x <- x[x > 0]
  n <- sum(x)
  pi <- x / n
  est <- -sum(pi * log(pi))
  return(est)
}
