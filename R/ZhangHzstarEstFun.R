ZhangHzstarEstFun <-
function(x) {
  x <- x[x > 0]
  n <- sum(x)
  forward <- sum(x / n * (digamma(n) - digamma(x)))
  k <- 1 : (n-1)
  back <- ObservedEstFun(x) - 
    sum(sapply(k, function(k) 1 / (n * k) * sum(x * (1 - x / n)^k)))
  return(forward + back)
}
