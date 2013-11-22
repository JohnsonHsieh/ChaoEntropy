ChaoEntropyEstFun <-
function(x) {
  x <- x[x > 0]
  n <- sum(x)
  f <- function(i) { sum(x == i) }
  A <- AFun(x)
  temp1 <- sum(x / n * (digamma(n) - digamma(x)))
  if(A == 1){
    temp2 <- 0
  } else {
    l <- 1:(n-1)
    temp2 <- f(1) / n * (1 - A)^(1-n) * (-log(A) - sum(1 / l * (1-A)^l))
  }
  est <- (temp1 + temp2)
  return(est)
}
