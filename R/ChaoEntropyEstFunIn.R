ChaoEntropyEstFunIn <-
function(y, t) {
  y <- y[y > 0]
  U <- sum(y)
  Q1 <- sum(y == 1)
  tmp <- A.C.Q0FunIn(y, t)
  A <- tmp[1]
  temp1 <- sum(y / t * (digamma(t) - digamma(y)))
  if (A == 1) {
    temp2 <- 0
  } else {
    r <- 1 : (t-1)
    temp2 <- Q1 / t * (1 - A)^(1-t) * (-log(A) - sum(1 / r * (1 - A)^r))
  }
  est <- t / U * (temp1 + temp2) + log(sum(y) / t)
  return(est)
}
