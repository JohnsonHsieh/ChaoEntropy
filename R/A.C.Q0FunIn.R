A.C.Q0FunIn <-
function(y, t) {
  Q1 <- sum(y == 1)
  Q2 <- sum(y == 2)
  U <- sum(y)
  if (Q2 > 0) {
    A <- 2 * Q2 / ((t - 1) * Q1 + 2 * Q2)
    C <- 1 - Q1 / U * ((t - 1) * Q1 / ((t - 1) * Q1 + 2 * Q2))
    Q0 <- (t - 1) / t * Q1^2 / (2 * Q2)
  } else if (Q2 == 0 & Q1 != 0) {
    A <- 2 / ((t - 1) * (Q1 - 1) + 2)
    C <- 1 - Q1 / U * ((t - 1) * (Q1 - 1) / ((t - 1) * (Q1 - 1) + 2))
    Q0 <- (t - 1) / t * Q1 * (Q1 - 1) / 2
  } else {
    A <- 1
    C <- 1
    Q0 <- (t - 1) / t * Q1 * (Q1 - 1) / 2
  }
  Q0 <- ceiling(Q0)
  return(c(A, C, Q0))
}
