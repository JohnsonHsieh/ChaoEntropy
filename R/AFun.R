AFun <-
function(x) {
  x <- x[x > 0]
  n <- sum(x)
  f <- function(i) { sum(x == i) }
  f1 <- f(1); f2 <- f(2)
  if (f2 > 0) {
    A <- 2 * f2 / ((n - 1) * f1 + 2 * f2)
  } else if (f2 == 0 & f1 != 0) {
    A <- 2 * (f2 + 1) / ((n - 1) * (f1 - 1)+ 2 * (f2 + 1))
  } else {
    A <- 1
  }
  return(A)
}
