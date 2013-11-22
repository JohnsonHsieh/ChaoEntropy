BootstrapFun <-
function(x, B, FunName) {
  n <- sum(x)
  f <- function(i) { sum(x == i) }
  f1 = f(1); f2 = f(2)
  tmp <- Candf0Fun(f1, f2, n)
  Chat <- tmp[1] ; f0 <- tmp[2]
  lambda <- (1 - Chat) / sum(x / n * (1 - x / n)^n)
  pi <- x / n * (1 - lambda * (1 - x /n)^n)
  pi.star <- c(pi, rep((1 - Chat) / f0, f0))
  #   set.seed(1234)
  X <- rmultinom(B, n, pi.star)
  se <- sd(apply(X, 2, function(x) FunName(x)))
  return(se)
}
