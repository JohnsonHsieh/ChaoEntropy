BootstrapFunIn <-
function(y, t, B, FunName) {
  U <- sum(y)
  tmp <- A.C.Q0FunIn(y, t)
  Chat <- tmp[2]
  Q0 <- tmp[3]
  tau <- U / t * (1 - Chat) / sum(y / t * (1 - y / t)^t)
  pi <- y / t * (1 - tau * (1 - y / t)^t)
  pi.star <- c(pi, rep(U / t * (1 - Chat) / Q0, Q0))
  #   set.seed(456)
  y1 <- matrix(rbinom(length(pi.star) * B, t, pi.star), ncol = B)
  se <- sd(apply(y1, 2, function(y2) FunName(y2, t)))
  return(se)
}
