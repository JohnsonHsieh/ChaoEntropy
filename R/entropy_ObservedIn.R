entropy_ObservedIn <-
function(data, B = 200, conf = 0.95) {
  t <- data[1]
  y <- data[-1]
  est <- ObservedEstFunIn(y, t)
  se <- BootstrapFunIn(y, t, B, ObservedEstFunIn)
  z <- qnorm(1 - (1 - conf)/2)
  CI <- c(max(est - z * se, 0), est + z * se)
  out <- matrix(c(est, se, CI), nrow = 1)
  rownames(out) <- c("Observed entropy")
  colnames(out) <- c("Estimator", "Bootstrap s.e.",
                     paste(conf*100, "% Lower"), paste(conf*100, "% Upper"))
  return(out) 
}
