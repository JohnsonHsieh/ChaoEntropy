entropy_Jackknife <-
function(data, B = 200, conf = 0.95) {
  if (is.matrix(data) == TRUE) {
    if (ncol(data) != 1 & nrow(data) != 1)
      stop("Error: The data format is wrong.")
    if (ncol(data) == 1) {
      data <- data[, 1]
    } else {
      data <- data[1, ]
    }
  }
  
  est <- ZahlJackEstFun(data)
  se <- BootstrapFun(data, B, ZahlJackEstFun)
  z <- qnorm(1 - (1 - conf)/2)
  CI <- c(max(est - z * se, 0), est + z * se)
  out <- matrix(c(est, se, CI), nrow = 1)
  rownames(out) <- c("Zahl (1977) Jackknife")
  colnames(out) <- c("Estimator", "Bootstrap s.e.",
                     paste(conf*100, "% Lower"), paste(conf*100, "% Upper"))
  return(out)
}
