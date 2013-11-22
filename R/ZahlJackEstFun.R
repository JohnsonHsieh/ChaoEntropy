ZahlJackEstFun <-
function(x) { 
  x <- x[x > 0]
  n <- sum(x)
  sumjk <- -1 * sum((n - x) * x / (n - 1) * log(x / (n - 1))) - 
    sum((x[x > 1] - 1) * (x[x > 1]) / (n - 1) * log((x[x > 1] - 1) / (n - 1)))
  n * ObservedEstFun(x) - (n - 1) / n * sumjk
}
