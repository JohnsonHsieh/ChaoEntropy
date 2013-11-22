ChaoEntropy.Sam <-
function(data, method = c("all", "Chao", "Observed"),
                            B = 200, conf = 0.95, se = TRUE) {
  method <- match.arg(method)
  if (method == "all") {
    a <- entropy_ChaoIn(data, B, conf)
    f <- entropy_ObservedIn(data, B, conf)
    out <- rbind(a, f)
  }
  if (method == "Chao")
    out <- entropy_ChaoIn(data, B, conf)
  if (method == "Observed")
    out <- entropy_ObservedIn(data, B, conf)
  if (se == FALSE) 
    out <- data.frame(Estimator = out[, 1], row.names = rownames(out))
  
  return(out)
}
