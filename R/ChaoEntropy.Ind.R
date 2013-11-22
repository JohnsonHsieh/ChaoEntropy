ChaoEntropy.Ind <-
function(data, method = c("all", "Chao", "ChaoShen", "Grassberger", 
                                             "Jackknife", "Zhang", "Observed"), 
                            B = 200, conf = 0.95, se = TRUE) {
  method <- match.arg(method)
  if (method == "all") {
    a <- entropy_Chao(data, B, conf)
    b <- entropy_ChaoShen(data, B, conf)
    c <- entropy_Grassberger(data, B, conf)
    d <- entropy_Jackknife(data, B, conf)
    e <- entropy_Zhang(data, B, conf)
    f <- entropy_Observed(data, B, conf)
    out <- rbind(a, b, c, d, e, f)
  }
  if (method == "Chao")
    out <- entropy_Chao(data, B, conf)
  if (method == "ChaoShen")
    out <- entropy_ChaoShen(data, B, conf)
  if (method == "Grassberger")
    out <- entropy_Grassberger(data, B, conf)
  if (method == "Jackknife") 
    out <- entropy_Jackknife(data, B, conf)
  if (method == "Zhang")
    out <- entropy_Zhang(data, B, conf)
  if (method == "Observed")
    out <- entropy_Observed(data, B, conf)
  
  if (se == FALSE) 
    out <- data.frame(Estimator = out[, 1], row.names = rownames(out))
  
  return(out)
}
