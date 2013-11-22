Count2Inci <-
function(count.data, t) {
  data <- c(t, rep(count.data[, 1], count.data[, 2]))
  return(data)
}
