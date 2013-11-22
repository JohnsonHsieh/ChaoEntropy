Matrix2Inci <-
function(mat.data) {
  data <- c(length(mat.data[1, ]), apply(mat.data, 1, sum))
  return(data)
}
