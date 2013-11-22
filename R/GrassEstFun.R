GrassEstFun <-
function(dat) { 
  n <- sum(dat)
  p.hat <- dat/n
  a = 0
  for(i in 1:length(dat)){
    if(p.hat[i] > 0){
      integrand <- function(x) { x^(dat[i] - 1) / (1 + x) }
      a <- a + p.hat[i] * (digamma(n) - digamma(dat[i]) - (-1)^dat[i]
                           * (integrate(integrand, lower = 0, upper = 1))$value)
    }
  }
  return(a)
}
