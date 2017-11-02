####Question 2
Monod <- function(p,x,y){
  S = p[1]
  K = p[2]
  sigma = exp(p[3])
  expected = (S*(x/(x+K)))
  mnd= -sum(dnorm(x=y, mean=expected, sd = sigma, log = TRUE))
  return(mnd)
}
params2 <- c(2,2,0)
fit4 <- optim(par = params2, fn=Monod, y=MmarinumGrowth$u, x = MmarinumGrowth$S)
fit4