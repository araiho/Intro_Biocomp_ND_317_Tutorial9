nnlike0 <- function(p,x,y){
  B0 = p[1]
  sigma = exp(p[2])
  expected = B0
  nll= -sum(dnorm(x=y, mean=expected, sd = sigma, log = TRUE))
  return(nll)
}

nnlike1 <- function(p,x,y){
  B0 = p[1]
  B1 = p[2]
  sigma = exp(p[3])
  expected = B0+B1*x
nll1= -sum(dnorm(x=y, mean=expected, sd = sigma, log = TRUE))
return(nll1)
}


params1 <- c(0.5,0.5,0.5)
params0 <- c(0.5,0.5)

fit0 = optim(par = params0, fn=nnlike0, y=subset(ponzr1, mutation == "WT")$ponzr1Counts, x = 0)
fit1 = optim(par = params1, fn=nnlike1, y=subset(ponzr1, mutation == "M124K")$ponzr1Counts, x = 1)
fit2 = optim(par = params1, fn=nnlike1, y=subset(ponzr1, mutation == "V456D")$ponzr1Counts, x = 1)
fit3 = optim(par = params1, fn=nnlike1, y=subset(ponzr1, mutation == "I213N")$ponzr1Counts, x = 1)


M124K <- pchisq(q = abs((fit0$value-fit1$value)), df = 1, lower.tail = FALSE)
V456D <- pchisq(q = abs((fit0$value-fit2$value)), df = 1, lower.tail = FALSE)
I213N <- pchisq(q = abs((fit0$value-fit3$value)), df = 1, lower.tail = FALSE)

M124K
V456D
I213N

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
