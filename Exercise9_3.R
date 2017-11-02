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
nnlike2 <- function(p,x,y){
  B0 = p[1]
  B1 = p[2]
  B2 = p[3]
  sigma = exp(p[4])
  expected = B0+B1*x + B2*x^2
  nll= -sum(dnorm(x=y, mean=expected, sd = sigma, log = TRUE))
  return(nll)
}

params0 <- c(200,100)
params1 <- c(200,15,50)
params2 <- c(200,10,-0.2,1)

fit0 = optim(par = params0, fn=nnlike0, x=leafDecomp$Ms, y = leafDecomp$decomp)
fit1 = optim(par = params1, fn=nnlike1, x=leafDecomp$Ms, y = leafDecomp$decomp)
fit2 = optim(par = params2, fn=nnlike2, x=leafDecomp$Ms, y = leafDecomp$decomp)

fit0$par
fit1$par
fit2$par

test0 <- pchisq(q = abs((fit0$value-fit1$value)), df = 1, lower.tail = FALSE)
test1 <- pchisq(q = abs((fit0$value-fit2$value)), df = 2, lower.tail = FALSE)
test2 <- pchisq(q = abs((fit1$value-fit2$value)), df = 1, lower.tail = FALSE)

test0
test1
test2