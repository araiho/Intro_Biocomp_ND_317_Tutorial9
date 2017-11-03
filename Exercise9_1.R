#### Exercise 9 Question 1 ####
## Import Data ##
ponzr1 <- read.csv("~/Desktop/data-shell/Intro_Biocomp_ND_317_Tutorial9/ponzr1.csv",sep = ',',header=TRUE)

## Custom Function that will determine 
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

## Parameters can be found by plotting the data ##
params1 <- c(2000,-1000,10)
params0 <- c(2000,10)

fit0 = optim(par = params0, fn=nnlike0, y=ponzr1[(ponzr1$mutation%in%c("WT","M124K")),]$ponzr1Counts, x = 0)
fit1 = optim(par = params1, fn=nnlike1, y=ponzr1[(ponzr1$mutation%in%c("WT","M124K")),]$ponzr1Counts, x = 1)
fit2 = optim(par = params1, fn=nnlike1, y=ponzr1[(ponzr1$mutation%in%c("WT","V456D")),]$ponzr1Counts, x = 1)
fit3 = optim(par = params1, fn=nnlike1, y=ponzr1[(ponzr1$mutation%in%c("WT","I213N")),]$ponzr1Counts, x = 1)

M124K <- pchisq(q = 2*abs((fit0$value-fit1$value)), df = 1, lower.tail = FALSE)
V456D <- pchisq(q = 2*abs((fit0$value-fit2$value)), df = 1, lower.tail = FALSE)
I213N <- pchisq(q = 2*abs((fit0$value-fit3$value)), df = 1, lower.tail = FALSE)

M124K
V456D
I213N
