#Housekeeping
rm(list=ls())

library(stringr)


#Loading Data
data = read.csv(file = "leafDecomp.csv", header = TRUE, sep = ",")

NullNllike<-function(p,y){
  B0=p[1]
  sigma=p[2]
  
  expected=B0
  
  nll=-sum(dnorm(x=y, mean=expected, sd=sigma, log = TRUE))
  return((nll))}

LinNllike<-function(p,x,y){
  B0=p[1]
  B1=p[2]
  sigma=p[3]
  
  expected=B0+B1*x
  
  nll=-sum(dnorm(x=y,mean=expected,sd=sigma,log=TRUE))
  return(nll)
}

QuadNllike<-function(p,x,y){
  B0=p[1]
  B1=p[2]
  B2=p[3]
  sigma=p[4]
  
  expected=B0+B1*x+B2*x^2
  
  nll=-sum(dnorm(x=y,mean=expected,sd=sigma,log=TRUE))
  return(nll)
}

NullGuess=c(600,164)
NullFit=optim(NullGuess, NullNllike, y=data$decomp)

LinGuess=c(318,6.3,54)
LinFit=optim(LinGuess, LinNllike, x=data$Ms, y=data$decomp)

QuadGuess=c(180,15.7,-0.11,10.7)
QuadFit=optim(QuadGuess, QuadNllike, x=data$Ms, y=data$decomp)

DLinNull=-2*(LinFit$value-NullFit$value)
pLinNull=pchisq(q=DLinNull, df=1, lower.tail=FALSE)

DQuadLin = 2*(LinFit$value-QuadFit$value)
pQuadLin=pchisq(q=DQuadLin, df=1, lower.tail=FALSE)

DQuadNull = -2*(QuadFit$value - NullFit$value)
pQuadNull=pchisq(q=DQuadNull,df=2,lower.tail=FALSE)
