###Exercise9_Question2###

#Load data
setwd("~/Desktop/data-shell/Exercise9/Intro_Biocomp_ND_317_Tutorial9/")
M <- read.csv("MmarinumGrowth.csv", header=TRUE)

#Look at the plot to check it's shape. It is an exponetial growth curve.
plot(M$S, M$u, xlab = "Concentration", ylab = "Growth Rate")

#Custom Function
nllike <- function(p,S,y){
  umax=p[1]
  Ks=p[2]
  sigma=exp(p[3])
  
  expected=umax*(S/(S+Ks))
    
  
  nll <- -sum(dnorm(M$u,mean=expected,sd=sigma,log=TRUE))
  return(nll)
}


#Estimate parameters by minimizing the negative log likihood.
intialGuess <- c(100,10,0)
fit <- optim(par=intialGuess,fn=nllike,S=M$S,y=M$u)

print(fit)





######Exercise9_Question3######

setwd("~/Desktop/data-shell/Exercise9/Intro_Biocomp_ND_317_Tutorial9/")
L <- read.csv("leafDecomp.csv", header=TRUE)

#Look at the plot to check it's shape. It is an exponetial growth curve.
plot(L$Ms, L$decomp, xlab = "Soil Moisture", ylab = "Decomposition Rate")

#Custom Function 1: Constant Rate (d = a)
nllike1 <- function(p,x,y){
  a=p[1]
  sigma=exp(p[2])
  
  expected=a
  
  nll1 <- -sum(dnorm(x=0,mean=expected,sd=sigma,log=TRUE))
  return(nll1)
}

#Estimate parameters by minimizing the negative log likihood.
intialGuess <- c(200,1)
fit1 <- optim(par=intialGuess,fn=nllike1,x=L$Ms,y=L$decomp)

print(fit1)


#Custom Function 2: Linear Rate (d = a + bMs)
nllike2 <- function(p,x,y){
  b=p[1]
  a=p[2]
  sigma=exp(p[3])
  
  expected=a+b*x
  
  nll2 <- -sum(dnorm(x=y,mean=expected,sd=sigma,log=TRUE))
  return(nll2)
}

#Estimate parameters by minimizing the negative log likihood.
intialGuess <- c(100,200,1)
fit2 <- optim(par=intialGuess,fn=nllike2,x=L$Ms,y=L$decomp)

print(fit2)

#Custom Function 1: Quadratic (d = a + bMs + c(Ms)^2)
nllike3 <- function(p,x,y){
  b=p[1]
  a=p[2]
  c=p[3]
  sigma=exp(p[4])
  
  expected= a + b*x + c*((x)^2)
  
  nll3 <- -sum(dnorm(L$Ms,mean=expected,sd=sigma,log=TRUE))
  return(nll3)
}

#Estimate parameters by minimizing the negative log likihood.
intialGuess <- c(0.1,200,10,1)
fit3 <- optim(par=intialGuess,fn=nllike3,x=L$Ms,y=L$decomp)

print(fit3)

#Chi-squared tests
pchisq(q=D,df=1,lower.tail=FALSE)
