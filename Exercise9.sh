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





###Exercise9_Question3###

setwd("~/Desktop/data-shell/Exercise9/Intro_Biocomp_ND_317_Tutorial9/")
L <- read.csv("leafDecomp.csv", header=TRUE)

#Look at the plot to check it's shape. It is an exponetial growth curve.
plot(L$Ms, L$decomp, xlab = "Soil Moisture", ylab = "Decomposition Rate")

#Custom Function 1: Constant Rate (d = a)
nllike1 <- function(p,x,y){
  Bo=p[1]
  sigma=exp(p[2])
  
  expected=Bo
  
  nll1 <- -sum(dnorm(x=0,mean=expected,sd=sigma,log=TRUE))
  return(nll1)
}

#Estimate parameters by minimizing the negative log likihood.
intialGuess <- c(475,1)
fit1 <- optim(par=intialGuess,fn=nllike1,x=L$Ms,y=L$decomp)

print(fit1)


#Custom Function 2: Linear Rate (d = a + bMs)
nllike2 <- function(p,x,y){
  Bo=p[1]
  B1=p[2]
  sigma=exp(p[3])
  
  expected=Bo+B1*x
  
  nll2 <- -sum(dnorm(x=y,mean=expected,sd=sigma,log=TRUE))
  return(nll2)
}

#Estimate parameters by minimizing the negative log likihood.
intialGuess <- c(100,200,1)
fit2 <- optim(par=intialGuess,fn=nllike2,x=L$Ms,y=L$decomp)

print(fit2)

#Custom Function 1: Quadratic (d = a + bMs + c(Ms)^2)
nllike3 <- function(p,x,y){
  Bo=p[1]
  B1=p[2]
  B2=p[3]
  sigma=exp(p[4])
  
  expected=Bo + B1*x + B2*x^2
  
  nll3 <- -sum(dnorm(L$Ms,mean=expected,sd=sigma,log=TRUE))
  return(nll3)
}

#Estimate parameters by minimizing the negative log likihood.
intialGuess <- c(200,10,-0.2,1)
fit3 <- optim(par=intialGuess,fn=nllike3,x=L$Ms,y=L$decomp)

print(fit3)

#Chi-squared tests

#Compare 1 vs 2
Q <-2*(fit1$value-fit2$value)
pchisq(q=Q,df=1,lower.tail=FALSE)

#Compare 2 vs 3
D <-2*(fit2$value-fit3$value)
pchisq(q=D,df=1,lower.tail=FALSE)

#Compare 1 vs 3
H <-2*(fit2$value-fit3$value)
pchisq(q=H,df=2,lower.tail=FALSE)
