setwd("C:/Users/DAVIS/Desktop/shell-novice-data/exe9/Intro_Biocomp_ND_317_Tutorial9/")

##########QUESTION 1 ################
#load ponzr csv data file
Ponzr <- read.csv("ponzr1.csv", header=T)

#make custom likelihood function that specifies model structure (parameters, observations)
#-'unpack parameters'
#-assign parameters to expected values (make a place to feed in variables)
#-dnorm(x,mean,sd)
like1 <- function(p,x,y){
  B0=p[1]
  B1=p[2]
  sigma = exp(p[3])
  expected = B0
  
  nll = -sum(dnorm(y, mean=expected, sd=sigma, log=TRUE))
  return(nll)
}

#make second likelihood model
like2 <- function(p,x,y){
  B0=p[1]
  B1=p[2]
  sigma = exp(p[3])
  expected = B0+B1*x
  
  nll = -sum(dnorm(y, mean=expected, sd=sigma, log=TRUE))
  return(nll)
}

#change all WT to 0 and each mutation to 1 and split them up
Ponzr1 <- Ponzr[which(Ponzr$mutation=="WT"),]
Ponzr1[,1] <- 0
Ponzr2 <- Ponzr[which(Ponzr$mutation=="M124K"),]
Ponzr2[,1] <- 1
Ponzr3 <- Ponzr[which(Ponzr$mutation=="V456D"),]
Ponzr3[,1] <- 1
Ponzr4 <- Ponzr[which(Ponzr$mutation=="I213N"),]
Ponzr4[,1] <- 1

Ponzr1_2 <- rbind(Ponzr1,Ponzr2)
Ponzr1_3 <- rbind(Ponzr1,Ponzr3)
Ponzr1_4 <- rbind(Ponzr1,Ponzr4)


#optim function to look for maximum likelihood of our model
#estimate parameters by minimizing the NLL
#create a vector with initial guesses
#minimize log likelihood

#wt vs M124K
Guess = c(1,1,1)
fit1=optim(Guess, like1, x=Ponzr1$mutation, y=Ponzr1$ponzr1Counts)
fit2=optim(Guess, like2, x=Ponzr2$mutation, y=Ponzr2$ponzr1Counts)
fit1$value
fit2$value
D = 2*(fit1$value-fit2$value)
pchisq(D, df=1, lower.tail=F)
##Mutation M124K p-value=0.72, no effect of treatment

#wt vs. V456D
Guess = c(1,1,1)
fit1=optim(Guess, like1, x=Ponzr1$mutation, y=Ponzr1$ponzr1Counts)
fit2=optim(Guess, like2, x=Ponzr3$mutation, y=Ponzr3$ponzr1Counts)
fit1$value
fit2$value
D = 2*(fit1$value-fit2$value)
pchisq(D, df=1, lower.tail=F)
##Mutation V456D p-value=0.0000056 effect of treatment

#wt vs.I213N
Guess = c(1,1,1)
fit1=optim(Guess, like1, x=Ponzr1$mutation, y=Ponzr1$ponzr1Counts)
fit2=optim(Guess, like2, x=Ponzr4$mutation, y=Ponzr4$ponzr1Counts)
fit1$value
fit2$value
D = 2*(fit1$value-fit2$value)
pchisq(D, df=1, lower.tail=F)
##Mutation I213N p-value 0.88 no effect of treatment


##Mutation M124K p-value=0.72, no effect of treatment
##Mutation V456D p-value=0.0000056 effect of treatment
##Mutation I213N p-value 0.88 no effect of treatment
##Therefore, the V456D mutation significantly reduced the expression of ponzr1


##########QUESTION 2 ################
#load csv file for mmarinum
Mmarinum <- read.csv("MmarinumGrowth.csv", header=T)

#make custom likelihood function
nllike <- function(p,x,y){
  B0=p[1]
  B1=p[2]
  sigma = exp(p[3])
  expected = B0*(x/(B1+x))
  
  nll = -sum(dnorm(x=y, mean=expected, sd=sigma, log=TRUE))
  return(nll)
}

#optim function to find the max growth rate and half-saturatation constant
Guess = c(1,1,1)
fit=optim(Guess, nllike, x=Mmarinum$S, y=Mmarinum$u)
print(cbind("umax is ", fit$par[1], " and Ks is ", fit$par[2]))

##max growth rate (umax)=1.46
##half sat constant=42.6
##sigma=0.04



#############QUESTION 3 #############
#load data for decomposition of leaves
leafDecomp <- read.csv("leafDecomp.csv", header = T)

#create custom functions for the three models:
#constant fit model(null model)
constant <- function(p,x,y){
  sigma = exp(p[1])
  B0=p[2]
  
  expected = B0
  
  nll = -sum(dnorm(x=y, mean=expected, sd=sigma, log=TRUE))
  return(nll)
}

#linear model
linear <- function(p,x,y){
  sigma = exp(p[1])
  B0=p[2]
  B1=p[3]
  
  expected = B0 + B1*x
  
  nll = -sum(dnorm(x=y, mean=expected, sd=sigma, log=TRUE))
  return(nll)
}

#quadratic model
quadratic <- function(p,x,y){
  sigma = exp(p[1])
  B0=p[2]
  B1=p[3]
  B2=p[4]
  
  expected = B0 + B1*x + B2*x*x
  
  nll = -sum(dnorm(x=y, mean=expected, sd=sigma, log=TRUE))
  return(nll)
}

#want to minimize parameters so we find max
#likelihood to capture shape but with as few parameters
#start with initial guess which gives the correct number of parameters
constantGuess = c(1,1)
linearGuess = c(1,1,1)
quadraticGuess = c(1,1,1,1)

#optimize each likelihood function with different number of parameters
constantResult=optim(constantGuess, constant, x=leafDecomp$Ms, y=leafDecomp$decomp)
linearResult=optim(linearGuess, linear, x=leafDecomp$Ms, y=leafDecomp$decomp)
quadraticResult=optim(quadraticGuess, quadratic, x=leafDecomp$Ms, y=leafDecomp$decomp)

constantResult$value
linearResult$value
quadraticResult$value

constantResult$par
linearResult$par
quadraticResult$par

#do t-tests on all combinations of min NLL for 3 models:
D1_2 = 2*(constantResult$value - linearResult$value)
D1_3 = 2*(constantResult$value - quadraticResult$value)
D2_3 = 2*(quadraticResult$value - linearResult$value)

#set df to 1 or 2 depending on difference in parameters between two models
pchisq(D1_2, df=1, lower.tail=F)
pchisq(D1_3, df=2, lower.tail=F)
pchisq(D2_3, df=1, lower.tail=F)

#p-values for likelihood ratio tests are all about 0
#contant fit B0=589.7, sigma=164
#linear fit B0=318, B1=6.3, sigma=54
#quadratic fit B0=180, B1=15.7, B2=-0.11, signma=10.7

##quadratic model is the best, linear model is second best, null model is the worst