#### Pseudo-code ####

#1. Load Data

#2 Custom function(params, observations)
# unpack parameters
# Assign parameters to expected values
# dnorm(x,mean,sd)

#3 Assess max likelihood of model parameters with optim()
# Create a vector with initial guesses
# Minimization function with appopriate arguments

#T-test:
#pchisq(D,df=1,lower.tail=FALSE) where D = 2*difference of negative log likelihood
#Take output from both likelihood function to calculate p-value

library(ggplot2)

#### Problem 1 ####

#Load and and subset data

ponData=read.csv("ponzr1.csv", stringsAsFactors = FALSE)

subsetI231N=ponData[ponData$mutation%in%c('WT','I231N'),]
subsetM124K=ponData[ponData$mutation%in%c('WT','M124K'),]
subsetV456D=ponData[ponData$mutation%in%c('WT','V456D'),]

#Replace WT entries with 0 and mutant with 1

for(i in subsetI231N[,1]){
  if ('WT' %in% subsetI231N[i,1]){
    subsetI231N[i,1]=0
  } else{
    subsetI231N[i,1]=1
  }
}

for(i in subsetV456D[,1]){
  if ('WT' %in% subsetV456D[i,1]){
    subsetV456D[i,1]=0
  } else{
    subsetV456D[i,1]=1
  }
}

for(i in subsetM124K[,1]){
  if ('WT' %in% subsetM124K[i,1]){
    subsetM124K[i,1]=0
  } else{
    subsetM124K[i,1]=1
  }
}

#Custum negative log likelihood functions (null and with 
#mutation effects)

nlllike=function(p,x,y){
  B0=p[1]
  B1=p[2]
  sigma=exp(p[3])
  
  mRNAcount=B0+B1*x
  
  nll=-sum(dnorm(x=y, mean =mRNAcount, sd=sigma, log=TRUE))
  return(nll)
}

nlllikeNull=function(p,x,y){
  B0=p[1]
  sigma=exp(p[2])
  
  mRNAcount=B0
  
  nll=-sum(dnorm(x=y, mean =mRNAcount, sd=sigma, log=TRUE))
  return(nll)
}
#Calculating parameters; BWT=Avg mRNA count across WT, 
#B0=Avg mRNA Count across all fish, 
#B1 = Diff in WT avg and M124K avg, B2 = " and V456D, B3 = "
#and I213N
#
#Note that when eveluating null model, use B0, but when evaluating
#linear model, use BWT for the "B0" in the function scipt

BWT=mean(ponData$mutation%in% 'WT')
B0=mean(ponData[,2])
B1=(mean(ponData[ponData$mutation%in%'M124K',]))-(mean(ponData[ponData$mutation%in%'WT',]))
B2=(mean(ponData[ponData$mutation%in%'V456D',]))-(mean(ponData[ponData$mutation%in%'WT',]))
B3=(mean(ponData[ponData$mutation%in%'I213N',]))-(mean(ponData[ponData$mutation%in%'WT',]))


#M124K

fitNullM124K=optim(par=c(B0,1), fn=nllikeNull, y=subsetM124K[,2])
fitLinearM124K=optim(par=c(BWT,B1,1), fn=nllike, x=subsetM124K[,1], y=subsetM124K[,2])

pvalM124K=pchisq(q=(fitNullM124K-fitLinearM124K), df=1, lower.tail = FALSE)

#v456D

fitNullv456D=optim(par=c(B0,1), fn=nllikeNull, y=subsetv456D[,2])
fitLinearv456D=optim(par=c(BWT,B2,1), fn=nllike, x=subsetv456D[,1], y=subsetv456D[,2])

pvalv456D=pchisq(q=(fitNullv456D-fitLinearv456D), df=1, lower.tail = FALSE)

#I213N

fitNullI213N=optim(par=c(B0,1), fn=nllikeNull, y=subsetI213N[,2])
fitLinearI213N=optim(par=c(BWT,B2,1), fn=nllike, x=subsetI213N[,1], y=subsetI213N[,2])

pvalI213N=pchisq(q=(fitNullI213N-fitLinearI213N), df=1, lower.tail = FALSE)

#Presenting Results

Res=rbind(c('pvalM124K','pvalv456D','pvalI213N'), c(pvalM124K, pvalv456D, pvalI213N))
Results=data.matrix(Res)






#### Problem 2 ####

#Read in Data
data2 <- read.csv(file = "MmarinumGrowth.csv", header = TRUE)
attach(data2)
df2 <- data.frame(x = S, y = u)

#Graph of raw data
graph2a <- ggplot(df2, aes(x = S, y = u)) + geom_point() + theme_classic()

#Create liklikhood function
###Same as in-class example (continuous, normal distribution)
###BUT, no longer a linear model - change formula for expected
nllike2 <- function(p, x, y){
  umax = p[1]
  Ks = p[2]
  sigma = p[3] #Why is this different from in-class example?
  
  expected = umax*S/(S + Ks)
  
  nll = -sum(dnorm(x = y, mean = expected, sd = sigma, log = TRUE))
  return(nll)
}

#Estimate Umax, Ks, sigma
initialGuess = c(1,1,1)
fit2 = optim(par = initialGuess, fn=nllike2, x = S, y = u)

#Creating data to test the estimated parameters
testx <- seq(0, 1000, length.out = 100)
testy <- fit2$par[1]*testx/(testx+fit2$par[2])
df2b <- data.frame(testx,testy)

#Graphing the estimated model
graph2b <- ggplot(df2b, aes(x = testx, y = testy)) + geom_point() + theme_classic()

#### Problem 3 ####

#Read-in data
data3 <- read.csv(file="leafDecomp.csv", header= TRUE)
attach(data3)


#Three liklihood functions - only difference is expected (y)
nllike3a <- function(p, x, y){
  a = p[1]
  b = p[2]
  c = p[3]
  sigma = p[4] 
  
  expected = a
  
  nll = -sum(dnorm(x = y, mean = expected, sd = sigma, log = TRUE))
  return(nll)
}
nllike3b <- function(p, x, y){
  a = p[1]
  b = p[2]
  c = p[3]
  sigma = p[4] 
  
  expected = a + b*Ms
  
  nll = -sum(dnorm(x = y, mean = expected, sd = sigma, log = TRUE))
  return(nll)
}
nllike3c <- function(p, x, y){
  a = p[1]
  b = p[2]
  c = p[3]
  sigma = p[4] 
  
  expected = a + b*Ms + c*Ms^2
  
  nll = -sum(dnorm(x = y, mean = expected, sd = sigma, log = TRUE))
  return(nll)
}

#Fit the models 
#The order of values are B0, B1, B2, and sigma
initialGuess = c(100,100,100,100)
fit3a = optim(par = initialGuess, fn=nllike3a, x = Ms, y = decomp)
fit3b = optim(par = initialGuess, fn=nllike3b, x = Ms, y = decomp)
initialGuess = c(200,10,-0.2,1)
fit3c = optim(par = initialGuess, fn=nllike3c, x = Ms, y = decomp)

#Conduct goodness of fit tests 
D3a <- 2*fit3a$value
D3b <- 2*fit3b$value
D3c <- 2*fit3c$value

P3_ab <- pchisq(D3a-D3b,df=1, lower.tail = FALSE)
P3_ac <- pchisq(D3a-D3c,df=2, lower.tail = FALSE)
P3_bc <- pchisq(D3b-D3c,df=1, lower.tail = FALSE)
