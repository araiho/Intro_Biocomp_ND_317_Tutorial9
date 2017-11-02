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
