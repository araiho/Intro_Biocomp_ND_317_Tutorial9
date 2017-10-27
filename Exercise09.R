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
#pchisq(D,df=1) where D = 2*difference of negative log likelihood
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
  sigma = exp(p[3])
  
  expected = umax*S/(S + Ks)
  
  nll = -sum(dnorm(x = y, mean = expected, sd = sigma, log = TRUE))
  return(nll)
}

#Estimate Umax, Ks, sigma
initialGuess = c(1,1,1)
fit2 = optim(par = initialGuess, fn=nllike2, x = S, y = u)

#Creating data to test the estimated parameters
testx <- seq(0, 1000, length.out = 100)
testy <- fit$par[1]*testx/(testx+fit$par[2])
df2b <- data.frame(testx,testy)

#Graphing the estimated model
graph2b <- ggplot(df2b, aes(x = testx, y = testy)) + geom_point() + theme_classic()

D2 <- 2*fit2$value
pchisq(D2,df=1)