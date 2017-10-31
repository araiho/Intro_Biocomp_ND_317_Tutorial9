#Question 1

library(ggplot2)

ggplot(data)+
  geom_point(aes(x=mutation, y=ponzr1Counts))

ggplot(data)+
  geom_boxplot(aes(x=mutation, y=ponzr1Counts))

data <- read.csv("/Users/elizabethfortin12/Documents/ND First Year/Biocomputing/R_Programming/Exercise_9/ponzr1.csv")

# adding a column for mutated = 1, wt = 0
for (i in 1:nrow(data)){
  if (data$mutation[i] == 'WT'){
    data$mutorno[i] <- 0
  }
  else{
    data$mutorno[i] <- 1
  }
}

# data frames containing only one mutation
m124k <- data[data$mutation=='WT',]
m124k[11:20,] <- data[data$mutation=='M124K',]
v456d <- data[data$mutation=='WT',]
v456d[11:20,] <- data[data$mutation=='V456D',]
i213n <- data[data$mutation=='WT',]
i213n[11:20,] <- data[data$mutation=='I213N',]

results <- matrix(0,3,3)
colnames(results) <- c("null", "linear", "chisq")

#Null model
nllnull <- function(p,y){
  B0=p[1]
  sig=exp(p[2])
  expected=B0
  nll=-sum(dnorm(x=y, mean=expected, sd=sig,log=TRUE))
  return(nll)
}

#Linear Model
nlllinear <- function(p,x,y){
  B0=p[1]
  B1=p[2]
  sig=exp(p[3])
  expected=B0+B1*x
  nll=-sum(dnorm(x=y, mean=expected, sd=sig,log=TRUE))
  return(nll)
}

#Null model M124K
initialGuess <- c(2000,10)
fit <- optim(par=initialGuess,fn=nllnull,y=m124k$ponzr1Counts)
results[1,1] <- fit$value

#Linear model M124K
initialGuess <- c(2000,-1000,10)
fit <- optim(par=initialGuess,fn=nlllinear,x=m124k$mutorno,y=m124k$ponzr1Counts)
results[1,2] <- fit$value

#Null model V456D
initialGuess <- c(2000,10)
fit <- optim(par=initialGuess,fn=nllnull,y=v456d$ponzr1Counts)
results[2,1] <- fit$value

#Linear Model V456D
initialGuess <- c(2000,-1000,10)
fit <- optim(par=initialGuess,fn=nlllinear,x=v456d$mutorno,y=v456d$ponzr1Counts)
results[2,2] <- fit$value

#Null model i213n
initialGuess <- c(2000,10)
fit <- optim(par=initialGuess,fn=nllnull,y=i213n$ponzr1Counts)
results[3,1] <- fit$value

#Linear Model i213n
initialGuess <- c(2000,-1000,10)
fit <- optim(par=initialGuess,fn=nlllinear,x=i213n$mutorno,y=i213n$ponzr1Counts)
results[3,2] <- fit$value

#likelihood ratio test
A <- 2*(results[1,1]-results[1,2])
B <- 2*(results[2,1]-results[2,2])
C <- 2*(results[3,1]-results[3,2])
results[1,3] <- pchisq(q=A, df=1, lower.tail=FALSE) 
results[2,3] <- pchisq(q=B, df=1, lower.tail=FALSE) 
results[3,3] <- pchisq(q=C, df=1, lower.tail=FALSE) 

#__________________________________________________________

# Quesiton 2

# Load Data

data1=read.csv("/Users/chelseaweibel/Desktop/Biocomputing/R.data/Biocomputing Exercise 9.Use/MmarinumGrowth.csv")

# Custom Function

nllike = function(p,S,u){
  umax=p[1]
  Ks=p[2]
  sigma=exp(p[3])
  
  expected=umax*(S/(S+Ks))
  
  nll=-sum(dnorm(x=u,mean=expected,sd=sigma,log=TRUE))
  return(nll)
}

# Initial Guess

initialGuess=c(1,1,1)

# Assess max likelihood of model parameters

fit3=optim(par=initialGuess,fn=nllike,S=data1$S,u=data1$u)
print(fit3)

# P1 = max growth rate (umax) = 1.46
# P2 = half-saturation constant = 42.6
# P3 = -sigma squared = -3.13
# sigma = 0.04

#__________________________________________________________

#Question 3

decomp <- read.csv("/Users/elizabethfortin12/Documents/ND First Year/Biocomputing/R_Programming/Exercise_9/leafDecomp.csv")

############## Null model
nllnull <- function(p,d){
  a=p[1]
  sig=exp(p[2])
  expected=a
  nll=-sum(dnorm(x=d, mean=expected, sd=sig,log=TRUE))
  return(nll)
}

initialGuess <- c(500,1)
null <- optim(par=initialGuess,fn=nllnull,d=decomp$decomp)
null

#### Results
# a = 590
# sigma = 164

################ Linear Model

nllike_linear=function(p,Ms,d){
  a=p[1]
  b=p[2]
  sigma=exp(p[3])
  
  expected=a+b*Ms
  
  nll_linear=-sum(dnorm(x=d,mean=expected,sd=sigma,log=TRUE))
  return(nll_linear)
}

initialGuess=c(500,1,1)
fit_linear=optim(par=initialGuess,fn=nllike_linear,Ms=decomp$Ms,d=decomp$decomp)

print(fit_linear)

#### Results

# a = 317
# b = 6.3
# sigma = 54

################### Humpshaped model

nllhump <- function(p,Ms,d){
  a=p[1]
  b=p[2]
  c=p[3]
  sig=exp(p[4])
  expected=a+b*Ms+c*Ms^2
  nll=-sum(dnorm(x=d, mean=expected, sd=sig,log=TRUE))
  return(nll)
}

initialGuess <- c(200,10,1,1)
hump <- optim(par=initialGuess,fn=nllhump,Ms=decomp$Ms,d=decomp$decomp)
hump

### Results

# a = 207
# b = 15
# c = -0.11
# sigma = 20

##### The linear model is better than the null model, but the hump 
# shaped model is the best model for this data set.

############# likelihood ratio test
F<-2*(fit_linear$value-null$value)
ratiotest3 <- pchisq(q=F, df=1, lower.tail=FALSE)
ratiotest3

G<-2*(hump$value-fit_linear$value)
ratiotest4 <- pchisq(q=G, df=1, lower.tail=FALSE)
ratiotest4

E<-2*(hump$value-null$value)
ratiotest2 <- pchisq(q=E, df=2, lower.tail=FALSE)
ratiotest2
