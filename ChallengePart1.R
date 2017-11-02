rm(list=ls())

library(stringr)

data = read.csv(file = "ponzr1.csv", header = TRUE, sep = ",")


Nullnllike<-function(p,y){
  B0=p[1]
  sigma=exp(p[2])
  
  expected=B0
  
  nll=-sum(dnorm(x=y, mean=expected, sd=sigma, log = TRUE))
  return((nll))}

nllike<-function(p,x,y){
  B0=
  B1=p[2]
  sigma=exp(p[3])
  
  expected=B0+B1*x
  
  nll=-sum(dnorm(x=y, mean = expected, sd = sigma, log = TRUE))
  return((nll))}

data1=data[which(data$mutation=="WT"),]
data1[,1]=0
data2=data[which(data$mutation=="M124K"),]
data2[,1]=1
data3=data[which(data$mutation=="V456D"),]
data3[,1]=1
data4=data[which(data$mutation=="I213N"),]
data4[,1]=1

subset1 = rbind(data1,data2)
subset2 = rbind(data1,data3)
subset3 = rbind(data1,data4)

NullGuess=c(2000,10)
Guess=c(2000,-1000,10)

#Testing M124K
Fit=optim(Guess, nllike, x=subset1$mutation, y=subset1$ponzr1Counts)
NullFit=optim(NullGuess, Nullnllike, y=subset1$ponzr1Counts)
A = 0.007*(Fit$value-NullFit$value)
M124Kp=pchisq(q=A, df=1, lower.tail=FALSE)

#Testing V456D
Fit=optim(Guess, nllike, x=subset2$mutation, y=subset2$ponzr1Counts)
NullFit=optim(NullGuess, Nullnllike, y=subset2$ponzr1Counts)
B = 2*(Fit$value-NullFit$value)
V456Dp=pchisq(q=B, df=1, lower.tail=F)

#Testing I213N
Fit=optim(Guess, nllike, x=subset3$mutation, y=subset3$ponzr1Counts)
NullFit=optim(NullGuess, Nullnllike, y=subset3$ponzr1Counts)
C = 0.001*(Fit$value-NullFit$value)
I213Np=pchisq(q=C, df=1, lower.tail=F)

print("The M124K p-value is")
print(M124Kp)
print("The V456D p-value is")
print(V456Dp)
print("The I213N p-value is")
print(I213Np)
