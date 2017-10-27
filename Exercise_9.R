#Question 1

data <- read.csv("/Users/elizabethfortin12/Documents/ND First Year/Biocomputing/R_Programming/Exercise_9/ponzr1.csv")

data$WTorNOT <- 0

for (i in 1:nrow(data)){
  if (data[i,1] == "WT"){
    data$WTorNOT[i] <- 0
  }
  else{
    data$WTorNOT[i] <- 1
  }
}

#Null model
nll1 <- function(p,y){
  B0=p[1]
  sig=exp(p[2])
  expected=B0
  nll=-sum(dnorm(x=y, mean=expected, sd=sig,log=TRUE))
  return(nll)
}

initialGuess <- c(2000,1)
fit <- optim(par=initialGuess,fn=nll1,y=data$ponzr1Counts)
fit

#Linear Model
nll2 <- function(p,x,y){
  B0=p[1]
  B1=p[2]
  sig=exp(p[3])
  expected=B0+B1*x
  nll=-sum(dnorm(x=y, mean=expected, sd=sig,log=TRUE))
  return(nll)
}

initialGuess <- c(1,1,1)
fit2 <- optim(par=initialGuess,fn=nll2,x=data$WTorNOT,y=data$ponzr1Counts)
fit2

#likelihood ratio test
D<-2*(fit2$value - fit$value)
ratiotest <- pchisq(q=D, df=1, lower.tail=FALSE)
ratiotest




#Question 3

decomp <- read.csv("/Users/elizabethfortin12/Documents/ND First Year/Biocomputing/R_Programming/Exercise_9/leafDecomp.csv")

#Null model
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

#Humpshaped model
nllhump <- function(p,Ms,d){
  a=p[1]
  b=p[2]
  c=p[3]
  sig=exp(p[4])
  expected=a+b*Ms+c*Ms^2
  nll=-sum(dnorm(x=d, mean=expected, sd=sig,log=TRUE))
  return(nll)
}

initialGuess <- c(1,1,1,1)
hump <- optim(par=initialGuess,fn=nllhump,Ms=decomp$Ms,d=decomp$decomp)
hump

#likelihood ratio test
E<-2*(hump$value-null$value)
ratiotest2 <- pchisq(q=E, df=2, lower.tail=FALSE)
ratiotest2
