setwd("C:/Users/DAVIS/Desktop/shell-novice-data/exe9/Intro_Biocomp_ND_317_Tutorial9/")

#Question1 V456D is sig
Ponzr <- read.csv("ponzr1.csv", header=T)

like1 <- function(p,x,y){
  B0=p[1]
  B1=p[2]
  sigma = exp(p[3])
  expected = B0
  
  nll = -sum(dnorm(y, mean=expected, sd=sigma, log=TRUE))
  return(nll)
}

like2 <- function(p,x,y){
  B0=p[1]
  B1=p[2]
  sigma = exp(p[3])
  expected = B0+B1*x
  
  nll = -sum(dnorm(y, mean=expected, sd=sigma, log=TRUE))
  return(nll)
}


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

Guess = c(1,1,1)
fit1=optim(Guess, like1, x=Ponzr1$mutation, y=Ponzr1$ponzr1Counts)
fit2=optim(Guess, like2, x=Ponzr2$mutation, y=Ponzr2$ponzr1Counts)
fit1$value
fit2$value
D = 2*(fit1$value-fit2$value)
pchisq(D, df=1, lower.tail=F)


Guess = c(1,1,1)
fit1=optim(Guess, like1, x=Ponzr1$mutation, y=Ponzr1$ponzr1Counts)
fit2=optim(Guess, like2, x=Ponzr3$mutation, y=Ponzr3$ponzr1Counts)
fit1$value
fit2$value
D = 2*(fit1$value-fit2$value)
pchisq(D, df=1, lower.tail=F)

Guess = c(1,1,1)
fit1=optim(Guess, like1, x=Ponzr1$mutation, y=Ponzr1$ponzr1Counts)
fit2=optim(Guess, like2, x=Ponzr4$mutation, y=Ponzr4$ponzr1Counts)
fit1$value
fit2$value
D = 2*(fit1$value-fit2$value)
pchisq(D, df=1, lower.tail=F)

#Question2
Mmarinum <- read.csv("MmarinumGrowth.csv", header=T)
nllike <- function(p,x,y){
  B0=p[1]
  B1=p[2]
  sigma = exp(p[3])
  expected = B0*(x/(B1+x))
  
  nll = -sum(dnorm(x=y, mean=expected, sd=sigma, log=TRUE))
  return(nll)
}
Guess = c(1,1,1)
fit=optim(Guess, nllike, x=Mmarinum$S, y=Mmarinum$u)
print(cbind("umax is ", fit$par[1], " and Ks is ", fit$par[2]))

#Question3 (not finished yet)

leafDecomp <- read.csv("leafDecomp.csv", header = T)

one_var <- function(p,x,y){
  B0=p[1]
  B1=p[2]
  B2=p[3]
  sigma = exp(p[4])
  expected = B0
  
  nll = -sum(dnorm(x=y, mean=expected, sd=sigma, log=TRUE))
  return(nll)
}

two_var <- function(p,x,y){
  B0=p[1]
  B1=p[2]
  B2=p[3]
  sigma = exp(p[4])
  expected = B0 + B1*x
  
  nll = -sum(dnorm(x=y, mean=expected, sd=sigma, log=TRUE))
  return(nll)
}

three_var <- function(p,x,y){
  B0=p[1]
  B1=p[2]
  B2=p[3]
  sigma = exp(p[4])
  expected = B0 + B1*x + B2*x*x
  
  nll = -sum(dnorm(x=y, mean=expected, sd=sigma, log=TRUE))
  return(nll)
}

Guess = c(1,1,1,1)
one_var_result=optim(Guess, one_var, x=leafDecomp$Ms, y=leafDecomp$decomp)
two_var_result=optim(Guess, two_var, x=leafDecomp$Ms, y=leafDecomp$decomp)
three_var_result=optim(Guess, three_var, x=leafDecomp$Ms, y=leafDecomp$decomp)

one_var_result$value
two_var_result$value
three_var_result$value

one_var_result$par
two_var_result$par
three_var_result$par

D1_2 = 2*(one_var_result$value-two_var_result$value)
D1_3 = 2*(one_var_result$value-three_var_result$value)
D2_3 = 2*(three_var_result$value-two_var_result$value)

pchisq(D1_2, df=1, lower.tail=F)
pchisq(D1_3, df=2, lower.tail=F)
pchisq(D2_3, df=1, lower.tail=F)





