###Exercise9_Question2###

#Load data
setwd("~/Desktop/data-shell/Exercise9/Intro_Biocomp_ND_317_Tutorial9/")
M <- read.csv("MmarinumGrowth.csv", header=TRUE)

#Custom Function
nllike <- function(p,x,y){
  umax=p[1]
  S=p[2]
  Ka=p[3]
  sigma=exp(p[4])
  
  expected=umax*(S/(S+Ka))
    
  
  nll <- -sum(dnorm(0,mean=expected,sd=sigma,log=TRUE))
  return(nll)
}

###expected=umax(S/(S+Ka))   p[1]*(p[2]/(p[2]+p[3]))

#Estimate parameters by minimizing the negative log likihood.
intialGuess <- c(0.5,0.5,40,0.001)
fit <- optim(par=intialGuess,fn=nllike,x=M$x,y=M$y)

print(fit)

