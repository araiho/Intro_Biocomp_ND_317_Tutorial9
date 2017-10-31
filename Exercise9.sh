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

