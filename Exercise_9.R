# Quesiton 2

# Load Data

data=read.csv("/Users/chelseaweibel/Desktop/Biocomputing/R.data/Biocomputing Exercise 9.Use/MmarinumGrowth.csv")

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

fit=optim(par=initialGuess,fn=nllike,S=data$S,u=data$u)

print(fit)
