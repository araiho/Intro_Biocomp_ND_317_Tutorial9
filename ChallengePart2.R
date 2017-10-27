#Pseudo-Code

#1) Load Data
#
#2) custom function(params, observation)
# - unpack parameters
# - assign parameters to expected values
# - dnorm(x,mean,sd)
#
#3) Asses max likelihood of model parameters with optim()
# - create vector with initial guess
# - minimization function with appropriate arguments
#
#pchisq(D,df=1)
#D = 2times difference of negative log-likelihood

#Loading Data
data = read.csv(file = "MmarinumGrowth.csv", header = TRUE, sep = ",")

s=data[,1]
u=data[,2]

nllike<-function(p,x,y){
  Umax=p[1]
  Ks=p[2]
  sigma=p[3]
  
  expected=Umax*s/(s+Ks)
  
  nll=-sum(dnorm(x=y,mean=expected,sd=sigma,log=TRUE))
  return(nll)
}

initialGuess=c(1,1,1)
fit=optim(par=initialGuess,fn=nllike,x=s,y=u)

print(fit)
