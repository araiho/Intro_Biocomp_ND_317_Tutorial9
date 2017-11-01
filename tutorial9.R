# load data
data = read.csv(file = "ponzr1.csv", header = TRUE, sep = ",")
# custom function(params,observations)
# unpack parameters
#assign variables
WT <- data[ which(data$mutation=='WT'), 2]
M124K <- data[ which(data$mutation==' M124K'), 2]
V456D<- data[ which(data$mutation==' V456D'), 2]
I213N <- data[ which(data$mutation=='I213N'), 2]
# assign parameters to expected values
nllike<-function(p,x,y){
  B0=p[1]
  E=p[2]
  sigma=exp(p[3])
  expected=B0+E
  nll=-sum(dnorm(x=y,mean=expected,sd=sigma,log=TRUE))
  return(nll)
}
# dnorm(x,mean,sd)
# create vector with inital guesses
# assess max likelihood of model parameters with optim()
initialGuess=c(1,1,1)
fit=optim(par=initialGuess,fn=nllike,x=WT,y=V456D)
# value=0

  nllike<-function(p,x,y){
    B0=p[1]
    B1=p[2]
    sigma=exp(p[3])
     E=p[4]
    expected=B0+B1*x+E
    nll=-sum(dnorm(x=y,mean=expected,sd=sigma,log=TRUE))
    return(nll)
  }
  initialGuess=c(1,1,1)
  fit=optim(par=initialGuess,fn=nllike,x=WT,y=V456D)
  value=0
  
  # minimize function with appropriate arguments
  # pchisq(D,df=1)
  # D = 2 times the difference of the negative log likelihoods
  p=1
  
  nllike<-function(p,x,y){
    B0=p[1]
    E=p[2]
    sigma=exp(p[3])
    expected=B0+E
    nll=-sum(dnorm(x=y,mean=expected,sd=sigma,log=TRUE))
    return(nll)
  }
    
    initialGuess=c(1,1,1)
    fit=optim(par=initialGuess,fn=nllike,x=WT,y=M124K)
    
    nllike<-function(p,x,y){
      B0=p[1]
      B1=p[2]
      sigma=exp(p[3])
      E=p[4]
      expected=B0+B1*x+E
      nll=-sum(dnorm(x=y,mean=expected,sd=sigma,log=TRUE))
      return(nll)
    }
    initialGuess=c(1,1,1,1)
    fit=optim(par=initialGuess,fn=nllike,x=WT,y=M124K)
    
    nllike<-function(p,x,y){
      B0=p[1]
      E=p[2]
      sigma=exp(p[3])
      expected=B0+E
      nll=-sum(dnorm(x=w, mean = expected, log = TRUE))
      return((nll))}
    nll1=258655.1
    
    nllike<-function(p,x,y){
      B0=p[1]
      B1=p[2]
      sigma=exp(p[3])
      E=p[4]
      expected=B0+B1*x+E
      nll=-sum(dnorm(x=y, mean = expected, sd = sigma, log = TRUE))
      return((nll))}
    
    initialGuess=c(1,1,1,1)
    fit=optim(par=initialGuess,fn=nllike,x=WT,y=I213N)
    nll2=73.96622
    nlldf=nll1-nll2
    D=517162.3
p=0