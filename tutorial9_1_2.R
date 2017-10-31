nllike<-function(p,x,z){
  B0=p[1]
  expected=B0
  nll=-sum(dnorm(x=z, mean = expected, sd = sigma, log = TRUE))
  return((nll))}nllike<-function(p,x,z){
  B0=p[1]
  B1=p[2]
  sigma=exp(p[3])
  expected=B0+B1*x
  nll=-sum(dnorm(x=z, mean = expected, sd = sigma, log = TRUE))
  return((nll))}
initialGuess=c(1,1,1)
fit=optim(par=initialGuess,fn=nllike,x=df$x,z=df$x)
nll3=102.0236
nll4=-330.6372
nlld2=nll3-nll4=432.6608
p=3.393539e-190