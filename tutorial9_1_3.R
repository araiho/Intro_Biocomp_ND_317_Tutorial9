nllike<-function(p,x,w){
  B0=p[1]
  expected=B0
  nll=-sum(dnorm(x=w, mean = expected, log = TRUE))
  return((nll))}
nllike<-function(p,x,w){
  B0=p[1]
  B1=p[2]
  sigma=exp(p[3])
  expected=B0+B1*x
  nll=-sum(dnorm(x=w, mean = expected, sd = sigma, log = TRUE))
  return((nll))}
initialGuess=c(1,1,1)
fit=optim(par=initialGuess,fn=nllike,x=df$x,w=df$x)
nll6=2588.043
nll7= -330.6372
nlld3= nll6-nll7
D3=2*nlld3
pchisq(D3, df=1, lower.tail = FALSE)
=1