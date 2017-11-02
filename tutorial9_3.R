<<<<<<< HEAD
data=read.csv(file = 'leafDecomp.csv', header = TRUE, sep = ",")
Ms <- data[,1]
d <- data[,2]
nllike<-function(p,x,y){
  d=p[1]
  a=p[2]
  expected=(d=a)
  nll=-sum(dnorm(x=y, mean = expected, log = TRUE))
  return((nll))}
nllike<-function(p,x,y){
  d=p[1]
  a=p[2]
  b=p[3]
  sigma=exp(p[4])
  expected=(d=a*Ms+b*Ms^2)
  nll=-sum(dnorm(x=y, mean = expected, sd = sigma, log = TRUE))
  return((nll))}
initialGuess=c(1,1,1,1)
fit=optim(par=initialGuess,fn=nllike,x=d,y=Ms)
print(fit)
value1=10723.22
value2=96.4776
v=value1-value2
=======
data=read.csv(file = 'leafDecomp.csv', header = TRUE, sep = ",")
Ms <- data[,1]
d <- data[,2]
nllike<-function(p,x,y){
  d=p[1]
  a=p[2]
  expected=(d=a)
  nll=-sum(dnorm(x=y, mean = expected, log = TRUE))
  return((nll))}
nllike<-function(p,x,y){
  d=p[1]
  a=p[2]
  b=p[3]
  sigma=exp(p[4])
  expected=(d=a*Ms+b*Ms^2)
  nll=-sum(dnorm(x=y, mean = expected, sd = sigma, log = TRUE))
  return((nll))}
initialGuess=c(1,1,1,1)
fit=optim(par=initialGuess,fn=nllike,x=d,y=Ms)
print(fit)
value1=10723.22
value2=96.4776
v=value1-value2
p=0
>>>>>>> 9cf11ca599d407d6eb4bd775a40442e4c1891c8c
p=0