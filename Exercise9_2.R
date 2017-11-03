#### Exercise 9 Question 2 ####
## Upload Data ##
Data <- read.csv("~/Desktop/data-shell/Intro_Biocomp_ND_317_Tutorial9/MmarinumGrowth.csv",sep = ',',header=TRUE)
View(Data)
## Define a function ##
Monod <- function(p,x,y){
  S = p[1]
  K = p[2]
  sigma = exp(p[3])
  expected = (S*(x/(x+K)))
  mnd= -sum(dnorm(x=y, mean=expected, sd = sigma, log = TRUE))
  return(mnd)
}
## Define parameters ##
params2 <- c(2,2,0)
fit4 <- optim(par = params2, fn=Monod, x = Data$S, y=Data$u)
fit4