#Load Libraries
library(stringr)

#Load and and subset data

ponData=read.csv("ponzr1.csv", stringsAsFactors = FALSE)

subsetI213N=ponData[ponData$mutation%in%c('WT','I213N'),]
subsetM124K=ponData[ponData$mutation%in%c('WT','M124K'),]
subsetV456D=ponData[ponData$mutation%in%c('WT','V456D'),]

#Replace WT entries with 0 and mutant with 1

a=str_replace(subsetI213N[,1], "WT", "0")
aNEW=str_replace(a, "I213N", "1")
QI213N=data.frame(aNEW, stringsAsFactors = FALSE)             
I213N=data.matrix(QI213N)

b=str_replace(subsetM124K[,1], "WT", "0")
bNEW=str_replace(b, "M124K", "1")
QM124K=data.frame(bNEW, stringsAsFactors = FALSE)
M124K=data.matrix(QM124K)

c=str_replace(subsetV456D[,1], "WT", "0")
cNEW=str_replace(c, "V456D", "1")
QV456D=data.frame(cNEW, stringsAsFactors = FALSE)
V456D=data.matrix(QV456D)

#Custum negative log likelihood functions (null and with 
#mutation effects)

nllike=function(p,x,y){
  B0=p[1]
  B1=p[2]
  sigma=exp(p[3])
  
  mRNAcount=B0+B1*x
  
  nll=-sum(dnorm(x=y, mean =mRNAcount, sd=sigma, log=TRUE))
  return(nll)
}

nllikeNull=function(p,x,y){
  B0=p[1]
  sigma=exp(p[2])
  
  mRNAcount=B0
  
  nll=-sum(dnorm(x=y, mean =mRNAcount, sd=sigma, log=TRUE))
  return(nll)
}
#Calculating parameters; BWT=Avg mRNA count across WT, 
#B0=Avg mRNA Count across all fish, 
#B1 = Diff in WT avg and M124K avg, B2 = " and V456D, B3 = "
#and I213N
#
#Note that when eveluating null model, use B0, but when evaluating
#linear model, use BWT for the "B0" in the function scipt

BWT=mean(ponData[ponData$mutation%in% 'WT',2])
B0=mean(ponData[,2])
B1=abs((mean(ponData[ponData$mutation %in% 'M124K',2]))-(mean(ponData[ponData$mutation%in%'WT',2])))
B2=abs((mean(ponData[ponData$mutation %in% 'V456D',2]))-(mean(ponData[ponData$mutation%in%'WT',2])))
B3=abs((mean(ponData[ponData$mutation %in% 'I213N',2]))-(mean(ponData[ponData$mutation%in%'WT',2])))


#M124K

fitNullM124K=optim(par=c(B0,2), fn=nllikeNull, y=subsetM124K[,2])
fitLinearM124K=optim(par=c(BWT,B1,3), fn=nllike, x=M124K[1,], y=subsetM124K[,2])

Da=2*fitNullM124K$value
Da1=2*fitLinearM124K$value

pvalM124K=pchisq(q=Da-Da1, df=1, lower.tail = FALSE)

#v456D

fitNullV456D=optim(par=c(B0,2), fn=nllikeNull, y=subsetV456D[,2])
fitLinearV456D=optim(par=c(BWT,B2,2), fn=nllike, x=V456D[1,], y=subsetV456D[,2])

Db=2*fitNullV456D$value
Db1=2*fitLinearV456D$value

pvalV456D=pchisq(q=Db-Db1, df=1, lower.tail = FALSE)

#I213N

fitNullI213N=optim(par=c(B0,2), fn=nllikeNull, y=subsetI213N[,2])
fitLinearI213N=optim(par=c(BWT,B3,2), fn=nllike, x=I213N[1,], y=subsetI213N[,2])

Dc=2*fitNullI213N$value
Dc1=2*fitLinearI213N$value

pvalI213N=pchisq(q=Dc-Dc1, df=1, lower.tail = FALSE)

#Presenting Results

Res=rbind(c('pvalM124K','pvalV456D','pvalI213N'), c(pvalM124K, pvalv456D, pvalI213N))
Results=data.matrix(Res)

