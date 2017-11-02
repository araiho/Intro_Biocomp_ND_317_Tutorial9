#Load and and subset data

ponData=read.csv(ponzr1.csv)

subsetI231N=data[data$mutation%in%c('WT','I231N'),]
subsetM124K=data[data$mutation%in%c('WT','M124K'),]
subsetV456D=data[data$mutation%in%c('WT','V456D'),]

#Replace WT entries with 0 and mutant with 1

for(i in subsetI231N[,1]){
  if (subsetI231N[i,1]=='WT'){
    subsetI231N[i,1]=0
  } else{
    subsetI231N[i,1]=1
  }
}

for(i in subsetV456D[,1]){
  if (subsetV456D[i,1]=='WT'){
    subsetV456D[i,1]=0
  } else{
    subsetV456D[i,1]=1
  }
}

for(i in subsetM124K[,1]){
  if (subsetM124K[i,1]=='WT'){
    subsetM124K[i,1]=0
  } else{
    subsetM124K[i,1]=1
  }
}

#Custum negative log likelihood functions (null and with 
#mutation effects)

nlllike=function(p,x,y){
  B0=p[1]
  B1=p[2]
  sigma=exp(p[3])
  
  mRNAcount=B0+B1*x
  
  nll=-sum(dnorm(x=y, mean =mRNAcount, sd=sigma, log=TRUE))
  return(nll)
}

nlllikeNull=function(p,x,y){
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

BWT=mean[ponData$mutation%in%c('WT'),]
B0=mean(ponData[,2])
B1=(mean[ponData$mutation%in%c('M124K'),])-(mean[ponData$mutation%in%c('WT'),])
B2=(mean[ponData$mutation%in%c('V456D'),])-(mean[ponData$mutation%in%c('WT'),])
B3=(mean[ponData$mutation%in%c('I213N'),])-(mean[ponData$mutation%in%c('WT'),])

#M124K

fitNullM124K=optim(par=c(B0,1), fn=nllikeNull, y=subsetM124K[,2])
fitLinearM124K=optim(par=c(BWT,B1,1), fn=nllike, x=subsetM124K[,1], y=subsetM124K[,2])

pvalM124K=pchisq(q=(fitNullM124K-fitLinearM124K), df=1, lower.tail = FALSE)

#v456D

fitNullv456D=optim(par=c(B0,1), fn=nllikeNull, y=subsetv456D[,2])
fitLinearv456D=optim(par=c(BWT,B2,1), fn=nllike, x=subsetv456D[,1], y=subsetv456D[,2])

pvalv456D=pchisq(q=(fitNullv456D-fitLinearv456D), df=1, lower.tail = FALSE)

#I213N

fitNullI213N=optim(par=c(B0,1), fn=nllikeNull, y=subsetI213N[,2])
fitLinearI213N=optim(par=c(BWT,B2,1), fn=nllike, x=subsetI213N[,1], y=subsetI213N[,2])

pvalI213N=pchisq(q=(fitNullI213N-fitLinearI213N), df=1, lower.tail = FALSE)

#Presenting Results

Res=rbind(c('pvalM124K','pvalv456D','pvalI213N'), c(pvalM124K, pvalv456D, pvalI213N))
Results=data.matrix(Res)

