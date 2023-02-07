rm( list=ls() )   
source("../functions.R")

skewness = function( a,b ) {
  (2*(b-a)*sqrt(a+b+1)) / ((a+b+2)*sqrt(a*b))
}

amin=3
amax=20
b=3
va=seq(from=amin,to=amax,length=100)
vb = seq(from=b,to=b,length=100)


ntraining = 250
nvalidation = 25
ntest= 50
nsample = 500
nDistEval=100
nbh=10
hlim=c(0.1,0.3)
blim=c(0.1,0.3)

mtraining = matrix(nrow=ntraining,ncol=nsample,data=-1)
ytraining = rep(0,ntraining)
atraining = runif(n=ntraining,min=amin,max=amax)
for (i in 1:ntraining) {
  a = atraining[i]
  mtraining[i,] = rbeta(n=nsample,shape1=a,shape2=b)
  ytraining[i] = skewness(a,b)
}


mvalidation = matrix(nrow=nvalidation,ncol=nsample,data=-1)
yvalidation = rep(0,nvalidation)
avalidation = runif(n=nvalidation,min=amin,max=amax)
for (i in 1:nvalidation) {
  a = avalidation[i]
  mvalidation[i,] = rbeta(n=nsample,shape1=a,shape2=b)
  yvalidation[i] = skewness(a,b)
}

mtest = matrix(nrow=ntest,ncol=nsample,data=-1)
ytest = rep(0,ntest)
atest = runif(n=ntest,min=amin,max=amax)
for (i in 1:ntest) {
  a = atest[i]
  mtest[i,] = rbeta(n=nsample,shape1=a,shape2=b)
  ytest[i] = skewness(a,b)
}


predtest = kernel_regression(rbind(mtraining,mvalidation),c(ytraining,yvalidation),
                             mtest,h=0.1,b=0.1,nDistEval=nDistEval)

plot(va,skewness(va,vb),type="l",xlab="",ylab="")
points(atest,predtest,type="p",pch=8,xlab="",ylab="")

bh = select_bh_kernel_regression(mtraining = mtraining,ytraining = ytraining,
                                 mvalidation = mvalidation,yvalidation = yvalidation,
                                 nDistEval=nDistEval,nbh=nbh,hlim=hlim,blim=blim)

predtest = kernel_regression(rbind(mtraining,mvalidation),c(ytraining,yvalidation),
                             mtest,h=bh[2],b=bh[1],nDistEval=nDistEval)

plot(va,skewness(va,vb),type="l",xlab="",ylab="")
points(atest,predtest,type="p",pch=8,xlab="",ylab="")