rm( list=ls() )   
source("../functions.R")

#############################################
#Settings
#############################################
set.seed(1)
skewness = function( a,b ) {
  (2*(b-a)*sqrt(a+b+1)) / ((a+b+2)*sqrt(a*b))
}

amin=3
amax=20
b=3

ntraining = 250
nvalidation = 25
ntest= 50
nsample = 500
nDistEval=4096
nbh=100
hlim=c(0,1)
blim=c(0,1)
name = "skewness_beta_distribution_regression"

#############################################
#Execution
#############################################
va = seq(from = amin,to=amax,length=100)
ya = skewness(va,b)

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



bh = select_bh_kernel_regression(mtraining = mtraining,ytraining = ytraining,
                                 mvalidation = mvalidation,yvalidation = yvalidation,
                                 nDistEval=nDistEval,nbh=nbh,hlim=hlim,blim=blim)

predtest = kernel_regression(rbind(mtraining,mvalidation),c(ytraining,yvalidation),
                             mtest,h=bh[2],b=bh[1],nDistEval=nDistEval)


##################################
#Save of the results
##################################

save(file=paste0(name,".Rdata"),bh,atest,predtest,va,ya,mtest,mvalidation,
     yvalidation,avalidation,mtraining,atraining,ytraining)

