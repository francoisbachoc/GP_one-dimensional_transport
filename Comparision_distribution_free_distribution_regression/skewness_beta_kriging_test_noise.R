rm( list=ls() )   
source("../functions.R")

#######################################################
#settings 
#######################################################
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
nugget = 10^(-8)
theta_inf = c(0.1,0.01,0.25,10^(-8))
theta_sup = c(100,10,1,0.1)
init_theta = c(1,1,1,10^(-4))
ind_active=c(1,2,3,4)
cov_matrix_function = function(theta,mq) { cov_matrix_exp_quantile(theta[1],theta[2],theta[3],mq,mq) + theta[4]*diag(dim(mq)[1]) }
cross_cov_matrix_function = function(theta,mq1,mq2) { cov_matrix_exp_quantile(theta[1],theta[2],theta[3],mq1,mq2) }
n_gradient_descent=10
name="skewness_beta_kriging_test_noise"

#######################################################
#execution
#######################################################
va = seq(from = amin,to=amax,length=100)
ya = skewness(va,b)

mtraining = matrix(nrow=ntraining,ncol=nsample,data=-1)
ytraining = rep(0,ntraining)
atraining = runif(n=ntraining,min=amin,max=amax)
for (i in 1:ntraining) {
  a = atraining[i]
  mtraining[i,] = rbeta(n=nsample,shape1=a,shape2=b)
  ytraining[i] = skewness(a,b) + rnorm(1)
}


mvalidation = matrix(nrow=nvalidation,ncol=nsample,data=-1)
yvalidation = rep(0,nvalidation)
avalidation = runif(n=nvalidation,min=amin,max=amax)
for (i in 1:nvalidation) {
  a = avalidation[i]
  mvalidation[i,] = rbeta(n=nsample,shape1=a,shape2=b)
  yvalidation[i] = skewness(a,b) + rnorm(1)
}

mtest = matrix(nrow=ntest,ncol=nsample,data=-1)
ytest = rep(0,ntest)
atest = runif(n=ntest,min=amin,max=amax)
for (i in 1:ntest) {
  a = atest[i]
  mtest[i,] = rbeta(n=nsample,shape1=a,shape2=b)
  ytest[i] = skewness(a,b) + rnorm(1)
}


#maximum likelihood
mq_obs = rbind(mtraining,mvalidation) 
for (i in 1:dim(mq_obs)[1]) {
  mq_obs[i,] = sort(mq_obs[i,])
}

y_obs=c(ytraining,yvalidation)

res_ML = maximum_likelihood(init_theta=init_theta,theta_inf=theta_inf,theta_sup=theta_sup,
                               ind_active=ind_active,nugget=nugget,cov_matrix_function=cov_matrix_function,
                               mq_obs=mq_obs,y_obs=y_obs,k=n_gradient_descent)  
hat_theta = res_ML$hat_theta
nugget = max(nugget,hat_theta[4])

mq_pred = mtest 
for (i in 1:dim(mq_pred)[1]) {
  mq_pred[i,] = sort(mq_pred[i,])
}

R_obs = cov_matrix_function(hat_theta,mq_obs)
R_pred = cov_matrix_function(hat_theta,mq_pred)
R_obs_pred = cross_cov_matrix_function(hat_theta,mq_obs,mq_pred)
predtest = prediction(R_obs,R_obs_pred,y_obs,nugget)
hat_sigma_test = predictive_sd(R_obs,R_pred,R_obs_pred,nugget)



#######################################################
#Save of the results
#######################################################
save(file=paste0(name,".Rdata"),atest,predtest,va,ya,mtest,mvalidation,
     yvalidation,avalidation,mtraining,atraining,ytraining,res_ML,hat_theta,hat_sigma_test)








