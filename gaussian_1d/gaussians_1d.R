rm( list=ls() )   
set.seed(1)  
source("../functions.R")

#######################################################
#settings
#######################################################
sigma_obs=c(0.5,0.75,1,1.5,2)
sigma_pred = seq(from=0.5,to=2,length=100)
sigma=10
lc=0.2
H=1
d = 100
nugget=10^(-8)

#######################################################
#execution
#######################################################
n = length(sigma_obs)
n_pred = length(sigma_pred) 

M_obs = matrix(nrow=n,ncol=d,data=NA)
for (i in 1:n) {
  M_obs[i,] = qnorm(p=seq(from=1/d,to=1-1/d,length=d),mean=0,sd=sigma_obs[i])
}

M_pred = matrix(nrow=n_pred,ncol=d,data=NA)
for (i in 1:n_pred) {
  M_pred[i,] = qnorm(p=seq(from=1/d,to=1-1/d,length=d),mean=0,sd=sigma_pred[i])
}

y_obs = seq(from=0,to=0,length=n)
for (i in 1:n) {
  y_obs[i] = 3 * sigma_obs[i]^4 - 5*sigma_obs[i]^2
}

y_pred = seq(from=0,to=0,length=n)
for (i in 1:n_pred) {
  y_pred[i] = 3 * sigma_pred[i]^4 - 5*sigma_pred[i]^2
}

R_obs = cov_matrix_exp_quantile(sigma,lc,H,M_obs,M_obs)
R_pred = cov_matrix_exp_quantile(sigma,lc,H,M_pred,M_pred)
R_obs_pred = cov_matrix_exp_quantile(sigma,lc,H,M_obs,M_pred)
hat_y_pred = prediction(R_obs,R_obs_pred,y_obs,nugget)
hat_sigma_pred = predictive_sd(R_obs,R_pred,R_obs_pred,nugget)

plot(x=sigma_obs,y=y_obs,type="p",col="red",lty = 1,xlim=c(0.5,2),
     ylim = c(-3,48),cex.axis=1.8,cex=2.8,lwd=1.8,font=2,font.axis=2,cex.lab=1.8)
points(x=sigma_pred,y=y_pred,type="l",col="black",lty=2,lwd=1.8)
points(x=sigma_pred,y=hat_y_pred,type="l",col="blue",lty=2,lwd=1.8)
points(x=sigma_pred,y=hat_y_pred + 1.96*hat_sigma_pred,type="l",col="green",lty=2,lwd=1.8)
points(x=sigma_pred,y=hat_y_pred - 1.96*hat_sigma_pred,type="l",col="green",lty=2,lwd=1.8)











