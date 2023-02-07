rm( list=ls() )   
set.seed(1)  
source("../functions.R")

#######################################################
#settings
#######################################################
sigma_obs=c(0.5,0.75,1,1.5,2)
sigma_pred = seq(from=0.5,to=2,length=100)
sigma_inf=0.1
sigma_sup=100
lc_inf=0.05
lc_sup=5
H_inf=0.25
H_sup=1
init_theta = c(1,10,0.5)
ind_active=c(1,2,3)
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

cov_matrix_function = function(theta,mq) { cov_matrix_exp_quantile(theta[1],theta[2],theta[3],mq,mq) }
hat_theta = maximum_likelihood(init_theta=init_theta,theta_inf=c(sigma_inf,lc_inf,H_inf),theta_sup=c(sigma_sup,lc_sup,H_sup),
                               ind_active=ind_active,nugget=nugget,cov_matrix_function=cov_matrix_function,
                               mq_obs=M_obs,y_obs=y_obs)  
sigma=hat_theta[1]
lc=hat_theta[2]
H=hat_theta[3]


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











