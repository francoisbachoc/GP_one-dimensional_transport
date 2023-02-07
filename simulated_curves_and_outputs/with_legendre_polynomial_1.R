rm( list=ls() )   
source("../functions.R")

#######################################################
#settings (change accross simulations)
#######################################################
seed=1
name="res_legendre"

#######################################################
#settings (do  not change accross simulations)
#######################################################
n = 100
n_pred = 500
order=5
init_theta = c(1,rep(1,order+1),1,0)
theta_inf = c(0.1,rep(0.01,order+1),0.25,0)
theta_sup = c(100,rep(10,order+1),1,0.1)
ind_active=seq(from=1,to=order+3,length=order+3)
d=100
nugget=10^(-6)
sigma_input = 2
lc_input = 0.4
nu_input = 3/2
threshold=0.5
f_code = function(v_pdf) {
x = seq(from=0,to=1,length=length(v_pdf))
m = mean( v_pdf * cos( 5*x ) )
m^2 - m + exp(m)-1
} 
cov_matrix_function = function(theta,X) { 
n = dim(X)[1]
d = dim(X)[2]
  for (i in 1:d) {
    X[,i] = X[,i] / theta[i+1]
  }
cov_matrix_exp_quantile(theta[1],1,theta[d+2],X,X) + theta[d+3]*diag(n)
}
cross_cov_matrix_function = function(theta,X,Y) { 
  n = dim(X)[1]
  d = dim(X)[2]
  for (i in 1:d) {
    X[,i] = X[,i] / theta[i+1]
    Y[,i] = Y[,i] / theta[i+1]
  }
  cov_matrix_exp_quantile(theta[1],1,theta[d+2],X,Y)
}
n_gradient_descent=10


#######################################################
#execution
#######################################################
set.seed(seed)  

#Generation of the input curves
m_pdf_obs = generate_pdf_matern(sigma_input,lc_input,nu_input,n,d,threshold)
matplot(t(m_pdf_obs[50:60,]))
m_proj_obs = from_pdf_to_vector_legendre(m_pdf_obs,order)
matplot(t(m_proj_obs[50:60,]))

#Generation of the prediction curves
m_pdf_pred = generate_pdf_matern(sigma_input,lc_input,nu_input,n_pred,d,threshold)
matplot(t(m_pdf_pred[50:60,]))
m_proj_pred = from_pdf_to_vector_legendre(m_pdf_pred,order)
matplot(t(m_proj_pred[50:60,]))

#Generation of the observation outputs
y_obs = seq(from=0,to=0,length=n)
for (i in 1:n) {
y_obs[i] = f_code(m_pdf_obs[i,])
}

#Generation of the predictand outputs
y_pred = seq(from=0,to=0,length=n_pred)
for (i in 1:n_pred) {
y_pred[i] = f_code(m_pdf_pred[i,])
}

#maximum likelihood
res_ML = maximum_likelihood(init_theta=init_theta,theta_inf=theta_inf,theta_sup=theta_sup,
                             ind_active=ind_active,nugget=nugget,cov_matrix_function=cov_matrix_function,
                             mq_obs=m_proj_obs,y_obs=y_obs,k=n_gradient_descent)  
hat_theta = res_ML$hat_theta


R_obs = cov_matrix_function(hat_theta,m_proj_obs)
R_pred = cov_matrix_function(hat_theta,m_proj_pred)
R_obs_pred = cross_cov_matrix_function(hat_theta,m_proj_obs,m_proj_pred)
hat_y_pred = prediction(R_obs,R_obs_pred,y_obs,nugget)
hat_sigma_pred = predictive_sd(R_obs,R_pred,R_obs_pred,nugget)




#######################################################
#Save of the results
#######################################################
save(m_pdf_obs,m_proj_obs,m_pdf_pred,m_proj_pred,y_obs,y_pred,res_ML,hat_theta,hat_y_pred,hat_sigma_pred,
     file=paste0(name,"_",seed,".Rdata"))







