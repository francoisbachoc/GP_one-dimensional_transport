rm( list=ls() )   
source("../functions.R")

#######################################################
#settings (change accross simulations)
#######################################################
seed=1
name="res_PCA"

#######################################################
#settings (do  not change accross simulations)
#######################################################
n = 100
n_pred = 500
order=15
init_theta = c(1,rep(1,order),1,0)
theta_inf = c(0.1,rep(0.01,order),0.25,0)
theta_sup = c(100,rep(20,order),1,0.1)
ind_active=seq(from=1,to=order+3,length=order+3)
d = 100
nugget=10^(-6)
f_code = f_code_m_sd
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
generate_pdf = generate_pdf_m_sd_pert
param_pdf = c(0.3,0.7,0.001,0.2,1,0.2,5/2)
n_gradient_descent=30

#######################################################
#execution
#######################################################
set.seed(seed)  

#Generation of the input curves
m_pdf_obs = generate_pdf(n,d,param_pdf)
matplot(t(m_pdf_obs[1:30,]),type="l")
m_input_kernel_obs = from_pdf_to_vector_PCA(m_pdf_obs,order)$m_proj
m_component = from_pdf_to_vector_PCA(m_pdf_obs,order)$m_component
matplot(t(m_input_kernel_obs[1:30,]),type="l")

#Generation of the prediction curves
m_pdf_pred = generate_pdf(n_pred,d,param_pdf)
matplot(t(m_pdf_pred[1:30,]),type="l")
m_input_kernel_pred = (1/d) * m_pdf_pred %*% t( m_component )
matplot(t(m_input_kernel_pred[1:30,]),type="l")

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
                            mq_obs=m_input_kernel_obs,y_obs=y_obs,k=n_gradient_descent)  
hat_theta = res_ML$hat_theta



R_obs = cov_matrix_function(hat_theta,m_input_kernel_obs)
R_pred = cov_matrix_function(hat_theta,m_input_kernel_pred)
R_obs_pred = cross_cov_matrix_function(hat_theta,m_input_kernel_obs,m_input_kernel_pred)
hat_y_pred = prediction(R_obs,R_obs_pred,y_obs,nugget)
hat_sigma_pred = predictive_sd(R_obs,R_pred,R_obs_pred,nugget)



#######################################################
#Save of the results
#######################################################
save(m_pdf_obs,m_input_kernel_obs,m_pdf_pred,m_input_kernel_pred,y_obs,y_pred,res_ML,hat_theta,hat_y_pred,hat_sigma_pred,
     file=paste0(name,"_order",order,"_",seed,".Rdata"))







