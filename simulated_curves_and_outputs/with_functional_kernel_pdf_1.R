rm( list=ls() )   
source("../functions.R")

#######################################################
#settings (change accross simulations)
#######################################################
seed=1
name="res_functional_kernel_pdf"

#######################################################
#settings (do  not change accross simulations)
#######################################################
n = 100
n_pred = 500
sigma_inf=0.1
sigma_sup=100
lc_inf=0.01
lc_sup=10
H_inf=0.25
H_sup=1
init_theta = c(1,1,1)
ind_active=c(1,2,3)
d = 100
nugget=10^(-8)
sigma_input = 2
lc_input = 0.4
nu_input = 3/2
threshold=0.5
f_code = function(v_pdf) {
  x = seq(from=0,to=1,length=length(v_pdf))
  m = mean( v_pdf * cos( 5*x ) )
  m^2 - m + exp(m)-1
} 
cov_matrix_function = function(theta,mq) { cov_matrix_exp_quantile(theta[1],theta[2],theta[3],mq,mq) }
n_gradient_descent=10


#######################################################
#execution
#######################################################
set.seed(seed)  

#Generation of the input curves
m_pdf_obs = generate_pdf_matern(sigma_input,lc_input,nu_input,n,d,threshold)
matplot(t(m_pdf_obs[50:60,]),type="l")
m_q_obs = from_pdf_to_quantile(m_pdf_obs)
matplot(t(m_q_obs[50:60,]))

#Generation of the prediction curves
m_pdf_pred = generate_pdf_matern(sigma_input,lc_input,nu_input,n_pred,d,threshold)
matplot(t(m_pdf_pred[50:60,]))
m_q_pred = from_pdf_to_quantile(m_pdf_pred)
matplot(t(m_q_pred[50:60,]))

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
res_ML = maximum_likelihood(init_theta=init_theta,theta_inf=c(sigma_inf,lc_inf,H_inf),theta_sup=c(sigma_sup,lc_sup,H_sup),
                               ind_active=ind_active,nugget=nugget,cov_matrix_function=cov_matrix_function,
                               mq_obs=m_pdf_obs,y_obs=y_obs,k=n_gradient_descent)  
hat_theta = res_ML$hat_theta
sigma=hat_theta[1]
lc=hat_theta[2]
H=hat_theta[3]


R_obs = cov_matrix_exp_quantile(sigma,lc,H,m_pdf_obs,m_pdf_obs)
R_pred = cov_matrix_exp_quantile(sigma,lc,H,m_pdf_pred,m_pdf_pred)
R_obs_pred = cov_matrix_exp_quantile(sigma,lc,H,m_pdf_obs,m_pdf_pred)
hat_y_pred = prediction(R_obs,R_obs_pred,y_obs,nugget)
hat_sigma_pred = predictive_sd(R_obs,R_pred,R_obs_pred,nugget)




#######################################################
#Save of the results
#######################################################
save(m_pdf_obs,m_q_obs,m_pdf_pred,m_q_pred,y_obs,y_pred,res_ML,hat_theta,hat_y_pred,hat_sigma_pred,
     file=paste0(name,"_",seed,".Rdata"))







