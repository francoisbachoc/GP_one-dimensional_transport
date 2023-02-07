rm( list=ls() )   
source("../functions.R")




#######################################################
#settings
#######################################################
seed=1
f_code = f_code_m_sd
ntraining = 75
nvalidation = 25
ntest= 500
nsample = 500
nugget = 10^(-8)
theta_inf = c(0.1,0.01,0.25,10^(-8))
theta_sup = c(100,10,1,0.1)
init_theta = c(1,1,1,10^(-4))
ind_active=c(1,2,3,4)
d = 1000
cov_matrix_function = function(theta,mq) { cov_matrix_exp_quantile(theta[1],theta[2],theta[3],mq,mq) + theta[4]*diag(dim(mq)[1]) }
cross_cov_matrix_function = function(theta,mq1,mq2) { cov_matrix_exp_quantile(theta[1],theta[2],theta[3],mq1,mq2) }
n_gradient_descent=10
generate_pdf = generate_pdf_m_sd_pert
param_pdf = c(0.3,0.7,0.001,0.2,1,0.2,5/2)
name="TEST_TIME_gaussian_2d_pert_kriging"


#######################################################
#execution
#######################################################
set.seed(seed)  

n = ntraining+nvalidation

#Generation of the input curves
m_pdf_obs = generate_pdf(n,d,param_pdf)
matplot(t(m_pdf_obs[1:10,]),type="l")
m_quantile_obs = from_pdf_to_quantile(m_pdf_obs)
matplot(t(m_quantile_obs[1:10,]),type="l")

#Generation of the prediction curves
m_pdf_pred = generate_pdf(ntest,d,param_pdf)
matplot(t(m_pdf_pred[1:10,]),type="l")
m_quantile_pred = from_pdf_to_quantile(m_pdf_pred)
matplot(t(m_quantile_pred[1:10,]),type="l")

#Generation of the input samples
m_input_kernel_obs = from_pdf_to_quantile(m_pdf_obs)

#Generation of the prediction samples
m_input_kernel_pred = from_pdf_to_quantile(m_pdf_pred)

#Generation of the observation outputs
y_obs = seq(from=0,to=0,length=n)
for (i in 1:n) {
  y_obs[i] = f_code(m_pdf_obs[i,])
}

#Generation of the predictand outputs
y_pred = seq(from=0,to=0,length=ntest)
for (i in 1:ntest) {
  y_pred[i] = f_code(m_pdf_pred[i,])
}

#Generation of the observation samples
mobservation = matrix(nrow=n,ncol=nsample,data=-1)
for (i in 1:n) {
  mobservation[i,] = m_quantile_obs[i,sample(x=1:d,size = nsample,replace=TRUE)]
}

#Generation of the prediction samples
mtest = matrix(nrow=ntest,ncol=nsample,data=-1)
for (i in 1:ntest) {
  mtest[i,] = m_quantile_pred[i,sample(x=1:d,size = nsample,replace=TRUE)]
}

a = proc.time()

#maximum likelihood
m_qsample_obs = mobservation
for (i in 1:dim(m_qsample_obs)[1]) {
  m_qsample_obs[i,] = sort(mobservation[i,])
}

res_ML = maximum_likelihood(init_theta=init_theta,theta_inf=theta_inf,theta_sup=theta_sup,
                            ind_active=ind_active,nugget=nugget,cov_matrix_function=cov_matrix_function,
                            mq_obs=m_qsample_obs,y_obs=y_obs,k=n_gradient_descent)  
hat_theta = res_ML$hat_theta

m_qsample_pred = mtest
for (i in 1:dim(m_qsample_pred)[1]) {
  m_qsample_pred[i,] = sort(mtest[i,])
}

R_obs = cov_matrix_function(hat_theta,m_qsample_obs)
R_pred = cov_matrix_function(hat_theta,m_qsample_pred)
R_obs_pred = cross_cov_matrix_function(hat_theta,m_qsample_obs,m_qsample_pred)
hat_y_pred = prediction(R_obs,R_obs_pred,y_obs,nugget)
hat_sigma_pred = predictive_sd(R_obs,R_pred,R_obs_pred,nugget)

b = proc.time()

#######################################################
#Save of the results
#######################################################
#save(m_pdf_obs,m_input_kernel_obs,m_pdf_pred,m_input_kernel_pred,y_obs,y_pred,res_ML,hat_theta,hat_y_pred,hat_sigma_pred,
#     file=paste0(name,".Rdata"))







