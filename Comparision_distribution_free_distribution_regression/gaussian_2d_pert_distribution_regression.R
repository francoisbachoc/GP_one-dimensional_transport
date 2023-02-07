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
nDistEval=4096
nbh=100
hlim=c(0,1)
blim=c(0,1)
d = 1000
generate_pdf = generate_pdf_m_sd_pert
param_pdf = c(0.3,0.7,0.001,0.2,1,0.2,5/2)
name="gaussian_2d_pert_distribution_regression"


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

mtraining = mobservation[1:ntraining,]
mvalidation = mobservation[(ntraining+1):(ntraining+nvalidation),]
ytraining = y_obs[1:ntraining]
yvalidation = y_obs[(ntraining+1):(ntraining+nvalidation)]

bh = select_bh_kernel_regression(mtraining = mtraining,ytraining = ytraining,
                                 mvalidation = mvalidation,yvalidation = yvalidation,
                                 nDistEval=nDistEval,nbh=nbh,hlim=hlim,blim=blim)

hat_y_pred = kernel_regression(rbind(mtraining,mvalidation),c(ytraining,yvalidation),
                             mtest,h=bh[2],b=bh[1],nDistEval=nDistEval)

#######################################################
#Save of the results
#######################################################
save(m_pdf_obs,m_input_kernel_obs,m_pdf_pred,m_input_kernel_pred,y_obs,y_pred,hat_y_pred,
     file=paste0(name,".Rdata"))







