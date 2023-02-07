rm( list=ls() )   
source("../functions.R")

v_name = c("res_functional_kernel_1.Rdata","res_legendre_order4_1.Rdata","res_legendre_order8_1.Rdata",
           "res_legendre_order9_1.Rdata",
           "res_legendre_order12_1.Rdata","res_legendre_order20_1.Rdata","res_legendre_order30_1.Rdata",
           "res_PCA_order4_1.Rdata","res_PCA_order6_1.Rdata","res_PCA_order12_1.Rdata")

for (name in v_name) {
  cat("*********************************** \n")
  cat(name," \n")
  cat("*********************************** \n")
  
  load(file=name)
  
  cat("RMSE = ",RMSE(y_pred,hat_y_pred),"\n")
  cat("Q2 = ",Q2(y_pred,hat_y_pred),"\n")
  cat("log_predictive_probability = ",log_predictive_probability(y_pred,hat_y_pred,hat_sigma_pred),"\n")
  cat("confidence interval coverage 0.90 = ",confidence_interval_coverage(y_pred,hat_y_pred,hat_sigma_pred,0.90),"\n")
  
  plot(y_pred,hat_y_pred,type="p")
  points(y_pred,hat_y_pred-1.64*hat_sigma_pred,type="p",col="green")
  points(y_pred,hat_y_pred+1.64*hat_sigma_pred,type="p",col="green")
}







