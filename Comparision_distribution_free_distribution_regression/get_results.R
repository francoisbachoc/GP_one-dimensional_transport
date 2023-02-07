rm( list=ls() )   
source("../functions.R")

vName = c("skewness_beta_distribution_regression","skewness_beta_distribution_regression_sample5000",
          "skewness_beta_kriging","skewness_beta_kriging_sample5000")
for (name in vName) {
  load(paste0(name,".Rdata"))
  pdf(file=paste0(name,".pdf"))
  plot(va,ya,type="l",xlab="",ylab="",cex.axis=2.5)
  points(atest,predtest,type="p",pch=8,xlab="",ylab="")
  dev.off()
}


vName = c("gaussian_2d_pert_kriging")
for (name in vName) {
  cat("*********************************** \n")
  cat(name," \n")
  cat("*********************************** \n")
  
  load(file=paste0(name,".Rdata"))
  
  cat("RMSE = ",RMSE(y_pred,hat_y_pred),"\n")
  cat("Q2 = ",Q2(y_pred,hat_y_pred),"\n")
  cat("log_predictive_probability = ",log_predictive_probability(y_pred,hat_y_pred,hat_sigma_pred),"\n")
  cat("confidence interval coverage 0.90 = ",confidence_interval_coverage(y_pred,hat_y_pred,hat_sigma_pred,0.90),"\n")
  
  plot(y_pred,hat_y_pred,type="p")
  points(y_pred,hat_y_pred-1.64*hat_sigma_pred,type="p",col="green")
  points(y_pred,hat_y_pred+1.64*hat_sigma_pred,type="p",col="green")
}


vName = c("gaussian_2d_pert_distribution_regression")
for (name in vName) {
  cat("*********************************** \n")
  cat(name," \n")
  cat("*********************************** \n")
  
  load(file=paste0(name,".Rdata"))
  
  cat("RMSE = ",RMSE(y_pred,hat_y_pred),"\n")
  cat("Q2 = ",Q2(y_pred,hat_y_pred),"\n")
  
  plot(y_pred,hat_y_pred,type="p")
}
