rm( list=ls() )   
source("../functions.R")


cat("*********************************** \n")
cat("functional kernel quantile \n")
cat("*********************************** \n")
load(file="res_functional_kernel_quantile_1.Rdata")

cat("RMSE = ",RMSE(y_pred,hat_y_pred),"\n")
cat("Q2 = ",Q2(y_pred,hat_y_pred),"\n")
cat("log_predictive_probability = ",log_predictive_probability(y_pred,hat_y_pred,hat_sigma_pred),"\n")
cat("confidence interval coverage 0.90 = ",confidence_interval_coverage(y_pred,hat_y_pred,hat_sigma_pred,0.90),"\n")

plot(y_pred,hat_y_pred,type="p")
points(y_pred,hat_y_pred-1.64*hat_sigma_pred,type="p",col="green")
points(y_pred,hat_y_pred+1.64*hat_sigma_pred,type="p",col="green")

cat("*********************************** \n")
cat("functional kernel pdf \n")
cat("*********************************** \n")
load(file="res_functional_kernel_pdf_1.Rdata")

cat("RMSE = ",RMSE(y_pred,hat_y_pred),"\n")
cat("Q2 = ",Q2(y_pred,hat_y_pred),"\n")
cat("log_predictive_probability = ",log_predictive_probability(y_pred,hat_y_pred,hat_sigma_pred),"\n")
cat("confidence interval coverage 0.90 = ",confidence_interval_coverage(y_pred,hat_y_pred,hat_sigma_pred,0.90),"\n")

plot(y_pred,hat_y_pred,type="p")
points(y_pred,hat_y_pred-1.64*hat_sigma_pred,type="p",col="green")
points(y_pred,hat_y_pred+1.64*hat_sigma_pred,type="p",col="green")



cat("*********************************** \n")
cat("legendre \n")
cat("*********************************** \n")
load(file="res_legendre_1.Rdata")

cat("RMSE = ",RMSE(y_pred,hat_y_pred),"\n")
cat("Q2 = ",Q2(y_pred,hat_y_pred),"\n")
cat("log_predictive_probability = ",log_predictive_probability(y_pred,hat_y_pred,hat_sigma_pred),"\n")
cat("confidence interval coverage 0.90 = ",confidence_interval_coverage(y_pred,hat_y_pred,hat_sigma_pred,0.90),"\n")

plot(y_pred,hat_y_pred,type="p")
points(y_pred,hat_y_pred-1.64*hat_sigma_pred,type="p",col="green")
points(y_pred,hat_y_pred+1.64*hat_sigma_pred,type="p",col="green")




