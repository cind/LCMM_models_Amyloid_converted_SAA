################################################################################################  Bootstrapping/CI Libraries  #####################################################################################################
library(foreach)
library(doParallel)
library(purrr)

################################################################################################  Bootstrapping Function  #####################################################################################################
lcmm_bootstrap_ci <- function(new_data, n_iterations, lcmm_data, name_of_biomarker) {
  
  #prepping for bootstrapping
  n_sample <- nrow(lcmm_data)
  n_rows <- nrow(new_data) # newdata for predictions
  B <- n_iterations
  boot_pred <- numeric(B) #for the predictions

  #bootstrap
  registerDoParallel(16)
  bootstrapped_list_values <- foreach(ind = 1:B) %dopar% {
    sprintf("%i:%s", ind, name_of_biomarker)
    set.seed(ind)
    boot <- sample(n_sample, n_sample, replace = TRUE)
    
    lcmm_data$biomarker <- lcmm_data[ , name_of_biomarker]
    
    fit.b <- lcmm::lcmm(biomarker ~ adjusted_new_time + age + age*adjusted_new_time + PTGENDER + PTEDUCAT + apoe, #link = c("5-quantsplines"),
                        random = ~adjusted_new_time, subject="RID", maxiter = 300, data = lcmm_data[boot,], link = "3-equi-splines")
    
    boot_pred[ind] <- lcmm::predictY(fit.b, new_data, var.time = "adjusted_new_time", draws = TRUE)
    boot_pred_fit <- as.data.frame(boot_pred[ind])
    
    return(boot_pred_fit)
  }
  
  return(bootstrapped_list_values)
  
}

lcmm_bootstrap_ci_cdrsb <- function(new_data, n_iterations, lcmm_data, name_of_biomarker) {
  
  #prepping for bootstrapping
  n_sample <- nrow(lcmm_data)
  n_rows <- nrow(new_data) # newdata for predictions
  B <- n_iterations
  boot_pred <- numeric(B) #for the predictions

  #bootstrap
  registerDoParallel(16)
  bootstrapped_list_values <- foreach(ind = 1:B) %dopar% {
    sprintf("%i:%s", ind, name_of_biomarker)
    set.seed(ind)
    boot <- sample(n_sample, n_sample, replace = TRUE)
    
    lcmm_data$biomarker <- lcmm_data[ , name_of_biomarker]
    
    fit.b <- lcmm::lcmm(biomarker ~ adjusted_new_time + age + age*adjusted_new_time + PTGENDER + PTEDUCAT + apoe + CDGLOBAL, #link = c("5-quantsplines"),
                        random = ~adjusted_new_time, subject="RID", maxiter = 300, data = lcmm_data[boot,], link = "3-equi-splines")
    
    boot_pred[ind] <- lcmm::predictY(fit.b, new_data, var.time = "adjusted_new_time", draws = TRUE)
    boot_pred_fit <- as.data.frame(boot_pred[ind])
    
    return(boot_pred_fit)
  }
  
  return(bootstrapped_list_values)
  
}

#function to estimate confidence intervals
confidence_interval <- function(vector, interval) {
  
  # Standard deviation of sample
  vec_sd <- sd(vector)
  
  # Sample size
  n <- length(vector)
  
  # Mean of sample
  vec_mean <- mean(vector)
  
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
  
  return(result)
  
}
