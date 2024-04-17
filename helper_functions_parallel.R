################################################################################################  Bootstrapping Function  #####################################################################################################
library(foreach)
library(doParallel)

lcmm_helper <- function(ind, new_data, n_sample, lcmm_data, name_of_biomarker, boot_pred, boot_derivs) {
    sprintf("%i:%s", ind, name_of_biomarker)
    boot <- sample(n_sample, n_sample, replace = TRUE)

    # lcmm_data$biomarker <- lcmm_data[,"name_of_biomarker"]
    #lcmm_data$biomarker <- lcmm_data[ , as.character(substitute(name_of_biomarker))]
    lcmm_data$biomarker <- lcmm_data[ , name_of_biomarker]

    fit.b <- lcmm(biomarker ~ adjusted_new_time + age + age*adjusted_new_time + PTGENDER + PTEDUCAT + apoe, #link = c("5-quantsplines"),
                  random = ~adjusted_new_time, subject="RID", maxiter = 300, data = lcmm_data[boot,], link = "3-equi-splines")

    boot_pred[ind] <- predictY(fit.b, new_data, var.time = "adjusted_new_time", draws = TRUE)
    boot_pred_fit <- as.data.frame(boot_pred[ind])
    boot_derivs[,ind] <- diff(boot_pred_fit$Ypred_50)/diff(new_data$adjusted_new_time)

    # Do we need to return data here ???
  }

lcmm_bootstrap_ci <- function(new_data, n_iterations, lcmm_data, name_of_biomarker) {

  #prepping for bootstrapping
  n_sample <- nrow(lcmm_data)
  n_rows <- nrow(new_data) # newdata for predictions
  B <- n_iterations
  boot_derivs <- matrix(nrow = n_rows-1, ncol=B) #for the finite difference
  boot_derivs <- as.data.frame(boot_derivs)
  boot_pred <- numeric(B) #for the predictions

  set.seed(121)
  # No parallel, test the code.
  #i <- 1
  #lcmm_helper(i, new_data, n_sample, lcmm_data, name_of_biomarker, boot_pred, boot_derivs)

  registerDoParallel(16)
  print("running loop")
  ind <- 1
  foreach (i<-1:B) %dopar% {
    lcmm_helper(ind, new_data, n_sample, lcmm_data, name_of_biomarker, boot_pred, boot_derivs)
    ind <- ind + 1
    sprintf("%i", ind)
  }

  boot_pred_data <- as.data.frame(boot_pred)
  data_2.5 <- boot_pred_data[ , grepl( "Ypred_2.5" , names( boot_pred_data ) ) ]
  data_97.5 <- boot_pred_data[ , grepl( "Ypred_97.5" , names( boot_pred_data ) ) ]
  ci_2.5 <- rowMeans(data_2.5)
  ci_97.5 <- rowMeans(data_97.5)

  #getting bootstrapped ci's from predictions
  data_preds <- boot_pred_data[ , grepl( "Ypred_50" , names( boot_pred_data ) ) ]
  median_preds_ci_2.5 <- apply(as.matrix(data_preds), 1, function(x){median(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})[1,]
  median_preds_ci_97.5 <- apply(as.matrix(data_preds), 1, function(x){median(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})[2,]
  mean_pred_ci_2.5 <- apply(as.matrix(data_preds), 1, function(x){mean(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})[1,]
  mean_pred_ci_97.5 <- apply(as.matrix(data_preds), 1, function(x){mean(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})[2,]


  #now for the slopes
  ci_diffs_2.5 <- append(apply(as.matrix(boot_derivs), 1, function(x){mean(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})[1,], NA)
  ci_diffs_97.5 <- append(apply(as.matrix(boot_derivs), 1, function(x){mean(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})[2,], NA)

  #now putting into a data frame
  bootstrapped_data <- data.frame(new_data, ci_2.5, ci_97.5, ci_diffs_2.5, ci_diffs_97.5, median_preds_ci_2.5, median_preds_ci_97.5, mean_pred_ci_2.5, mean_pred_ci_97.5)

  bootstrapped_data$signif <- case_when(bootstrapped_data$ci_diffs_2.5 > 0 & bootstrapped_data$ci_diffs_97.5 > 0 ~ "sig",
                                        bootstrapped_data$ci_diffs_2.5 < 0 & bootstrapped_data$ci_diffs_97.5 < 0 ~ "sig",
                                        TRUE ~ "not sig")

  bootstrapped_data$biomarker <- name_of_biomarker

  return(bootstrapped_data)

}

################################################################################################  Smoothing Function  #####################################################################################################

smoothing_bootstrapped_data <- function(bootstrapped_data, lcmm_predictions, span_parameter) {

  bootstrapped_data$predictions <- as.data.frame(lcmm_predictions$pred)$Ypred_50

  #first smoothing the bootstrapped ci's and predictions and CI's from curve
  bootstrapped_data$smooth_lower_ci <- predict(loess(bootstrapped_data$mean_pred_ci_2.5 ~ as.numeric(bootstrapped_data$adjusted_new_time), span = span_parameter))

  bootstrapped_data$smooth_upper_ci <- predict(loess(bootstrapped_data$mean_pred_ci_97.5 ~ as.numeric(bootstrapped_data$adjusted_new_time), span = span_parameter))

  bootstrapped_data$smoothed_predictions <- predict(loess(as.data.frame(lcmm_predictions$pred)$Ypred_50 ~ as.numeric(bootstrapped_data$adjusted_new_time), span = span_parameter))

  # now smoothing the CI's provided by LCMM
  bootstrapped_data$smooth_lower_lcmm_ci <- predict(loess(bootstrapped_data$ci_2.5 ~ as.numeric(bootstrapped_data$adjusted_new_time), span = span_parameter))

  bootstrapped_data$smooth_upper_lcmm_ci <- predict(loess(bootstrapped_data$ci_97.5 ~ as.numeric(bootstrapped_data$adjusted_new_time), span = span_parameter))

  #now scaling for plotting purposes
  bootstrapped_data$lower_scaled = (bootstrapped_data$smooth_lower_ci - mean(bootstrapped_data$smoothed_predictions) ) / sd(bootstrapped_data$smoothed_predictions)

  bootstrapped_data$upper_scaled = (bootstrapped_data$smooth_upper_ci - mean(bootstrapped_data$smoothed_predictions) ) / sd(bootstrapped_data$smoothed_predictions)

  bootstrapped_data$lower_lcmm_scaled = (bootstrapped_data$smooth_lower_lcmm_ci - mean(bootstrapped_data$smoothed_predictions) ) / sd(bootstrapped_data$smoothed_predictions)

  bootstrapped_data$upper_lcmm_scaled = (bootstrapped_data$smooth_upper_lcmm_ci - mean(bootstrapped_data$smoothed_predictions) ) / sd(bootstrapped_data$smoothed_predictions)

  bootstrapped_data$predictions_scaled <- scale(bootstrapped_data$smoothed_predictions)

  #now smoothing the finite differences confidence intervals

  bootstrapped_data$smoothed_diff_lower_ci <- append(predict(loess(bootstrapped_data$ci_diffs_2.5 ~ as.numeric(bootstrapped_data$adjusted_new_time), span = span_parameter)), NA)

  bootstrapped_data$smoothed_diff_upper_ci <- append(predict(loess(bootstrapped_data$ci_diffs_97.5 ~ as.numeric(bootstrapped_data$adjusted_new_time), span = span_parameter)), NA)

  bootstrapped_data$signif_diff_smoothed <- case_when(bootstrapped_data$smoothed_diff_lower_ci > 0 & bootstrapped_data$smoothed_diff_upper_ci > 0 ~ "sig",
                                                      bootstrapped_data$smoothed_diff_lower_ci < 0 & bootstrapped_data$smoothed_diff_upper_ci < 0 ~ "sig",
                                                      TRUE ~ "not sig")

  #now making first instance of significance
  bootstrapped_data$sign <- case_when(bootstrapped_data$smoothed_diff_lower_ci > 0 & bootstrapped_data$smoothed_diff_upper_ci > 0 ~ "pos",
                                      TRUE ~ "not pos")
  # temp_data <- bootstrapped_data[bootstrapped_data$sign == "pos",][1, ]
  #
  # bootstrapped_data$y_line <- temp_data$adjusted_new_time
  if(nrow(bootstrapped_data[bootstrapped_data$sign == "pos",]) > 0){

    temp_data <- bootstrapped_data[bootstrapped_data$sign == "pos",][1, ]

    bootstrapped_data$y_line <- temp_data$adjusted_new_time

    } else {

    bootstrapped_data$y_line <- NA

    }

  return(bootstrapped_data)
}
