library(mgcv)
library(gratia)
library(ggseg)
library(stringr)
library(ggplot2)

long_combat_gam <- function(batch, all_features, all_covariates, ID, model_formula){ #batch is the batch variable, all_features is the names of the variables that you plan to harmonize, also include your ID variable (in my case, RID)
  
  # number of batches
  m <- nlevels(batch)
  # row IDs for each batch 
  batches <- lapply(levels(batch), function(x) which(batch==x))
  # number of observations for each batch
  ni <- sapply(batches, length)
  
  #getting the names of all the features
  all_featurenames <- names(all_features) #all_features
  
  # number of features
  V <- length(all_featurenames)
  
  # getting the names of all the covariates
  all_covariatenames <- names(all_covariates) #all_covariates
  
  #combining the data
  all_data       <- as.data.frame(cbind(all_features, all_covariates))
  colnames(all_data) <- c(all_featurenames, all_covariatenames)
  
  # total number of observations
  L <- nrow(all_data) #or is this nrow(all_test)
  
  #now making models for all features
  modlist    <- list()
  
  sigma_estimates <- rep(NA, V)
  predicted <- matrix(nrow=L, ncol=V)
  batch_effects <- matrix(nrow=(m-1), ncol=V)
  
  for (j in 1:length(all_featurenames)) {
    if(TRUE) {
      cat("[ComGam] ", "Fit model", j, "out of", length(all_featurenames), "\n")
    }
    gammod              <- gam(as.formula(paste(all_featurenames[j], model_formula, sep = "")) ,   
                               method = "REML", data = all_data)
    
    modlist[[j]]        <- gammod
    
    #now getting the sigma (standard deviation) estimates for later standardization
    variance <- as.data.frame(variance_comp(gammod)) #looking at scale std dev
    sigma_estimates[j] <- variance[variance$component == 'scale', 'std_dev']
    
    # getting batch effects for each variable by isolating the coefficient for each batch
    #temp_batch <- paste("batch", as.character(batch), sep = '')
    coef_info <- as.data.frame(summary(modlist[[j]])[1]$p.coeff) %>%
      tibble::rownames_to_column(var = "coefficient_name") %>%
      filter(grepl("batch", coefficient_name))
    
    batch_effects[,j] <- coef_info$`summary(modlist[[j]])[1]$p.coeff`
    
    # getting the predicted values
    predicted[,j] <- predict.gam(modlist[[j]])
  } # end loop over features
  
  # create a L*V matrix of sigma estimates
  sigmas <- matrix(rep(sigma_estimates, each=L), nrow=L, ncol=V)
  
  # calculate the gamma1 hats  (estimates of additive batch effects)
  gamma1hat <- -(ni[2:m] %*% batch_effects)/L
  
  # add gamma1hat to the rest of the batch effect table
  batch_effects_adjusted <- sweep(batch_effects, 2, gamma1hat, FUN='+')
  
  # add gamma1hat as the top row
  batch_effects_adjusted <- rbind(gamma1hat, batch_effects_adjusted)
  
  # expand the adjusted batch effects to all timepoints
  batch_effects_expanded <- matrix(nrow=L, ncol=V)
  
  for(i in 1:m){ # begin loop over batches
    batch_effects_expanded[batches[[i]],] <- matrix(
      rep(batch_effects_adjusted[i,],length(batches[[i]])),
      ncol=V, byrow=TRUE) 
  } # end loop over batches
  
  # standardizing the data with these new estimates/predictions/etc.
  data_std <- (all_data[,all_featurenames] - predicted + batch_effects_expanded) / sigmas
  
  ##############################
  # method of moments to estimate hyperparameters (Hyperparameters  are estimated from standardized data using the method of moments)
  ##############################
  gammahat <- matrix(nrow=m, ncol=V) #location parameters before Empirical Bayes
  delta2hat <- matrix(nrow=m, ncol=V) #shift parameters before Empirical Bayes
  
  
  for (i in 1:m){ # begin loop over batches
    gammahat[i,] <- colMeans(data_std[batches[[i]],], na.rm = TRUE)
    delta2hat[i,] <- apply(data_std[batches[[i]],], 2, var, na.rm = TRUE)
  } # end loop over batches
  gammabar <- rowMeans(gammahat)
  tau2bar <- apply(gammahat, 1, var)
  Dbar <- rowMeans(delta2hat)
  S2bar <- apply(delta2hat, 1, var)
  # inverse gamma parameters
  lambdabar <- (Dbar^2 + 2*S2bar) / S2bar
  thetabar <- (Dbar^3 + Dbar*S2bar) / S2bar
  
  ##############################
  # empirical Bayes to estimate batch effects
  ##############################
  #if (verbose) cat('[longCombat] using empirical Bayes to estimate batch effects...\n')
  #if (verbose) cat('[longCombat] initializing...\n')
  # get initial estimates
  gammastarhat0 <- matrix(nrow=m, ncol=V) #-0.4675299, 0.1807892, 0.5481943t66
  for (v in 1:V){ # begin loop over features
    gammastarhat0[,v] <- ((ni * tau2bar * gammahat[,v]) + (delta2hat[,v] * gammabar))/((ni * tau2bar) + delta2hat[,v])
  } # end loop over features
  delta2starhat0 <- matrix(nrow=m, ncol=V) #0.7304484, 0.9766520, 0.6830730
  for (v in 1:V){ # begin loop over features
    for(i in 1:m){ # begin loop over batches
      zminusgammastarhat2 <- sum((data_std[batches[[i]],v] - gammastarhat0[i,v])^2, na.rm =TRUE)
      delta2starhat0[i,v] <- (thetabar[i] + 0.5*zminusgammastarhat2) / (ni[i]/2 + lambdabar[i] - 1)
    } # end loop over features
  } # end loop over batches
  #setting a default number of iterations to 30
  niter = 30
  
  # iterate
  gammastarhat <- array(dim=c(m, V, (niter+1)))
  gammastarhat[,,1] <- gammastarhat0
  delta2starhat <- array(dim=c(m, V, (niter+1)))
  delta2starhat[,,1] <- delta2starhat0
  for(b in 2:(niter+1)){ # begin loop over iterations
    #if (verbose) cat(paste0('[longCombat] starting EM algorithm iteration ', (b-1), '\n')) 
    for (v in 1:V){ # begin loop over features
      gammastarhat[,v,b] <- ((ni * tau2bar * gammahat[,v]) + (delta2starhat[,v,(b-1)] * gammabar))/((ni * tau2bar) + delta2starhat[,v,(b-1)])
      for(i in 1:m){ # begin loop over batches
        zminusgammastarhat2 <- sum((data_std[batches[[i]],v] - gammastarhat[i,v,(b-1)])^2, na.rm = TRUE)
        delta2starhat[i,v,b] <- (thetabar[i] + 0.5*zminusgammastarhat2) / (ni[i]/2 + lambdabar[i] - 1)
      } # end loop over batches
    } # end loop over features
  } # end loop over iterations
  # save final result
  gammastarhat_final <- gammastarhat[,,niter+1] #location parameters after Empirical Bays
  delta2starhat_final <- delta2starhat[,,niter+1] #shift parameters after Empirical Bayes
  
  ##############################
  # adjust data for batch effects
  ##############################
  #if (verbose) cat('[longCombat] adjusting data for batch effects\n')
  # repeat each row the correct number of times
  gammastarhat_expanded <- matrix(nrow=L, ncol=V)
  delta2starhat_expanded <- matrix(nrow=L, ncol=V)
  for(i in 1:m){ # loop over batches
    gammastarhat_expanded[batches[[i]],] <- matrix(
      rep(gammastarhat_final[i,],length(batches[[i]])),
      ncol=V, byrow=TRUE) 
    delta2starhat_expanded[batches[[i]],] <- matrix(
      rep(delta2starhat_final[i,],length(batches[[i]])),
      ncol=V, byrow=TRUE) 
  } # end loop over batches
  # do ComBat 
  all_data_combat <- (sigmas/sqrt(delta2starhat_expanded))*(data_std - gammastarhat_expanded) + predicted - batch_effects_expanded
  
  ##############################
  # label the data - have to figure out how to get this part adjusted (do not want to )
  ##############################
  #labeling post-combat variables
  colnames(all_data_combat) <-  paste0(colnames(all_data_combat), '.combat') 
  #combining the data with previous covariates
  all_data_combat_covariates <- cbind(all_covariates, all_data_combat) #all_covariates
  
  colnames(gammahat) <- all_featurenames
  colnames(delta2hat) <- all_featurenames
  colnames(gammastarhat_final) <- all_featurenames
  colnames(delta2starhat_final) <- all_featurenames
  rownames(gammahat) <- levels(batch)
  rownames(delta2hat) <- levels(batch)
  rownames(gammastarhat_final) <- levels(batch)
  rownames(delta2starhat_final) <- levels(batch)
  
  return(all_data_combat_covariates, gammastarhat_final, delta2starhat_final) #gammastarhat is the location, deltastarhat is the shift
  
}
