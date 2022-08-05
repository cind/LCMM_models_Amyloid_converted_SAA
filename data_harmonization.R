library(mgcv)
library(gratia)
library(ggseg)
library(stringr)
library(ggplot2)

test <- filtered_merge

#getting age for each RID
dem <- read.csv("C:\\Documents\\Data Tables\\Medical_History\\PTDEMOG.csv") %>%
  select(RID, PTDOBMM, PTDOBYY, PTGENDER) %>%
  mutate(PTDOBDD = 1) %>%
  mutate(DOB = paste(PTDOBYY, PTDOBMM, PTDOBDD, sep = "-")) %>%
  mutate(DOB = as.POSIXct(DOB, format = "%Y-%m-%d")) %>%
  unique()

# first checking for missing data and fixing the missing data for gender
test <- test %>%
  select(!(c(Visit_code.x, LONIID.x, STATUS, OVERALLQC, TEMPQC, FRONTQC, PARQC, INSULAQC, OCCQC, BGQC, CWMQC, VENTQC, VISCODE2, IMAGETYPE, oldest_date, newest_date, ADNI_years, update_stamp, COLPROT, HIPPOQC, IMAGEID, LHIPQC, RHIPQC))) %>%
  merge(dem, by = "RID", all.x = TRUE) %>%
  mutate(Gender = case_when(PTGENDER == 2 ~ "F",
                            PTGENDER == 1 ~ "M")) %>%
  group_by(RID) %>% 
  fill(Gender, .direction = "downup") %>%
  ungroup() %>%
  select(-c(Sex)) %>%
  filter(!(is.na(DOB))) %>%
  filter(!is.na(FIELD_STRENGTH)) #this cuts down on quite a few people, so maybe this data can still be found somewhere?
sum(is.na(test)) #if missing data exists, then impute/delete rows/columns with missing data and then try again

test$age <- time_length(difftime(test$EXAMDATE, test$DOB), "years")

test <- test %>%
  select(-c(PTDOBMM, PTDOBYY, PTDOBDD, PTGENDER, DOB))

#for now getting rid of manufacturer, model name, and software version, but might want to run a seperate analysis with these
test <- test %>%
  select(-c(Manufacturer, ManufacturerModelName, SoftwareVersion))

#get rid of na values in test
# getting the counts of NA for each region  
na_count <- sapply(test, function(y) sum(length(which(is.na(y)))))

#getting rid of variables with mostly NA values
test <- test[, which(colMeans(!is.na(test)) > 0.9)] # only keeping variables with at least 90% non-na values

#now can get rid of NA's in the data - now get rid of the rows that contribute only some missing values. If I did this at an earlier step, there would be way too few rows
test <- test %>%
  na.omit() %>%
  unique() #now it is down to 1031 cases

# making batch a factor
test$batch <- as.factor(test$FIELD_STRENGTH)
batch <- droplevels(as.factor(test$batch)) # drop levels in case there are still factors from previous versions of the dataset
# check for batches with only one observation
min(table(batch)) # check to see that there is at least 2 observations of each batch 

# number of batches
m <- nlevels(batch)
# row IDs for each batch 
batches <- lapply(levels(batch), function(x) which(batch==x))
# number of observations for each batch
ni <- sapply(batches, length)

#making all other factors that need to be factors, factors
test$RID <- as.factor(test$RID)

# isolating the numeric features (the brain regions)
all_features <- test[, c(3:314)]

# feature names
featurenames <- names(all_features)

# number of features
V <- length(featurenames)
# total number of observations
L <- nrow(test) #or is this nrow(test)

# defining the covariates
all_covariates <- test[, c("RID", "age", "batch", "EXAMDATE", "Gender")] #"DXSUM" excluded
covariatenames <- names(all_covariates)


# defining the formula
model_formula = "~ batch + s(age) + s(RID, bs = 're')" #should gender be in here?

data <- as.data.frame(cbind(all_features, all_covariates))
colnames(data) <- c(featurenames, covariatenames)

#potentially here check for normality of columns we want harmonized?

#now making models for all features
modlist    <- list()

sigma_estimates <- rep(NA, V)
predicted <- matrix(nrow=L, ncol=V)
batch_effects <- matrix(nrow=(m-1), ncol=V)

for (j in 1:length(featurenames)) {
  if(TRUE) {
    cat("[ComGam] ", "Fit model", j, "out of", length(featurenames), "\n")
  }
  gammod              <- gam(as.formula(paste(featurenames[j], model_formula, sep = "")) ,   
                             method = "REML", data = data)
  
  modlist[[j]]        <- gammod
  
  #now getting the sigma (standard deviation) estimates for later standardization
  variance <- as.data.frame(variance_comp(gammod)) #looking at scale std dev
  sigma_estimates[j] <- variance[variance$component == 'scale', 'std_dev']
  
  # getting batch effects for each variable by isolating the coefficient for each batch
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
data_std <- (data[,featurenames] - predicted + batch_effects_expanded) / sigmas

######################################################################################################################

##############################
# method of moments to estimate hyperparameters (Hyperparameters  are estimated from standardized data using the method of moments)
##############################
gammahat <- matrix(nrow=m, ncol=V) #location parameters before Empirical Bayes
delta2hat <- matrix(nrow=m, ncol=V) #shift parameters before Empirical Bayes
for (i in 1:m){ # begin loop over batches
  gammahat[i,] <- colMeans(data_std[batches[[i]],])
  delta2hat[i,] <- apply(data_std[batches[[i]],], 2, var)
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
gammastarhat0 <- matrix(nrow=m, ncol=V)
for (v in 1:V){ # begin loop over features
  gammastarhat0[,v] <- ((ni * tau2bar * gammahat[,v]) + (delta2hat[,v] * gammabar))/((ni * tau2bar) + delta2hat[,v])
} # end loop over features
delta2starhat0 <- matrix(nrow=m, ncol=V)
for (v in 1:V){ # begin loop over features
  for(i in 1:m){ # begin loop over batches
    zminusgammastarhat2 <- sum((data_std[batches[[i]],v] - gammastarhat0[i,v])^2)
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
      zminusgammastarhat2 <- sum((data_std[batches[[i]],v] - gammastarhat[i,v,(b-1)])^2)
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
data_combat <- (sigmas/sqrt(delta2starhat_expanded))*(data_std - gammastarhat_expanded) + predicted - batch_effects_expanded

##############################
# label the data
##############################
# add IDs, time variable, and batch variable to data_combat
data_combat <- cbind(data[,c("RID", "EXAMDATE", "age", "batch", "Gender")], data_combat) #"DXSUM"
# add names
colnames(data_combat) <- c("RID", "EXAMDATE", "age", "batch", "Gender", paste0(featurenames, '.combat')) #"DXSUM"
colnames(gammahat) <- featurenames
colnames(delta2hat) <- featurenames
colnames(gammastarhat_final) <- featurenames
colnames(delta2starhat_final) <- featurenames
rownames(gammahat) <- levels(batch)
rownames(delta2hat) <- levels(batch)
rownames(gammastarhat_final) <- levels(batch)
rownames(delta2starhat_final) <- levels(batch)
