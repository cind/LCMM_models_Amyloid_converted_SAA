################################################################################################  Bootstrapping Function  #####################################################################################################

lcmm_bootstrap_ci <- function(new_data, n_iterations, lcmm_data, name_of_biomarker) {
  
  #prepping for bootstrapping
  n_sample <- nrow(lcmm_data)
  n_rows <- nrow(new_data) # newdata for predictions
  B <- n_iterations
  boot_derivs <- matrix(nrow = n_rows-1, ncol=B) #for the finite difference
  boot_derivs <- as.data.frame(boot_derivs)
  boot_pred <- numeric(B) #for the predictions
  
  #making a for loop for the iterations
  set.seed(121)
  for (i in 1:B) {
    boot <- sample(n_sample, n_sample, replace = TRUE)
    
    # lcmm_data$biomarker <- lcmm_data[,"name_of_biomarker"]
    lcmm_data$biomarker <- lcmm_data[,as.character(substitute(name_of_biomarker))]
    
    
    fit.b <- lcmm(biomarker ~ adjusted_new_time + age + age*adjusted_new_time + PTGENDER + PTEDUCAT + apoe, #link = c("5-quantsplines"), 
                  random = ~adjusted_new_time, subject="RID", maxiter = 300, data = lcmm_data[boot,], link = "3-equi-splines")
    
    # boot_predict_splines <- predictlink(fit.b)
    
    boot_pred[i] <- predictY(fit.b, new_data, var.time = "adjusted_new_time", draws = TRUE)
    
    boot_pred_fit <- as.data.frame(boot_pred[i])
    
    boot_derivs[,i] <- diff(boot_pred_fit$Ypred_50)/diff(new_data$adjusted_new_time)
    
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
                                  
################################################################################################  Harmonization Functions  #####################################################################################################
library(ggplot2)
library(mgcv)
library(dplyr)
library(reshape2)
library(devtools)
library(matrixStats)
# devtools::install_github("dgrtwo/ebbr")
library(ebbr)
# BiocManager::install("BiocParallel")
library(BiocParallel)
library(furniture)

################################################################################################    #####################################################################################################
FindIndices <- function(data, trainvar, trainlevel) {
  #
  # Constructs vector of rownames or indices based on subsetting parameters
  # Helpful for GroupModel if interested in training models based
  #  on subset of data
  #
  # Args:
  #  data: full dataframe
  #  trainvar: string of variable name (must be factor variable) 
  #  trainlevel: string of factor level of interest
  #
  # Returns:
  #  Vector of rownames or indices
  #
  index      <- which(data[[trainvar]] == trainlevel, 
                      arr.ind = TRUE)
  index      <- as.vector(index)
  return(index)
}

################################################################################################    #####################################################################################################
BuildDict <- function(covar.data, ref.batch=NULL) {
  #
  # Extract and store batch characteristics for later use
  # Returns dictionary of batch characteristics
  #
  # Args:
  #  covar.data: dataframe of covariates/batch assosciation 
  #
  # Returns:
  #  batches: row indicies split by batch assosciation 
  #  n.batch: number of batches
  #  n.array: number of observations
  #  n.batches: number of observations per batch
  #
  dict    <- list()
  
  batches <- lapply(levels(covar.data[["STUDY"]]), 
                    function(x)which(covar.data[["STUDY"]] == x))
  
  dict[["batches"]]     <- batches 
  dict[["n.batch"]]     <- nlevels(covar.data[["STUDY"]])
  dict[["n.array"]]     <- nrow(covar.data)
  dict[["n.batches"]]   <- sapply(batches, length)
  return(dict)
}

################################################################################################    #####################################################################################################

BuildFormula <- function(covar.data, smooth.terms = NULL, k.val = NULL) {
  #
  # Constructs formula for GAM model
  # Allows specification of smooth predictors along with smoothing parameters
  #
  # Args:
  #  covar.data: Dataframe of covariate data
  #  smooth.terms: Vector of column names for smoothing
  #  k.val:  vector of knot values for smooth column
  # (note) length of k.val must equal length of smooth.terms
  #
  # Returns:
  #  A string to be coerced into formula class for GAM regression
  #
  formstring <- " ~  -1 + STUDY"
  cols <- colnames(covar.data)
  cols <- cols[! cols == "STUDY" ]
  if(!is.null(smooth.terms)) {
    for(i in 1:length(smooth.terms)) {
      smth         <- paste("s", "(",  smooth.terms[i], ",", "k=", 
                            toString(k.val[i]),  ")", sep = "" )
      formstring   <- paste(formstring, smth, sep  = " + ")
    }
    
    cols <- cols[!cols %in% smooth.terms]
  }
  for(i in 1:length(cols)) {
    formstring <- paste(formstring, cols[i], sep= " + " )
  }
  return(formstring)
}

################################################################################################    #####################################################################################################


FitModel <- function(feature.data, covar.data, model.formula, verbose = FALSE) {
  #
  # Fits a GAM model for each feature
  #
  # Args:
  #  feature.data: Imaging data
  #  covar.data: Corresponding covariate data
  #  training.indices: (OPTIONAL) Vector of indices to train models
  #  verbose: (OPTIONAL) Extra output
  #
  # Returns:
  #  A list of GAM models 
  #
  imcols  <- colnames(feature.data)
  covcols <- colnames(covar.data)
  data       <- as.data.frame(cbind(feature.data, covar.data))
  colnames(data) <- c(imcols, covcols)
  feature.data <- as.data.frame(feature.data)
  covar.data <-as.data.frame(covar.data)
  checkzero  <- which(feature.data <= 0, arr.ind = TRUE)[,1]
  
  if(length(checkzero) > 0  ) {
    if(verbose) {
      cat( "[ComGam] ", "Removed", length(checkzero), 
           "rows with non-positive response values from training" , "\n")
    }
    data     <- data[ -checkzero, ]
    
  }
  modlist    <- list()
  for (j in 1:length(imcols)) {
    if(verbose) {
      cat("[ComGam] ", "Fit model", j, "out of", length(imcols), "\n")
    }
    gammod              <- gam(as.formula(paste(imcols[j], model.formula, sep = "")) ,   
                               method = "REML", data = data)
    modlist[[j]]        <- gammod
  }
  names(modlist)        <- imcols
  
  return(modlist)
}

################################################################################################    #####################################################################################################

StanAcrossFeatures <- function(feature.data, covar.data, models.list, data.dict) {
  #
  # Takes data to be harmonized and standardizes the values
  # Removes and preserves biological variance to isolate just batch variance
  # 
  # Args:
  #  dat: Just features to be harmonized (e.x Imaging data/Biomarker data)
  #  covar.data: Dataframe of covariate data
  #  models.list: List of fitted GAM models
  #  data.dict: Dictionary with batch characteristics
  #
  # Returns:
  #  std.data: standardized features data
  #  design:   design matrix
  #  mod.names: names of harmonization features
  #  stand.mean: mean of each feature weighted by batches
  #  var.pooled: variance of each feature
  #
  n.batch   <- data.dict[["n.batch"]]
  n.batches <- data.dict[["n.batches"]]
  n.array   <- data.dict[["n.array"]]
  feature.data       <- t(feature.data)
  design    <- predict.gam(models.list[[1]], 
                           type = "lpmatrix", 
                           newdata = covar.data)
  
  ## dat feature names
  mod.names  <- names(models.list)
  
  ## covariate names
  coef.names <- names(models.list[[1]][["coefficients"]])
  
  B.hat <- data.frame() 
  for(i in 1:length(models.list)) {
    #extract coeffs
    B.hat.row <- models.list[[i]][["coefficients"]]
    B.hat     <- rbind(B.hat, B.hat.row)
  }
  ## build coeff matrix
  colnames(B.hat) <- coef.names
  rownames(B.hat) <- mod.names
  B.hat           <- as.matrix(B.hat)
  B.hat           <- t(B.hat)
  predicted       <- design %*% B.hat
  predicted       <- t(predicted)
  ## find mean value for features using weighted average
  
  grand.mean <- crossprod(n.batches / n.array, B.hat[1 : n.batch, ])
  stand.mean <- crossprod(grand.mean, t(rep(1, n.array)))
  stand.mean1 <- stand.mean
  
  # add bio variance to stand.mean
  design.tmp <- design
  design.tmp[ ,c(1 : n.batch)] <- 0
  bio.var    <- design.tmp %*% B.hat
  bio.var    <- t(bio.var)
  stand.mean <- stand.mean + bio.var
  
  ## variance
  var.pooled  <- (feature.data - predicted) ^ 2
  var.adj     <- as.matrix(rep(1, n.array)) / n.array
  var.pooled  <- var.pooled %*% var.adj
  var.pooled  <- sqrt(var.pooled)
  var.pooled1 <- var.pooled
  var.pooled  <- var.pooled %*% as.matrix(t(rep(1, n.array)))
  std.data    <- (feature.data - stand.mean) / var.pooled
  return(list("std.data"   = std.data, 
              "design"     = design, 
              "mod.names"  = mod.names,
              "stand.mean" = stand.mean,
              "var.pooled" = var.pooled,
              "grand.mean" = grand.mean,
              "var.pooled1" = var.pooled1))
}



################################################################################################    #####################################################################################################

CalcGammaDelta <- function(stan.dict, data.dict) {
  #
  # Calculates the Mean/Variance (Shift/Scale) of batch effect
  # 
  # Args:
  # stan.dict: Dictionary with standardized data and other values
  # Output of StanAcrossFeatures
  #
  # Returns:
  # gamma.hat: Shift value per batch per feature
  # delta.hat: Scale value per batch per feature
  #
  batches   <- data.dict[["batches"]]
  n.batch   <- data.dict[["n.batch"]]
  design    <- stan.dict[["design"]]
  std.data  <- stan.dict[["std.data"]]
  mod.names <- stan.dict[["mod.names"]]
  batch.mod <- design[, 1:n.batch]
  
  ## shift (mean)
  gamma.hat <- tcrossprod(solve(crossprod(batch.mod, batch.mod)), batch.mod)
  gamma.hat <- tcrossprod(gamma.hat, std.data)
  
  ## scale (variance)
  delta.hat <- data.frame()
  for(i in batches) {
    delta.hat <- rbind(delta.hat, 
                       rowVars(std.data, cols = i, na.rm = TRUE))
  }
  colnames(delta.hat) <- mod.names
  rownames(delta.hat) <- colnames(batch.mod)
  return(list("gamma.hat" = gamma.hat,
              "delta.hat" = delta.hat,
              "stan.data" = std.data))
}


################################################################################################    #####################################################################################################

ModelDiagnostics <- function(mod.list) {
  #
  # Pulls information on model fitting for each features smooth terms 
  # Useful for checking model fit accuracy in harmonizationn
  #
  # Args:
  #  mod.list: Output from FitModel
  #
  # Returns:
  #  matrix with diagnostic information 
  #
  diag.list  <-list()
  mod.names  <- names(mod.list)
  diagnos.df <- data.frame(matrix(ncol = 7))
  for (i in 1:length(mod.list)) {
    gam.o <- capture.output(gam.check(mod.list[[i]]))
    dev.off()
    for(j in 13:length(gam.o)) {
      o.vals.split <- c(strsplit(gam.o[j], split = " "))
      o.vals <- o.vals.split[[1]]
      o.vals <- o.vals[o.vals != ""]
      if("---" %in% o.vals | "Signif." %in% o.vals) {
        #skip
      } else {
        if(length(o.vals) == 6) {
          row.val <- c(mod.names[i], o.vals)
        } else {
          row.val <- c(mod.names[i], o.vals, "")
        }
      }
      
      diagnos.df <- rbind(diagnos.df, as.character(row.val))
    }
  }
  diagnos.df <- as.data.frame(diagnos.df)
  colnames(diagnos.df) <- c("feature",
                            "smooth",
                            "k",
                            "edf",
                            "k.index",
                            "p.value",
                            "signif code (alpha = 0.05")
  diagnos.df <- unique(diagnos.df[2:nrow(diagnos.df), ])
  return(diagnos.df)
}
################################################################################################    #####################################################################################################

ApplyGammaDelta <- function(stan.dict, site.params, data.dict) {
  #
  # Takes gamma and delta values and applies them to the data
  # Reintroduces biological variance after gamma/delta correction
  #
  # Args:
  #  stan.dict: StanAcrossFeatures output
  #  site.params: CalcGammaDelta output
  #  data.dict: BuildDict Output
  #
  # Returns:
  #  Harmonized data
  #
  ## extract all variables
  std.data   <- stan.dict[["std.data"]]
  design     <- stan.dict[["design"]]
  mod.names  <- stan.dict[["mod.names"]]
  stand.mean <- stan.dict[["stand.mean"]]
  var.pooled <- stan.dict[["var.pooled"]]
  batches    <- data.dict[["batches"]]
  n.batches  <- data.dict[["n.batches"]]
  n.batch    <- data.dict[["n.batch"]]
  gamma.hat  <- site.params[["gamma.hat"]]
  delta.hat  <- as.matrix(site.params[["delta.hat"]])
  batch.mod  <- design[, 1:n.batch]
  ## apply gamma/delta
  std.adj <- std.data
  j <- 1 
  for(i in batches){ 
    shft <- std.adj[ ,i] - t(batch.mod[i, ] %*% gamma.hat)
    scl  <- tcrossprod(sqrt(delta.hat[j, ]), rep(1, n.batches[j]))
    std.adj[,i] <- shft / scl
    j <- j + 1
  }
  
  #reintroduce biological mean and variance
  std.adj <- (std.adj * var.pooled) + stand.mean
  return(std.adj)
}

#######################################################################################
getEbEstimators <- function(naiveEstimators,
                            s.data, 
                            dataDict,
                            parametric=TRUE, 
                            mean.only=FALSE,
                            BPPARAM=bpparam("SerialParam")
){
  gamma.hat <- naiveEstimators[["gamma.hat"]]
  delta.hat <- naiveEstimators[["delta.hat"]]
  batches   <- dataDict$batches
  n.batch   <- dataDict$n.batch
  ref.batch <- dataDict$ref.batch
  ref <- dataDict$ref
  
  .getParametricEstimators <- function(){
    gamma.star <- delta.star <- NULL
    for (i in 1:n.batch){
      if (mean.only){
        gamma.star <- rbind(gamma.star, postmean(gamma.hat[i,], gamma.bar[i], 1, 1, t2[i]))
        delta.star <- rbind(delta.star, rep(1, nrow(s.data)))
      } else {
        temp <- it.sol(s.data[,batches[[i]]],
                       gamma.hat[i,],
                       delta.hat[i,],
                       gamma.bar[i],
                       t2[i],
                       a.prior[i],
                       b.prior[i])
        gamma.star <- rbind(gamma.star,temp[1,])
        delta.star <- rbind(delta.star,temp[2,])
      }
    }
    rownames(gamma.star) <- rownames(delta.star) <- names(batches)
    out <- list(gamma.star=gamma.star, delta.star=delta.star)
    return(out)
  }
  .getNonParametricEstimators <- function(BPPARAM=bpparam("SerialParam")){
    gamma.star <- delta.star <- NULL
    
    results <- bplapply(1:n.batch, function(i){
      if (mean.only){
        delta.hat[i, ] = 1
      }
      temp <- int.eprior(as.matrix(s.data[, batches[[i]]]),
                         gamma.hat[i,],
                         delta.hat[i,])
      return(temp)
    }, BPPARAM = BPPARAM)
    gamma.star <- lapply(results, function(x) x[1,])
    delta.star <- lapply(results, function(x) x[2,])
    gamma.star <- do.call("rbind",gamma.star)
    delta.star <- do.call("rbind",delta.star)
    #for (i in 1:n.batch){
    #    gamma.star <- rbind(gamma.star,temp[1,])
    #    delta.star <- rbind(delta.star,temp[2,])
    #}
    rownames(gamma.star) <- rownames(delta.star) <- names(batches)
    out <- list(gamma.star=gamma.star, delta.star=delta.star)
    return(out)
  }
  
  gamma.bar <- rowMeans(gamma.hat, na.rm=TRUE)
  t2 <- rowVars(gamma.hat, na.rm=TRUE)
  names(t2) <- rownames(gamma.hat)
  a.prior <- apriorMat(delta.hat)
  b.prior <- bpriorMat(delta.hat)
  if (parametric){
    temp <- .getParametricEstimators()
  } else {
    temp <- .getNonParametricEstimators(BPPARAM=BPPARAM)
  }
  if(!is.null(ref.batch)){
    temp[["gamma.star"]][ref,] <- 0  ## set reference batch mean equal to 0
    temp[["delta.star"]][ref,] <- 1  ## set reference batch variance equal to 1
  }
  out <- list()
  out[["gamma.star"]] <- temp[["gamma.star"]]
  out[["delta.star"]] <- temp[["delta.star"]]
  out[["gamma.bar"]] <- gamma.bar
  out[["t2"]] <- t2
  out[["a.prior"]] <- a.prior
  out[["b.prior"]] <- b.prior
  return(out)
}

################################################################################################    #####################################################################################################

################################################################################################    #####################################################################################################

aprior <- function(delta.hat){
  m=mean(delta.hat)
  s2=var(delta.hat)
  return((2*s2+m^2)/s2)
}

bprior <- function(delta.hat){
  m=mean(delta.hat)
  s2=var(delta.hat)
  return((m*s2+m^3)/s2)
}

apriorMat <- function(delta.hat) {
  m  <- rowMeans2(delta.hat)
  s2 <- rowVars(delta.hat)
  out <- (2*s2+m^2)/s2
  names(out) <- rownames(delta.hat)
  return(out)
}

bpriorMat <- function(delta.hat) {
  m <- rowMeans2(delta.hat)
  s2 <- rowVars(delta.hat)
  out <- (m*s2+m^3)/s2
  names(out) <- rownames(delta.hat)
  return(out)
}

postmean <- function(g.hat, g.bar, n, d.star, t2){
  (t2*n*g.hat+d.star*g.bar)/(t2*n+d.star)
}

postvar <- function(sum2, n, a, b){
  (.5*sum2+b)/(n/2+a-1)
}

# Helper function for parametric adjustements:
it.sol  <- function(sdat, g.hat, d.hat, g.bar, t2, a, b, conv=.0001){
  #n <- apply(!is.na(sdat),1,sum)
  n <- rowSums(!is.na(sdat))
  g.old  <- g.hat
  d.old  <- d.hat
  change <- 1
  count  <- 0
  ones <- rep(1,ncol(sdat))
  
  while(change>conv){
    g.new  <- postmean(g.hat,g.bar,n,d.old,t2)
    sum2   <- rowSums2((sdat-tcrossprod(g.new, ones))^2, na.rm=TRUE)
    d.new  <- postvar(sum2,n,a,b)
    change <- max(abs(g.new-g.old)/g.old,abs(d.new-d.old)/d.old)
    g.old <- g.new
    d.old <- d.new
    count <- count+1
  }
  adjust <- rbind(g.new, d.new)
  rownames(adjust) <- c("g.star","d.star")
  return(adjust)
}



# Helper function for non-parametric adjustements:
int.eprior <- function(sdat, g.hat, d.hat){
  g.star <- d.star <- NULL
  r <- nrow(sdat)
  for(i in 1:r){
    g <- g.hat[-i]
    d <- d.hat[-i]		
    x <- sdat[i,!is.na(sdat[i,])]
    n <- length(x)
    j <- numeric(n)+1
    dat <- matrix(as.numeric(x), length(g), n, byrow=TRUE)
    resid2 <- (dat-g)^2
    sum2 <- resid2 %*% j
    LH <- 1/(2*pi*d)^(n/2)*exp(-sum2/(2*d))
    LH[LH=="NaN"]=0
    g.star <- c(g.star, sum(g*LH)/sum(LH))
    d.star <- c(d.star, sum(d*LH)/sum(LH))
  }
  adjust <- rbind(g.star,d.star)
  rownames(adjust) <- c("g.star","d.star")
  return(adjust)
} 

###########################################################################

###########################################################################
ApplyHarm <- function(feature.data, covariate.data, comgam.out) {
  std.data           <- comgam.out$stan.dict
  site.params        <- comgam.out$shift.scale.params
  models.list        <- comgam.out$models.list
  datadict           <- BuildDict(covariate.data)
  batches            <- datadict[["batches"]]
  n.batches          <- datadict[["n.batches"]]
  n.batch            <- datadict[["n.batch"]]
  n.array            <- datadict[["n.array"]]
  grand.mean         <- std.data[["grand.mean"]]
  var.pooled1        <- std.data[["var.pooled1"]]
  feature.data       <- t(feature.data)
  design             <- predict.gam(models.list[[1]], 
                                    type = "lpmatrix", 
                                    newdata = covariate.data)
  
  mod.names  <- names(models.list)
  stand.mean <- crossprod(grand.mean, t(rep(1, n.array)))
  var.pooled <- var.pooled1 %*% as.matrix(t(rep(1, n.array)))
  ## covariate names
  coef.names <- names(models.list[[1]][["coefficients"]])
  B.hat <- data.frame() 
  for(i in 1:length(models.list)) {
    #extract coeffs
    B.hat.row <- models.list[[i]][["coefficients"]]
    B.hat     <- rbind(B.hat, B.hat.row)
  }
  ## build coeff matrix
  colnames(B.hat) <- coef.names
  rownames(B.hat) <- mod.names
  B.hat           <- as.matrix(B.hat)
  B.hat           <- t(B.hat)
  predicted       <- design %*% B.hat
  predicted       <- t(predicted)
  design.tmp      <- design
  design.tmp[ ,c(1 : n.batch)] <- 0
  bio.var    <- design.tmp %*% B.hat
  bio.var    <- t(bio.var)
  stand.mean <- stand.mean + bio.var
  std.data   <- (feature.data - stand.mean) / var.pooled
  gamma.hat  <- site.params[["gamma.hat"]]
  delta.hat  <- as.matrix(site.params[["delta.hat"]])
  batch.mod  <- design[, 1:n.batch]
  ## apply gamma/delta
  std.adj <- std.data
  j <- 1 
  for(i in batches) { 
    shft <- std.adj[ ,i] - t(batch.mod[i, ] %*% gamma.hat)
    scl  <- tcrossprod(sqrt(delta.hat[j, ]), rep(1, n.batches[j]))
    std.adj[,i] <- shft / scl
    j <- j + 1
  }
  #reintroduce biological mean and variance
  std.adj <- (std.adj * var.pooled) + stand.mean
  return(std.adj)
}


ComGamHarm <- function(feature.data, 
                       covar.data, 
                       eb                = FALSE,
                       parametric        = TRUE,
                       smooth.terms      = NULL,
                       k.val             = NULL,
                       verbose           = TRUE,
                       model.diagnostics = FALSE) {
  
  ### check object types and parameters ###
  if(!is.data.frame(feature.data)) {
    stop('ComGamHarm: feature.data is not of type "data.frame"')
  }
  if(!is.data.frame(covar.data)) {
    stop('ComGamHarm: covar.data is not of type "data.frame"')
  }
  if(!("STUDY" %in% colnames(covar.data))) {
    stop('ComGamHarm: covar.data must include column "STUDY"')
  }
  if(!is.factor(covar.data[["STUDY"]])) {
    stop('ComGamHarm: "STUDY" column must be of class "factor"')
  }
  if(any(is.na(feature.data)) | any(is.na(covar.data))) {
    stop('ComGamHarm: data frames cannot contain missing data')
  }
  if(!(nrow(feature.data) == nrow(covar.data))) {
    stop('ComGamHarm: # of rows inconsistent between data frames')
  }
  if(!is.null(smooth.terms) | !is.null(k.val)) {
    if(!(!is.null(smooth.terms) & !is.null(k.val))) {
      stop('ComGamHarm: both smooth.terms & k.val must be supplied if one is supplied')
    }
    if(!(length(smooth.terms == length(k.val)))) {
      stop('ComGamHarm: smooth.terms & k.val must be vectors of same length')
      
    }
  }
  data.dict        <-  BuildDict(covar.data =  covar.data)
  
  model.formula    <-  BuildFormula(covar.data   = covar.data,
                                    smooth.terms = smooth.terms,
                                    k.val        = k.val)
  
  models.list      <-  FitModel(feature.data     = feature.data,
                                covar.data       = covar.data,
                                model.formula    = model.formula,
                                verbose = verbose)
  
  stan.dict        <-  StanAcrossFeatures(feature.data = feature.data,
                                          covar.data   = covar.data,
                                          models.list  = models.list,
                                          data.dict    = data.dict)
  
  features.adj     <-  CalcGammaDelta(stan.dict =  stan.dict,
                                      data.dict = data.dict)
  
  features.adj.untouched <- features.adj
  
  if(eb) {
    gamma.hat <- as.matrix(features.adj$gamma.hat)
    delta.hat <- as.matrix(features.adj$delta.hat)
    
    naiveest <- list("gamma.hat" = gamma.hat,
                     "delta.hat" = delta.hat)
    stddata   <- stan.dict$std.data
    data.dict <- data.dict
    data.dict[["ref.batch"]] <- NULL
    EbEst <- getEbEstimators(naiveEstimators = naiveest,
                             s.data = stddata,
                             dataDict = data.dict,
                             parametric = parametric)
    
    features.adj <- list("gamma.hat" = EbEst[["gamma.star"]],
                         "delta.hat" = EbEst[["delta.star"]])
    
  }
  
  
  
  features.results <-  ApplyGammaDelta(stan.dict   = stan.dict,
                                       site.params = features.adj,
                                       data.dict   = data.dict)
  if(!eb) {
    priors <- NULL
  } else {
    priors <- EbEst
  }
  
  if(model.diagnostics) {
    
    mod.diags <- ModelDiagnostics(mod.list = models.list)
    
    
    return.list <- list("harm.results"       = features.results,
                        "stan.dict"          = stan.dict,
                        "shift.scale.params" = features.adj,
                        "models.list"        = models.list,
                        "model.formula"      = model.formula,
                        "data.dict"          = data.dict,
                        "model.diagnostics"  = mod.diags,
                        "priors"             = priors,
                        "priors_no_eb"       = features.adj.untouched,
                        "feature_data"       = feature.data,
                        "covs_data"          = covar.data)
  } else {
    
    return.list <- list("harm.results"       = features.results,
                        "stan.dict"          = stan.dict,
                        "shift.scale.params" = features.adj,
                        "models.list"        = models.list,
                        "model.formula"      = model.formula,
                        "data.dict"          = data.dict,
                        "priors"             = priors,
                        "priors_no_eb"       = features.adj.untouched,
                        "feature_data"       = feature.data,
                        "covs_data"          = covar.data) 
  }
  return(return.list)
}
