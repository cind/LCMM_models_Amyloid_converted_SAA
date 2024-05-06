library(mgcv)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(reshape2)
library(matrixStats)
library(ebbr)
library(BiocParallel)
library(patchwork)

source("helper_functions_parallel.R")
root = "lcmm_modelling_data"

######################################################################################################################################
centiloid_plot_data <- read.csv(file.path(root, "centiloid_lcmm_data.csv"))
tau_plot_data <- read.csv(file.path(root, "tau_lcmm_data.csv"))
ptau_plot_data <- read.csv(file.path(root, "ptau_lcmm_data.csv"))
abeta_plot_data <- read.csv(file.path(root, "abeta_lcmm_data.csv"))
npi_plot_data <- read.csv(file.path(root, "npi_lcmm_data.csv"))
mpacctrailsb_plot_data <- read.csv(file.path(root, "mpacctrailsb_lcmm_data.csv"))
cdrsb_plot_data <- read.csv(file.path(root, "cdrsb_lcmm_data.csv"))
mmse_plot_data <- read.csv(file.path(root, "mmse_lcmm_data.csv"))
fdg_plot_data <- read.csv(file.path(root, "fdg_lcmm_data.csv"))
meta_roi_plot_data <- read.csv(file.path(root, "meta_roi_lcmm_data_no_adni1.csv"))
hippocampal_volume_plot_data <- read.csv(file.path(root, "hippocampal_volume_lcmm_data_no_adni1.csv"))
adas13_plot_data <- read.csv(file.path(root, "adas13_lcmm_data.csv"))
ecog_s_plot_data <- read.csv(file.path(root, "ecog_s_lcmm_data.csv"))
ecog_p_plot_data <- read.csv(file.path(root, "ecog_p_lcmm_data.csv"))

test <- data.frame(adjusted_new_time = seq(-10, 10, by = 0.1), 
                   RID = rep(0, length(seq(-10, 10, by = 0.1))), 
                   PTGENDER = rep(2, length(seq(-10, 10, by = 0.1))),
                   PTEDUCAT = rep(16, length(seq(-10, 10, by = 0.1))),
                   age = rep(73, length(seq(-10, 10, by = 0.1))),
                   apoe = "E3", 
                   CDGLOBAL = rep(0, length(seq(-10, 10, by = 0.1))))

# can also use formula: mgcv::gam(outcome_var ~ s(adjusted_new_time, bs = "cs", k = 4, sp = 2) +
#                                   s(RID, bs = "re") +
#                                   age + PTEDUCAT + apoe + PTGENDER,
#                                 data = var_plot_data,
#                                 family = gaussian,
#                                 method = "REML")

#####################################################################################################
# Centiloid
#####################################################################################################

centiloid_plot_data$RID <- as.numeric(centiloid_plot_data$RID)

gam.model_centiloid <- mgcv::gam(formula = Centiloid ~ s(adjusted_new_time, bs = "cs", k = 4, fx = TRUE) +
                                   s(RID, bs = "re") + 
                                   age + PTEDUCAT + apoe + PTGENDER + CDGLOBAL,
                                 data = centiloid_plot_data,
                                 family = gaussian,
                                 method = "ML")

m2.d_centiloid <- Deriv(gam.model_centiloid, newdata = test) #Calculating derivatives of the model

m2.dci_centiloid <- confint(m2.d_centiloid, term = 'adjusted_new_time') #Calculating upper and lower confidence interval

m2.dsig_centiloid <- signifD(seq(-10, 10, by = 0.1), d = m2.d_centiloid[['adjusted_new_time']]$deriv, m2.dci_centiloid[['adjusted_new_time']]$upper, m2.dci_centiloid[['adjusted_new_time']]$lower) # calculating the significance of change

modelFit_centiloid <- data.frame(predict(gam.model_centiloid, se = TRUE,
                                         data.frame(adjusted_new_time = centiloid_plot_data$adjusted_new_time, 
                                                    RID = 0, 
                                                    age = 73, 
                                                    PTEDUCAT = 16, 
                                                    PTGENDER = 2,
                                                    CDGLOBAL = 0,
                                                    apoe = "E3")))

df.fit_centiloid <- data.frame(cbind(adjusted_new_time = centiloid_plot_data$adjusted_new_time,
                                     meas = as.numeric(centiloid_plot_data$Centiloid),
                                     fit = modelFit_centiloid$fit,
                                     upperBound = modelFit_centiloid$fit + 2 * modelFit_centiloid$se.fit,
                                     lowerBound = modelFit_centiloid$fit - 2 * modelFit_centiloid$se.fit))

df.fit_centiloid <- df.fit_centiloid[!is.na(df.fit_centiloid$meas),]

meas_fn_centiloid <- "PET Aβ Centiloid"

gam.p_centiloid <- ggplot(data=df.fit_centiloid, aes(adjusted_new_time,meas)) + ylim(c(min(df.fit_centiloid$lowerBound),max(df.fit_centiloid$upperBound))) +
  geom_line(aes(adjusted_new_time, fit), color = "black", size=0.25, data=df.fit_centiloid) +
  ggplot2::geom_line(aes(y=upperBound),se=FALSE, stat="smooth",linetype="dashed",alpha=.5) +
  ggplot2::geom_line(aes(y=lowerBound),se=FALSE, stat="smooth",linetype="dashed",alpha=.5) +
  theme_blank() + theme(axis.line = element_line(size = 0.5, colour = "black", linetype=1)) + 
  ggtitle("") + 
  theme(plot.title=element_text(hjust=0.5),text=element_text(size=15)) +
  theme(legend.position="none") +
  xlab("Years from Aβ+") + ylab(meas_fn_centiloid)

inter_centiloid <- na.contiguous(m2.dsig_centiloid$incr)

gam.p_centiloid <- gam.p_centiloid + 
  geom_line(aes(adjusted_new_time, fit), color = "black", size=1.5, data=df.fit_centiloid[df.fit_centiloid$adjusted_new_time>=min(inter_centiloid) & df.fit_centiloid$adjusted_new_time<=max(inter_centiloid),]) + 
  ggplot2::geom_vline(aes(xintercept=min(inter_centiloid)),linetype="dashed",color="black") + #inflection point
  ggplot2::annotate(geom="text",label=paste0(round(min(inter_centiloid),digits=1)),x=(min(inter_centiloid)-1),angle='90',y=10)

#####################################################################################################
# fdg
#####################################################################################################

fdg_plot_data$RID <- as.numeric(fdg_plot_data$RID)

gam.model_fdg <- mgcv::gam(formula = adjusted_Meta_ROI ~ s(adjusted_new_time, bs = "cs", k = 4, fx = TRUE) +
                             s(RID, bs = "re") + 
                             age + PTEDUCAT + apoe + PTGENDER + CDGLOBAL,
                           data = fdg_plot_data,
                           family = gaussian,
                           method = "ML")

m2.d_fdg <- Deriv(gam.model_fdg, newdata = test) #Calculating derivatives of the model


m2.dci_fdg <- confint(m2.d_fdg, term = 'adjusted_new_time') #Calculating upper and lower confidence interval

m2.dsig_fdg <- signifD(seq(-10,10, by = 0.1), d = m2.d_fdg[['adjusted_new_time']]$deriv, m2.dci_fdg[['adjusted_new_time']]$upper, m2.dci_fdg[['adjusted_new_time']]$lower) # calculating the significance of change

modelFit_fdg <- data.frame(predict(gam.model_fdg, se = TRUE,
                                   data.frame(adjusted_new_time = fdg_plot_data$adjusted_new_time, 
                                              RID = 0, 
                                              age = 73, 
                                              PTEDUCAT = 16,
                                              CDGLOBAL = 0,
                                              PTGENDER = 2,
                                              apoe = "E3")))

df.fit_fdg <- data.frame(cbind(adjusted_new_time = fdg_plot_data$adjusted_new_time,
                               meas = fdg_plot_data$adjusted_Meta_ROI,
                               fit = modelFit_fdg$fit,
                               upperBound = modelFit_fdg$fit + 2 * modelFit_fdg$se.fit,
                               lowerBound = modelFit_fdg$fit - 2 * modelFit_fdg$se.fit))

df.fit_fdg <- df.fit_fdg[!is.na(df.fit_fdg$meas),]

meas_fn_fdg <- "FDG PET Meta-ROI"

gam.p_fdg <- ggplot(data=df.fit_fdg, aes(adjusted_new_time,meas)) + ylim(c(min(df.fit_fdg$lowerBound),max(df.fit_fdg$upperBound))) +
  geom_line(aes(adjusted_new_time, fit), color = "black", size=0.25, data=df.fit_fdg) +
  ggplot2::geom_line(aes(y=upperBound),se=FALSE, stat="smooth",linetype="dashed",alpha=.5) +
  ggplot2::geom_line(aes(y=lowerBound),se=FALSE, stat="smooth",linetype="dashed",alpha=.5) +
  theme_blank() + theme(axis.line = element_line(size = 0.5, colour = "black", linetype=1)) + 
  ggtitle("") + 
  theme(plot.title=element_text(hjust=0.5),text=element_text(size=15)) +
  theme(legend.position="none") +
  xlab("Years from Aβ+") + ylab(meas_fn_fdg)

inter_fdg <- na.contiguous(m2.dsig_fdg$decr)

gam.p_fdg <- gam.p_fdg + 
  geom_line(aes(adjusted_new_time, fit), color = "black", size=1.5, data=df.fit_fdg[df.fit_fdg$adjusted_new_time>=min(inter_fdg) & df.fit_fdg$adjusted_new_time<=max(inter_fdg),]) + 
  ggplot2::geom_vline(aes(xintercept=min(inter_fdg)),linetype="dashed",color="black") + #inflection point
  ggplot2::annotate(geom="text",label=paste0(round(min(inter_fdg),digits=1)),x=(min(inter_fdg)-1),angle='90',y=1.28)

#####################################################################################################
# csf tau
#####################################################################################################

tau_plot_data$RID <- as.numeric(tau_plot_data$RID)

gam.model_tau <- mgcv::gam(formula = TAU ~ s(adjusted_new_time, bs = "cs", k = 4, fx = TRUE) + # k = 4????, fx = TRUE, bs = "cs"
                             s(RID, bs = "re") + 
                             age + PTEDUCAT + apoe + PTGENDER + CDGLOBAL,
                           data = tau_plot_data,
                           family = gaussian,
                           method = "ML")

m2.d_tau <- Deriv(gam.model_tau, newdata = test) #Calculating derivatives of the model

m2.dci_tau <- confint(m2.d_tau, term = 'adjusted_new_time') #Calculating upper and lower confidence interval

m2.dsig_tau <- signifD(seq(-10,10, by = 0.1), d = m2.d_tau[['adjusted_new_time']]$deriv, m2.dci_tau[['adjusted_new_time']]$upper, m2.dci_tau[['adjusted_new_time']]$lower) # calculating the significance of change

modelFit_tau <- data.frame(predict(gam.model_tau, se = TRUE,
                                   data.frame(adjusted_new_time = tau_plot_data$adjusted_new_time, 
                                              RID = 0, 
                                              age = 73, 
                                              PTEDUCAT = 16, 
                                              CDGLOBAL = 0,
                                              PTGENDER = 2,
                                              apoe = "E3"))) 

df.fit_tau <- data.frame(cbind(adjusted_new_time = tau_plot_data$adjusted_new_time,
                               meas = tau_plot_data$TAU,
                               fit = modelFit_tau$fit,
                               upperBound = modelFit_tau$fit + 2 * modelFit_tau$se.fit,
                               lowerBound = modelFit_tau$fit - 2 * modelFit_tau$se.fit))

df.fit_tau <- df.fit_tau[!is.na(df.fit_tau$meas),]

meas_fn_tau <- "CSF Total Tau"

gam.p_tau <- ggplot(data=df.fit_tau, aes(adjusted_new_time,meas)) + ylim(c(min(df.fit_tau$lowerBound),max(df.fit_tau$upperBound))) +
  geom_line(aes(adjusted_new_time, fit), color = "black", size=0.25, data=df.fit_tau) +
  ggplot2::geom_line(aes(y=upperBound),se=FALSE, stat="smooth",linetype="dashed",alpha=.5) +
  ggplot2::geom_line(aes(y=lowerBound),se=FALSE, stat="smooth",linetype="dashed",alpha=.5) +
  theme_blank() + theme(axis.line = element_line(size = 0.5, colour = "black", linetype=1)) + 
  ggtitle("") + 
  theme(plot.title=element_text(hjust=0.5),text=element_text(size=15)) +
  theme(legend.position="none") +
  xlab("Years from Aβ+") + ylab(meas_fn_tau)

inter_tau <- na.contiguous(m2.dsig_tau$incr)

gam.p_tau <- gam.p_tau + 
  geom_line(aes(adjusted_new_time, fit), color = "black", size=1.5, data=df.fit_tau[df.fit_tau$adjusted_new_time>=min(inter_tau) & df.fit_tau$adjusted_new_time<=max(inter_tau),]) + 
  ggplot2::geom_vline(aes(xintercept=min(inter_tau)),linetype="dashed",color="black") + #inflection point
  ggplot2::annotate(geom="text",label=paste0(round(min(inter_tau),digits=1)),x=(min(inter_tau)-1),angle='90',y=300)

#####################################################################################################
# csf ptau
#####################################################################################################

ptau_plot_data$RID <- as.numeric(ptau_plot_data$RID)

gam.model_ptau <- mgcv::gam(formula = PTAU ~ s(adjusted_new_time, bs = "cs", k = 4, fx = TRUE) +
                              s(RID, bs = "re") + 
                              age + PTEDUCAT + apoe + PTGENDER + CDGLOBAL,
                            data = ptau_plot_data,
                            family = gaussian,
                            method = "ML")

m2.d_ptau <- Deriv(gam.model_ptau, newdata = test) #Calculating derivatives of the model

m2.dci_ptau <- confint(m2.d_ptau, term = 'adjusted_new_time') #Calculating upper and lower confidence interval

m2.dsig_ptau <- signifD(seq(-10,10, by = 0.1), d = m2.d_ptau[['adjusted_new_time']]$deriv, m2.dci_ptau[['adjusted_new_time']]$upper, m2.dci_ptau[['adjusted_new_time']]$lower) # calculating the significance of change

modelFit_ptau <- data.frame(predict(gam.model_ptau, se = TRUE,
                                    data.frame(adjusted_new_time = ptau_plot_data$adjusted_new_time, 
                                               RID = 0, 
                                               age = 73, 
                                               PTEDUCAT = 16, 
                                               CDGLOBAL = 0,
                                               PTGENDER = 2,
                                               apoe = "E3"))) 

df.fit_ptau <- data.frame(cbind(adjusted_new_time = ptau_plot_data$adjusted_new_time,
                                meas = ptau_plot_data$PTAU,
                                fit = modelFit_ptau$fit,
                                upperBound = modelFit_ptau$fit + 2 * modelFit_ptau$se.fit,
                                lowerBound = modelFit_ptau$fit - 2 * modelFit_ptau$se.fit))

df.fit_ptau <- df.fit_ptau[!is.na(df.fit_ptau$meas),]

meas_fn_ptau <- "CSF P-Tau181"

gam.p_ptau <- ggplot(data=df.fit_ptau, aes(adjusted_new_time,meas)) + ylim(c(min(df.fit_ptau$lowerBound),max(df.fit_ptau$upperBound))) +
  geom_line(aes(adjusted_new_time, fit), color = "black", size=0.25, data=df.fit_ptau) +
  ggplot2::geom_line(aes(y=upperBound),se=FALSE, stat="smooth",linetype="dashed",alpha=.5) +
  ggplot2::geom_line(aes(y=lowerBound),se=FALSE, stat="smooth",linetype="dashed",alpha=.5) +
  theme_blank() + theme(axis.line = element_line(size = 0.5, colour = "black", linetype=1)) + 
  ggtitle("") + 
  theme(plot.title=element_text(hjust=0.5),text=element_text(size=15)) +
  theme(legend.position="none") +
  xlab("Years from Aβ+") + ylab(meas_fn_ptau)

inter_ptau <- na.contiguous(m2.dsig_ptau$incr)

gam.p_ptau <- gam.p_ptau + 
  geom_line(aes(adjusted_new_time, fit), color = "black", size=1.5, data=df.fit_ptau[df.fit_ptau$adjusted_new_time>=min(inter_ptau) & df.fit_ptau$adjusted_new_time<=max(inter_ptau),]) + 
  ggplot2::geom_vline(aes(xintercept=min(inter_ptau)),linetype="dashed",color="black") + #inflection point
  ggplot2::annotate(geom="text",label=paste0(round(min(inter_ptau),digits=1)),x=(min(inter_ptau)-1),angle='90',y=27)

#####################################################################################################
# csf Amyloid
#####################################################################################################

abeta_plot_data$RID <- as.numeric(abeta_plot_data$RID)

gam.model_abeta <- mgcv::gam(formula = ABETA ~ s(adjusted_new_time, bs = "cs", k = 4, fx = TRUE) +
                               s(RID, bs = "re") + 
                               age + PTEDUCAT + apoe + PTGENDER + CDGLOBAL,
                             data = abeta_plot_data,
                             family = gaussian,
                             method = "ML")

m2.d_abeta <- Deriv(gam.model_abeta, newdata = test) #Calculating derivatives of the model

m2.dci_abeta <- confint(m2.d_abeta, term = 'adjusted_new_time') #Calculating upper and lower confidence interval

m2.dsig_abeta <- signifD(seq(-10,10, by = 0.1), d = m2.d_abeta[['adjusted_new_time']]$deriv, m2.dci_abeta[['adjusted_new_time']]$upper, m2.dci_abeta[['adjusted_new_time']]$lower) # calculating the significance of change

modelFit_abeta <- data.frame(predict(gam.model_abeta, se = TRUE,
                                     data.frame(adjusted_new_time = abeta_plot_data$adjusted_new_time, 
                                                RID = 0, 
                                                age = 73, 
                                                PTEDUCAT = 16, 
                                                CDGLOBAL = 0,
                                                PTGENDER = 2,
                                                apoe = "E3")))

df.fit_abeta <- data.frame(cbind(adjusted_new_time = abeta_plot_data$adjusted_new_time,
                                 meas = abeta_plot_data$ABETA,
                                 fit = modelFit_abeta$fit,
                                 upperBound = modelFit_abeta$fit + 2 * modelFit_abeta$se.fit,
                                 lowerBound = modelFit_abeta$fit - 2 * modelFit_abeta$se.fit))

df.fit_abeta <- df.fit_abeta[!is.na(df.fit_abeta$meas),]

meas_fn_abeta <- "CSF Aβ42"

gam.p_abeta <- ggplot(data=df.fit_abeta, aes(adjusted_new_time,meas)) + ylim(c(min(df.fit_abeta$lowerBound),max(df.fit_abeta$upperBound))) +
  geom_line(aes(adjusted_new_time, fit), color = "black", size=0.25, data=df.fit_abeta) +
  ggplot2::geom_line(aes(y=upperBound),se=FALSE, stat="smooth",linetype="dashed",alpha=.5) +
  ggplot2::geom_line(aes(y=lowerBound),se=FALSE, stat="smooth",linetype="dashed",alpha=.5) +
  theme_blank() + theme(axis.line = element_line(size = 0.5, colour = "black", linetype=1)) + 
  ggtitle("") + 
  theme(plot.title=element_text(hjust=0.5),text=element_text(size=15)) +
  theme(legend.position="none") +
  xlab("Years from Aβ+") + ylab(meas_fn_abeta)

inter_abeta <- na.contiguous(m2.dsig_abeta$decr)

gam.p_abeta <- gam.p_abeta + 
  geom_line(aes(adjusted_new_time, fit), color = "black", size=1.5, data=df.fit_abeta[df.fit_abeta$adjusted_new_time>=min(inter_abeta) & df.fit_abeta$adjusted_new_time<=max(inter_abeta),]) + 
  ggplot2::geom_vline(aes(xintercept=min(inter_abeta)),linetype="dashed",color="black") + #inflection point
  ggplot2::annotate(geom="text",label=paste0(round(min(inter_abeta),digits=1)),x=(min(inter_abeta)-1),angle='90',y=1300)

#####################################################################################################
# meta-roi
#####################################################################################################

meta_roi_plot_data$RID <- as.numeric(meta_roi_plot_data$RID)

gam.model_meta_roi <- mgcv::gam(formula = meta_ROI ~ s(adjusted_new_time, bs = "cs", k = 4, fx = TRUE) +
                                  s(RID, bs = "re") + 
                                  age + PTEDUCAT + apoe + PTGENDER + CDGLOBAL,
                                data = meta_roi_plot_data,
                                family = gaussian,
                                method = "ML")

m2.d_meta_roi <- Deriv(gam.model_meta_roi, newdata = test) #Calculating derivatives of the model

m2.dci_meta_roi <- confint(m2.d_meta_roi, term = 'adjusted_new_time') #Calculating upper and lower confidence interval

m2.dsig_meta_roi <- signifD(seq(-10,10, by = 0.1), d = m2.d_meta_roi[['adjusted_new_time']]$deriv, m2.dci_meta_roi[['adjusted_new_time']]$upper, m2.dci_meta_roi[['adjusted_new_time']]$lower) # calculating the significance of change

modelFit_meta_roi <- data.frame(predict(gam.model_meta_roi, se = TRUE,
                                        data.frame(adjusted_new_time = meta_roi_plot_data$adjusted_new_time, 
                                                   RID = 0, 
                                                   age = 73, 
                                                   PTEDUCAT = 16, 
                                                   CDGLOBAL = 0,
                                                   PTGENDER = 2,
                                                   apoe = "E3")))

df.fit_meta_roi <- data.frame(cbind(adjusted_new_time = meta_roi_plot_data$adjusted_new_time,
                                    meas = meta_roi_plot_data$meta_ROI,
                                    fit = modelFit_meta_roi$fit,
                                    upperBound = modelFit_meta_roi$fit + 2 * modelFit_meta_roi$se.fit,
                                    lowerBound = modelFit_meta_roi$fit - 2 * modelFit_meta_roi$se.fit))

df.fit_meta_roi <- df.fit_meta_roi[!is.na(df.fit_meta_roi$meas),]

meas_fn_meta_roi <- "Meta-ROI"

gam.p_meta_roi <- ggplot(data=df.fit_meta_roi, aes(adjusted_new_time,meas)) + ylim(c(min(df.fit_meta_roi$lowerBound),max(df.fit_meta_roi$upperBound))) +
  geom_line(aes(adjusted_new_time, fit), color = "black", size=0.25, data=df.fit_meta_roi) +
  ggplot2::geom_line(aes(y=upperBound),se=FALSE, stat="smooth",linetype="dashed",alpha=.5) +
  ggplot2::geom_line(aes(y=lowerBound),se=FALSE, stat="smooth",linetype="dashed",alpha=.5) +
  theme_blank() + theme(axis.line = element_line(size = 0.5, colour = "black", linetype=1)) + 
  ggtitle("") + 
  theme(plot.title=element_text(hjust=0.5),text=element_text(size=15)) +
  theme(legend.position="none") +
  xlab("Years from Aβ+") + ylab(meas_fn_meta_roi)

inter_meta_roi <- na.contiguous(m2.dsig_meta_roi$decr)

# gam.p_meta_roi <- gam.p_meta_roi + 
#   geom_line(aes(adjusted_new_time, fit), color = "black", size=1.5, data=df.fit_meta_roi[df.fit_meta_roi$adjusted_new_time>=min(inter_meta_roi) & df.fit_meta_roi$adjusted_new_time<=max(inter_meta_roi),]) + 
#   ggplot2::geom_vline(aes(xintercept=min(inter_meta_roi)),linetype="dashed",color="black") + #inflection point
#   ggplot2::annotate(geom="text",label=paste0(round(min(inter_meta_roi),digits=1)),x=(min(inter_meta_roi)-1),angle='90',y=10)

#no significance

#####################################################################################################
# hippocampal volume
#####################################################################################################

hippocampal_volume_plot_data$RID <- as.numeric(hippocampal_volume_plot_data$RID)

gam.model_hippocampal_volume <- mgcv::gam(formula = hippocampal_volume ~ s(adjusted_new_time, bs = "cs", k = 4, fx = TRUE) +
                                            s(RID, bs = "re") + 
                                            age + PTEDUCAT + apoe + PTGENDER + CDGLOBAL,
                                          data = hippocampal_volume_plot_data,
                                          family = gaussian,
                                          method = "ML")

m2.d_hippocampal_volume <- Deriv(gam.model_hippocampal_volume, newdata = test) #Calculating derivatives of the model

m2.dci_hippocampal_volume <- confint(m2.d_hippocampal_volume, term = 'adjusted_new_time') #Calculating upper and lower confidence interval

m2.dsig_hippocampal_volume <- signifD(seq(-10,10, by = 0.1), d = m2.d_hippocampal_volume[['adjusted_new_time']]$deriv, m2.dci_hippocampal_volume[['adjusted_new_time']]$upper, m2.dci_hippocampal_volume[['adjusted_new_time']]$lower) # calculating the significance of change

modelFit_hippocampal_volume <- data.frame(predict(gam.model_hippocampal_volume, se = TRUE,
                                                  data.frame(adjusted_new_time = hippocampal_volume_plot_data$adjusted_new_time, 
                                                             RID = 0, 
                                                             age = 73, 
                                                             PTEDUCAT = 16, 
                                                             CDGLOBAL = 0,
                                                             PTGENDER = 2,
                                                             apoe = "E3")))

df.fit_hippocampal_volume <- data.frame(cbind(adjusted_new_time = hippocampal_volume_plot_data$adjusted_new_time,
                                              meas = hippocampal_volume_plot_data$hippocampal_volume,
                                              fit = modelFit_hippocampal_volume$fit,
                                              upperBound = modelFit_hippocampal_volume$fit + 2 * modelFit_hippocampal_volume$se.fit,
                                              lowerBound = modelFit_hippocampal_volume$fit - 2 * modelFit_hippocampal_volume$se.fit))

df.fit_hippocampal_volume <- df.fit_hippocampal_volume[!is.na(df.fit_hippocampal_volume$meas),]

meas_fn_hippocampal_volume <- "Hippocampal Volume"

gam.p_hippocampal_volume <- ggplot(data=df.fit_hippocampal_volume, aes(adjusted_new_time,meas)) + ylim(c(min(df.fit_hippocampal_volume$lowerBound),max(df.fit_hippocampal_volume$upperBound))) +
  geom_line(aes(adjusted_new_time, fit), color = "black", size=0.25, data=df.fit_hippocampal_volume) +
  ggplot2::geom_line(aes(y=upperBound),se=FALSE, stat="smooth",linetype="dashed",alpha=.5) +
  ggplot2::geom_line(aes(y=lowerBound),se=FALSE, stat="smooth",linetype="dashed",alpha=.5) +
  theme_blank() + theme(axis.line = element_line(size = 0.5, colour = "black", linetype=1)) + 
  ggtitle("") + 
  theme(plot.title=element_text(hjust=0.5),text=element_text(size=15)) +
  theme(legend.position="none") +
  xlab("Years from Aβ+") + ylab(meas_fn_hippocampal_volume)

inter_hippocampal_volume <- na.contiguous(m2.dsig_hippocampal_volume$decr)

gam.p_hippocampal_volume <- gam.p_hippocampal_volume + 
  geom_line(aes(adjusted_new_time, fit), color = "black", size=1.5, data=df.fit_hippocampal_volume[df.fit_hippocampal_volume$adjusted_new_time>=min(inter_hippocampal_volume) & df.fit_hippocampal_volume$adjusted_new_time<=max(inter_hippocampal_volume),]) + 
  ggplot2::geom_vline(aes(xintercept=min(inter_hippocampal_volume)),linetype="dashed",color="black") + #inflection point
  ggplot2::annotate(geom="text",label=paste0(round(min(inter_hippocampal_volume),digits=1)),x=(min(inter_hippocampal_volume)-1),angle='90',y=2200)

#####################################################################################################
# Ecog - Subject
#####################################################################################################

ecog_s_plot_data$RID <- as.numeric(ecog_s_plot_data$RID)

gam.model_ecog_s <- mgcv::gam(formula = EcogGlobal ~ s(adjusted_new_time, bs = "cs", k = 4, fx = TRUE) +
                                s(RID, bs = "re") + 
                                age + PTEDUCAT + apoe + PTGENDER + CDGLOBAL,
                              data = ecog_s_plot_data,
                              family = gaussian,
                              method = "ML")

m2.d_ecog_s <- Deriv(gam.model_ecog_s, newdata = test) #Calculating derivatives of the model

m2.dci_ecog_s <- confint(m2.d_ecog_s, term = 'adjusted_new_time') #Calculating upper and lower confidence interval

m2.dsig_ecog_s <- signifD(seq(-10,10, by = 0.1), d = m2.d_ecog_s[['adjusted_new_time']]$deriv, m2.dci_ecog_s[['adjusted_new_time']]$upper, m2.dci_ecog_s[['adjusted_new_time']]$lower) # calculating the significance of change

modelFit_ecog_s <- data.frame(predict(gam.model_ecog_s, se = TRUE,
                                      data.frame(adjusted_new_time = ecog_s_plot_data$adjusted_new_time, 
                                                 RID = 0, 
                                                 age = 73, 
                                                 PTEDUCAT = 16, 
                                                 CDGLOBAL = 0,
                                                 PTGENDER = 2,
                                                 apoe = "E3")))

df.fit_ecog_s <- data.frame(cbind(adjusted_new_time = ecog_s_plot_data$adjusted_new_time,
                                  meas = ecog_s_plot_data$EcogGlobal,
                                  fit = modelFit_ecog_s$fit,
                                  upperBound = modelFit_ecog_s$fit + 2 * modelFit_ecog_s$se.fit,
                                  lowerBound = modelFit_ecog_s$fit - 2 * modelFit_ecog_s$se.fit))

df.fit_ecog_s <- df.fit_ecog_s[!is.na(df.fit_ecog_s$meas),]

meas_fn_ecog_s <- "ECog - Subject"

gam.p_ecog_s <- ggplot(data=df.fit_ecog_s, aes(adjusted_new_time,meas)) + ylim(c(min(df.fit_ecog_s$lowerBound),max(df.fit_ecog_s$upperBound))) +
  geom_line(aes(adjusted_new_time, fit), color = "black", size=0.25, data=df.fit_ecog_s) +
  ggplot2::geom_line(aes(y=upperBound),se=FALSE, stat="smooth",linetype="dashed",alpha=.5) +
  ggplot2::geom_line(aes(y=lowerBound),se=FALSE, stat="smooth",linetype="dashed",alpha=.5) +
  theme_blank() + theme(axis.line = element_line(size = 0.5, colour = "black", linetype=1)) + 
  ggtitle("") + 
  theme(plot.title=element_text(hjust=0.5),text=element_text(size=15)) +
  theme(legend.position="none") +
  xlab("Years from Aβ+") + ylab(meas_fn_ecog_s)

inter_ecog_s <- na.contiguous(m2.dsig_ecog_s$incr)

gam.p_ecog_s <- gam.p_ecog_s + 
  geom_line(aes(adjusted_new_time, fit), color = "black", size=1.5, data=df.fit_ecog_s[df.fit_ecog_s$adjusted_new_time>=min(inter_ecog_s) & df.fit_ecog_s$adjusted_new_time<=max(inter_ecog_s),]) + 
  ggplot2::geom_vline(aes(xintercept=min(inter_ecog_s)),linetype="dashed",color="black") + #inflection point
  ggplot2::annotate(geom="text",label=paste0(round(min(inter_ecog_s),digits=1)),x=(min(inter_ecog_s)-1),angle='90',y=1.35)


#####################################################################################################
# Ecog - Study Partner
#####################################################################################################

ecog_p_plot_data$RID <- as.numeric(ecog_p_plot_data$RID)

gam.model_ecog_p <- mgcv::gam(formula = EcogGlobal ~ s(adjusted_new_time, bs = "cs", k = 4, fx = TRUE) +
                                s(RID, bs = "re") + 
                                age + PTEDUCAT + apoe + PTGENDER + CDGLOBAL,
                              data = ecog_p_plot_data,
                              family = gaussian,
                              method = "ML")

m2.d_ecog_p <- Deriv(gam.model_ecog_p, newdata = test) #Calculating derivatives of the model

m2.dci_ecog_p <- confint(m2.d_ecog_p, term = 'adjusted_new_time') #Calculating upper and lower confidence interval
m2.dsig_ecog_p <- signifD(seq(-10,10,by=0.1), d = m2.d_ecog_p[['adjusted_new_time']]$deriv, m2.dci_ecog_p[['adjusted_new_time']]$upper, m2.dci_ecog_p[['adjusted_new_time']]$lower) # calculating the significance of change

modelFit_ecog_p <- data.frame(predict(gam.model_ecog_p, se = TRUE,
                                      data.frame(adjusted_new_time = ecog_p_plot_data$adjusted_new_time, 
                                                 RID = 0, 
                                                 age = 73, 
                                                 PTEDUCAT = 16, 
                                                 CDGLOBAL = 0,
                                                 PTGENDER = 2,
                                                 apoe = "E3")))

df.fit_ecog_p <- data.frame(cbind(adjusted_new_time = ecog_p_plot_data$adjusted_new_time,
                                  meas = ecog_p_plot_data$EcogGlobal,
                                  fit = modelFit_ecog_p$fit,
                                  upperBound = modelFit_ecog_p$fit + 2 * modelFit_ecog_p$se.fit,
                                  lowerBound = modelFit_ecog_p$fit - 2 * modelFit_ecog_p$se.fit))

df.fit_ecog_p <- df.fit_ecog_p[!is.na(df.fit_ecog_p$meas),]

meas_fn_ecog_p <- "ECog - Study Partner"

gam.p_ecog_p <- ggplot(data=df.fit_ecog_p, aes(adjusted_new_time,meas)) + ylim(c(min(df.fit_ecog_p$lowerBound),max(df.fit_ecog_p$upperBound))) +
  geom_line(aes(adjusted_new_time, fit), color = "black", size=0.25, data=df.fit_ecog_p) +
  ggplot2::geom_line(aes(y=upperBound),se=FALSE, stat="smooth",linetype="dashed",alpha=.5) +
  ggplot2::geom_line(aes(y=lowerBound),se=FALSE, stat="smooth",linetype="dashed",alpha=.5) +
  theme_blank() + theme(axis.line = element_line(size = 0.5, colour = "black", linetype=1)) + 
  ggtitle("") + 
  theme(plot.title=element_text(hjust=0.5),text=element_text(size=15)) +
  theme(legend.position="none") +
  xlab("Years from Aβ+") + ylab(meas_fn_ecog_p)

inter_ecog_p <- na.contiguous(m2.dsig_ecog_p$incr)

gam.p_ecog_p <- gam.p_ecog_p + 
  geom_line(aes(adjusted_new_time, fit), color = "black", size=1.5, data=df.fit_ecog_p[df.fit_ecog_p$adjusted_new_time>=min(inter_ecog_p) & df.fit_ecog_p$adjusted_new_time<=max(inter_ecog_p),]) + 
  ggplot2::geom_vline(aes(xintercept=min(inter_ecog_p)),linetype="dashed",color="black") + #inflection point
  ggplot2::annotate(geom="text",label=paste0(round(min(inter_ecog_p),digits=1)),x=(min(inter_ecog_p)-1),angle='90',y=1.1)

#####################################################################################################
# ADAS13
#####################################################################################################
adas13_plot_data$RID <- as.numeric(adas13_plot_data$RID)

gam.model_adas13 <- mgcv::gam(formula = ADAS13 ~ s(adjusted_new_time, bs = "cs", k = 4, fx = TRUE) +
                                s(RID, bs = "re") + 
                                age + PTEDUCAT + apoe + PTGENDER + CDGLOBAL,
                              data = adas13_plot_data,
                              family = gaussian,
                              method = "ML")

m2.d_adas13 <- Deriv(gam.model_adas13, newdata = test) #Calculating derivatives of the model

m2.dci_adas13 <- confint(m2.d_adas13, term = 'adjusted_new_time') #Calculating upper and lower confidence interval

m2.dsig_adas13 <- signifD(seq(-10, 10, by = 0.1), d = m2.d_adas13[['adjusted_new_time']]$deriv, m2.dci_adas13[['adjusted_new_time']]$upper, m2.dci_adas13[['adjusted_new_time']]$lower) # calculating the significance of change

modelFit_adas13 <- data.frame(predict(gam.model_adas13, se = TRUE,
                                      data.frame(adjusted_new_time = adas13_plot_data$adjusted_new_time, 
                                                 RID = 0, 
                                                 age = 73, 
                                                 PTEDUCAT = 16, 
                                                 CDGLOBAL = 0,
                                                 PTGENDER = 2,
                                                 apoe = "E3")))

df.fit_adas13 <- data.frame(cbind(adjusted_new_time = adas13_plot_data$adjusted_new_time,
                                  meas = adas13_plot_data$ADAS13,
                                  fit = modelFit_adas13$fit,
                                  upperBound = modelFit_adas13$fit + 2 * modelFit_adas13$se.fit,
                                  lowerBound = modelFit_adas13$fit - 2 * modelFit_adas13$se.fit))

df.fit_adas13 <- df.fit_adas13[!is.na(df.fit_adas13$meas),]

meas_fn_adas13 <- "ADAS13"

gam.p_adas13 <- ggplot(data=df.fit_adas13, aes(adjusted_new_time,meas)) + ylim(c(min(df.fit_adas13$lowerBound),max(df.fit_adas13$upperBound))) +
  geom_line(aes(adjusted_new_time, fit), color = "black", size=0.25, data=df.fit_adas13) +
  ggplot2::geom_line(aes(y=upperBound),se=FALSE, stat="smooth",linetype="dashed",alpha=.5) +
  ggplot2::geom_line(aes(y=lowerBound),se=FALSE, stat="smooth",linetype="dashed",alpha=.5) +
  theme_blank() + theme(axis.line = element_line(size = 0.5, colour = "black", linetype=1)) + 
  ggtitle("") + 
  theme(plot.title=element_text(hjust=0.5),text=element_text(size=15)) +
  theme(legend.position="none") +
  xlab("Years from Aβ+") + ylab(meas_fn_adas13)

inter_adas13 <- na.contiguous(m2.dsig_adas13$incr)

gam.p_adas13 <- gam.p_adas13 + 
  geom_line(aes(adjusted_new_time, fit), color = "black", size=1.5, data=df.fit_adas13[df.fit_adas13$adjusted_new_time>=min(inter_adas13) & df.fit_adas13$adjusted_new_time<=max(inter_adas13),]) + 
  ggplot2::geom_vline(aes(xintercept=min(inter_adas13)),linetype="dashed",color="black") + #inflection point
  ggplot2::annotate(geom="text",label=paste0(round(min(inter_adas13),digits=1)),x=(min(inter_adas13)-1),angle='90',y=10)

#####################################################################################################
# MMSE
#####################################################################################################

mmse_plot_data$RID <- as.numeric(mmse_plot_data$RID)

gam.model_mmse <- mgcv::gam(formula = MMSE ~ s(adjusted_new_time, bs = "cs", k = 4, fx = TRUE) +
                              s(RID, bs = "re") + 
                              age + PTEDUCAT + apoe + PTGENDER- + CDGLOBAL,
                            data = mmse_plot_data,
                            family = gaussian,
                            method = "ML")

m2.d_mmse <- Deriv(gam.model_mmse, newdata = test) #Calculating derivatives of the model

m2.dci_mmse <- confint(m2.d_mmse, term = 'adjusted_new_time') #Calculating upper and lower confidence interval

m2.dsig_mmse <- signifD(seq(-10,10, by = 0.1), d = m2.d_mmse[['adjusted_new_time']]$deriv, m2.dci_mmse[['adjusted_new_time']]$upper, m2.dci_mmse[['adjusted_new_time']]$lower) # calculating the significance of change

modelFit_mmse <- data.frame(predict(gam.model_mmse, se = TRUE,
                                    data.frame(adjusted_new_time = mmse_plot_data$adjusted_new_time, 
                                               RID = 0, 
                                               age = 73, 
                                               PTEDUCAT = 16, 
                                               CDGLOBAL = 0,
                                               PTGENDER = 2,
                                               apoe = "E3")))

df.fit_mmse <- data.frame(cbind(adjusted_new_time = mmse_plot_data$adjusted_new_time,
                                meas = mmse_plot_data$MMSE,
                                fit = modelFit_mmse$fit,
                                upperBound = modelFit_mmse$fit + 2 * modelFit_mmse$se.fit,
                                lowerBound = modelFit_mmse$fit - 2 * modelFit_mmse$se.fit))

df.fit_mmse <- df.fit_mmse[!is.na(df.fit_mmse$meas),]

meas_fn_mmse <- "MMSE"

gam.p_mmse <- ggplot(data=df.fit_mmse, aes(adjusted_new_time,meas)) + ylim(c(min(df.fit_mmse$lowerBound),max(df.fit_mmse$upperBound))) +
  geom_line(aes(adjusted_new_time, fit), color = "black", size=0.25, data=df.fit_mmse) +
  ggplot2::geom_line(aes(y=upperBound),se=FALSE, stat="smooth",linetype="dashed",alpha=.5) +
  ggplot2::geom_line(aes(y=lowerBound),se=FALSE, stat="smooth",linetype="dashed",alpha=.5) +
  theme_blank() + theme(axis.line = element_line(size = 0.5, colour = "black", linetype=1)) + 
  ggtitle("") + 
  theme(plot.title=element_text(hjust=0.5),text=element_text(size=15)) +
  theme(legend.position="none") +
  xlab("Years from Aβ+") + ylab(meas_fn_mmse)

inter_mmse <- na.contiguous(m2.dsig_mmse$decr)

gam.p_mmse <- gam.p_mmse + 
  geom_line(aes(adjusted_new_time, fit), color = "black", size=1.5, data=df.fit_mmse[df.fit_mmse$adjusted_new_time>=min(inter_mmse) & df.fit_mmse$adjusted_new_time<=max(inter_mmse),]) + 
  ggplot2::geom_vline(aes(xintercept=min(inter_mmse)),linetype="dashed",color="black") + #inflection point
  ggplot2::annotate(geom="text",label=paste0(round(min(inter_mmse),digits=1)),x=(min(inter_mmse)-1),angle='90',y=27)

#####################################################################################################
# mPACCtrailsB
#####################################################################################################
mpacctrailsb_plot_data$RID <- as.numeric(mpacctrailsb_plot_data$RID)

gam.model_mpacctrailsb <- mgcv::gam(formula = mPACCtrailsB ~ s(adjusted_new_time, bs = "cs", k = 4, fx = TRUE) +
                                      s(RID, bs = "re") + 
                                      age + PTEDUCAT + apoe + PTGENDER + CDGLOBAL,
                                    data = mpacctrailsb_plot_data,
                                    family = gaussian,
                                    method = "ML")

m2.d_mpacctrailsb <- Deriv(gam.model_mpacctrailsb, newdata = test) #Calculating derivatives of the model

m2.dci_mpacctrailsb <- confint(m2.d_mpacctrailsb, term = 'adjusted_new_time') #Calculating upper and lower confidence interval

m2.dsig_mpacctrailsb <- signifD(seq(-10,10, by = 0.1), d = m2.d_mpacctrailsb[['adjusted_new_time']]$deriv, m2.dci_mpacctrailsb[['adjusted_new_time']]$upper, m2.dci_mpacctrailsb[['adjusted_new_time']]$lower) # calculating the significance of change

modelFit_mpacctrailsb <- data.frame(predict(gam.model_mpacctrailsb, se = TRUE,
                                            data.frame(adjusted_new_time = mpacctrailsb_plot_data$adjusted_new_time, 
                                                       RID = 0, 
                                                       age = 73, 
                                                       PTEDUCAT = 16, 
                                                       CDGLOBAL = 0,
                                                       PTGENDER = 2,
                                                       apoe = "E3")))

df.fit_mpacctrailsb <- data.frame(cbind(adjusted_new_time = mpacctrailsb_plot_data$adjusted_new_time,
                                        meas = mpacctrailsb_plot_data$mPACCtrailsB,
                                        fit = modelFit_mpacctrailsb$fit,
                                        upperBound = modelFit_mpacctrailsb$fit + 2 * modelFit_mpacctrailsb$se.fit,
                                        lowerBound = modelFit_mpacctrailsb$fit - 2 * modelFit_mpacctrailsb$se.fit))

df.fit_mpacctrailsb <- df.fit_mpacctrailsb[!is.na(df.fit_mpacctrailsb$meas),]

meas_fn_mpacctrailsb <- "mPACCtrailsB"

gam.p_mpacctrailsb <- ggplot(data=df.fit_mpacctrailsb, aes(adjusted_new_time,meas)) + ylim(c(min(df.fit_mpacctrailsb$lowerBound),max(df.fit_mpacctrailsb$upperBound))) +
  geom_line(aes(adjusted_new_time, fit), color = "black", size=0.25, data=df.fit_mpacctrailsb) +
  ggplot2::geom_line(aes(y=upperBound),se=FALSE, stat="smooth",linetype="dashed",alpha=.5) +
  ggplot2::geom_line(aes(y=lowerBound),se=FALSE, stat="smooth",linetype="dashed",alpha=.5) +
  theme_blank() + theme(axis.line = element_line(size = 0.5, colour = "black", linetype=1)) + 
  ggtitle("") + 
  theme(plot.title=element_text(hjust=0.5),text=element_text(size=15)) +
  theme(legend.position="none") +
  xlab("Years from Aβ+") + ylab(meas_fn_mpacctrailsb)

inter_mpacctrailsb <- na.contiguous(m2.dsig_mpacctrailsb$decr)

gam.p_mpacctrailsb <- gam.p_mpacctrailsb + 
  geom_line(aes(adjusted_new_time, fit), color = "black", size=1.5, data=df.fit_mpacctrailsb[df.fit_mpacctrailsb$adjusted_new_time>=min(inter_mpacctrailsb) & df.fit_mpacctrailsb$adjusted_new_time<=max(inter_mpacctrailsb),]) + 
  ggplot2::geom_vline(aes(xintercept=min(inter_mpacctrailsb)),linetype="dashed",color="black") + #inflection point
  ggplot2::annotate(geom="text",label=paste0(round(min(inter_mpacctrailsb),digits=1)),x=(min(inter_mpacctrailsb)-1),angle='90',y=-.3)


#####################################################################################################
# CDRSB
#####################################################################################################
cdrsb_plot_data$RID <- as.numeric(cdrsb_plot_data$RID)

gam.model_cdrsb <- mgcv::gam(formula = CDRSB ~ s(adjusted_new_time, bs = "cs", k = 4, fx = TRUE) +
                               s(RID, bs = "re") + 
                               age + PTEDUCAT + apoe + PTGENDER + CDGLOBAL,
                             data = cdrsb_plot_data,
                             family = gaussian,
                             method = "ML")

m2.d_cdrsb <- Deriv(gam.model_cdrsb, newdata = test) #Calculating derivatives of the model

m2.dci_cdrsb <- confint(m2.d_cdrsb, term = 'adjusted_new_time') #Calculating upper and lower confidence interval

m2.dsig_cdrsb <- signifD(seq(-10,10,by=0.1), d = m2.d_cdrsb[['adjusted_new_time']]$deriv, m2.dci_cdrsb[['adjusted_new_time']]$upper, m2.dci_cdrsb[['adjusted_new_time']]$lower) # calculating the significance of change

modelFit_cdrsb <- data.frame(predict(gam.model_cdrsb, se = TRUE,
                                     data.frame(adjusted_new_time = cdrsb_plot_data$adjusted_new_time, 
                                                RID = 0, 
                                                age = 73, 
                                                PTEDUCAT = 16, 
                                                PTGENDER = 2,
                                                apoe = "E3",
                                                CDGLOBAL = 0)))

df.fit_cdrsb <- data.frame(cbind(adjusted_new_time = cdrsb_plot_data$adjusted_new_time,
                                 meas = cdrsb_plot_data$CDRSB,
                                 fit = modelFit_cdrsb$fit,
                                 upperBound = modelFit_cdrsb$fit + 2 * modelFit_cdrsb$se.fit,
                                 lowerBound = modelFit_cdrsb$fit - 2 * modelFit_cdrsb$se.fit))

df.fit_cdrsb <- df.fit_cdrsb[!is.na(df.fit_cdrsb$meas),]

meas_fn_cdrsb <- "CDR-SB"

gam.p_cdrsb <- ggplot(data=df.fit_cdrsb, aes(adjusted_new_time,meas)) + ylim(c(min(df.fit_cdrsb$lowerBound),max(df.fit_cdrsb$upperBound))) +
  geom_line(aes(adjusted_new_time, fit), color = "black", size=0.25, data=df.fit_cdrsb) +
  ggplot2::geom_line(aes(y=upperBound),se=FALSE, stat="smooth",linetype="dashed",alpha=.5) +
  ggplot2::geom_line(aes(y=lowerBound),se=FALSE, stat="smooth",linetype="dashed",alpha=.5) +
  theme_blank() + theme(axis.line = element_line(size = 0.5, colour = "black", linetype=1)) + 
  ggtitle("") + 
  theme(plot.title=element_text(hjust=0.5),text=element_text(size=15)) +
  theme(legend.position="none") +
  xlab("Years from Aβ+") + ylab(meas_fn_cdrsb)

inter_cdrsb <- na.contiguous(m2.dsig_cdrsb$incr)

gam.p_cdrsb <- gam.p_cdrsb + 
  geom_line(aes(adjusted_new_time, fit), color = "black", size=1.5, data=df.fit_cdrsb[df.fit_cdrsb$adjusted_new_time>=min(inter_cdrsb) & df.fit_cdrsb$adjusted_new_time<=max(inter_cdrsb),]) + 
  ggplot2::geom_vline(aes(xintercept=min(inter_cdrsb)),linetype="dashed",color="black") + #inflection point
  ggplot2::annotate(geom="text",label=paste0(round(min(inter_cdrsb),digits=1)),x=(min(inter_cdrsb)-1),angle='90',y= 0.4)

#####################################################################################################
# plotting
#####################################################################################################

csf_centiloid <- gam.p_centiloid + gam.p_abeta + gam.p_ptau + gam.p_tau + plot_layout(nrow = 2)

scans <- gam.p_hippocampal_volume + gam.p_meta_roi + gam.p_fdg + plot_layout(nrow = 1)

(csf_centiloid) / (scans) + plot_layout(heights = c(2, 1)) + plot_annotation(tag_levels = 'a') &  theme(plot.tag = element_text(face = 'bold'))

gam.p_cdrsb + gam.p_adas13 + gam.p_mmse + gam.p_mpacctrailsb + gam.p_ecog_p + gam.p_ecog_s + plot_layout(nrow = 2) + plot_annotation(tag_levels = 'a') &  theme(plot.tag = element_text(face = 'bold'))
