#can also use SUVR 1.42 as a cutoff for amyloid
library(dplyr)
library(ggplot2)
library(segmented)
library(lcmm)
library(ggrepel)

#####################################################################################
# making a new dataset to make predictions on
#####################################################################################
set.seed(123)
datnew   <- data.frame(adjusted_new_time = seq(-8, 8, length = 200),
                       age = round(seq(55.7, 95.4, length = 100), 1),
                       apoe = sample(c("E2", "E3", "E4"), 100, prob = c(0.06, 0.56, 0.38), replace = TRUE),
                       PTGENDER = sample(c(1, 2), 100, prob = c(0.45, 0.55), replace = TRUE),
                       PTEDUCAT = sample(c(12, 13, 14, 15, 16, 17, 18, 19, 20), 100, prob = c(0.11, 0.03, 0.04, 0.05, 0.26, 0.05, 0.22, 0.08, 0.14), replace = TRUE)
)

#####################################################################################
#getting the participant ID's that have changed from amyloid negative to amyloid positive
#####################################################################################

amyloid_pet <- readr::read_delim("C:\\Work Folder\\Amprion\\Source Data\\Source Data\\UCBERKELEY_AMY_6MM_30Sep2023.csv") [,1:14]
amyloid_pet <- amyloid_pet %>% dplyr::rename(suvr_summary=SUMMARY_SUVR,
                                             Centiloid=CENTILOIDS,
                                             AmyloidPosPET=AMYLOID_STATUS,
                                             EXAMDATE_pet = SCANDATE) %>%
  dplyr::mutate(EXAMDATE_pet = as.Date(EXAMDATE_pet),
                RID = as.character(RID)) %>%
  dplyr::select(RID, EXAMDATE_pet, AmyloidPosPET, Centiloid, suvr_summary)

change_in_amyloid <- amyloid_pet %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(all(c(0, 1) %in% AmyloidPosPET)) %>%
  # dplyr::select(-contains("_bl")) %>%
  dplyr::distinct()

amyloid_RIDs <- unique(change_in_amyloid$RID)

#getting important dates for all ids
last_a_negative_date_pet <- amyloid_pet %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(AmyloidPosPET == 0) %>%
  dplyr::mutate(last_a_neg_date_pet = max(EXAMDATE_pet)) %>%
  dplyr::select(RID, last_a_neg_date_pet) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()

#PET - getting first A+ dates by RID
first_a_positive_date_pet <- amyloid_pet %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(AmyloidPosPET == 1) %>%
  dplyr::mutate(first_a_pos_date_pet = min(EXAMDATE_pet)) %>%
  dplyr::select(RID, first_a_pos_date_pet) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()

change_in_amyloid <- merge(change_in_amyloid, last_a_negative_date_pet, by = "RID")
change_in_amyloid <- merge(change_in_amyloid, first_a_positive_date_pet, by = "RID")

change_in_amyloid_filtered <- change_in_amyloid %>%
  dplyr::filter(EXAMDATE_pet == last_a_neg_date_pet | EXAMDATE_pet == first_a_pos_date_pet)

#####################################################################################
# building a date to frame new x-axis
#####################################################################################
lm_data <- data.frame(matrix(ncol = 2, nrow = 0))
names(lm_data) <- c("RID", "predicted_date")

for (id_num in 1:length(unique(change_in_amyloid_filtered$RID))) {
  
  id <- unique(change_in_amyloid_filtered$RID)[id_num]
  
  temp_data_model <- change_in_amyloid_filtered %>%
    dplyr::filter(RID == id)
  
  temp_lm <- lm(EXAMDATE_pet ~ Centiloid, data = temp_data_model)
  
  predicted_value <- predict(temp_lm, data.frame(Centiloid = 20))
  
  temp_data_predictions <- data.frame(matrix(ncol = 2, nrow = 1))
  names(temp_data_predictions) <- c("RID", "predicted_date")
  
  temp_data_predictions$RID <- id
  temp_data_predictions$predicted_date <- as.Date(predicted_value)
  
  # predicted_value <- as.Date(predicted_value)
  
  # predicted_values_lm <- cbind(predicted_values_lm, as.Date(predicted_value))
  lm_data <- rbind(lm_data, temp_data_predictions)
  
}

#####################################################################################
# looking at change in amyloid status
#####################################################################################
centiloid_plot_data <- merge(change_in_amyloid, lm_data, all.x = TRUE) %>%
  dplyr::mutate(adjusted_new_time = lubridate::time_length(difftime(EXAMDATE_pet, predicted_date), "years"))

#getting rid of the RID's that have weird Centiloid patterns
wonky_RIDs <- centiloid_plot_data %>%
  filter(adjusted_new_time < 0 & Centiloid > 20)

centiloid_plot_data <- centiloid_plot_data %>%
  dplyr::filter(!(RID %in% unique(wonky_RIDs$RID) | RID == "1261" | RID == "6234"))

#############################################################################
# getting APOE, education, and birth year into a dataset
#############################################################################
#getting APOE type
apoeres <- read.csv("C:\\Work Folder\\paper data longitudinal phases\\APOERES.csv") %>% #OR Local: "C:\\Documents\\paper data longitudinal phases\\APOERES.csv"
  dplyr::filter(RID %in% unique(centiloid_plot_data$RID)) %>%
  dplyr::select(RID, APGEN1, APGEN2) %>%
  dplyr::mutate(apoe = paste(paste("E", APGEN1, sep = ""), paste("E", APGEN2, sep = ""), sep = "/"),
                RID = as.character(RID)) %>%
  dplyr::distinct()

#getting information to calculate age
dem <- read.csv("C:\\Work Folder\\paper data longitudinal phases\\demographics.csv") %>%
  dplyr::mutate(day = 1,
                birthdate = as.character(paste(day, PTDOB, sep = "/")))

dem <- dem %>%
  dplyr::mutate(birthdate = as.Date(birthdate, format = "%d/%m/%Y")) %>%
  dplyr::select(RID, PTGENDER, PTEDUCAT, birthdate) %>%
  dplyr::distinct() %>%
  dplyr::filter(!is.na(birthdate))

dem <- merge(dem, apoeres, all = TRUE) %>%
  dplyr::mutate(apoe = case_when(apoe == "E2/E3" | apoe == "E2/E2" ~ "E2",
                                 apoe == "E3/E3" ~ "E3",
                                 apoe == "E3/E4" | apoe == "E4/E4" ~ "E4"),
                RID = as.character(RID))

dem_all <- dem

dem <- dem %>%
  dplyr::filter(RID %in% unique(centiloid_plot_data$RID))
  

#############################################################################
#
#############################################################################

#adding diagnoses as CN, MCI, or AD
diagnoses <- read.csv("C:\\Work Folder\\Data Tables\\Diagnosis\\DXSUM_PDXCONV_ADNIALL.csv") # OR Local: "C:\\Documents\\Data Tables\\Diagnosis\\DXSUM_PDXCONV_ADNIALL.csv"
diagnoses[diagnoses == ""] <- NA

diagnoses_all <- diagnoses %>%
  dplyr::select(RID, EXAMDATE, DIAGNOSIS, DXCHANGE, DXCURREN) %>%
  dplyr::mutate(diags = dplyr::case_when(DIAGNOSIS==1 | DXCHANGE==1 | DXCHANGE==7 | DXCHANGE==9 | DXCURREN == 1 ~ 'CN',
                                         DIAGNOSIS==2 | DXCHANGE==2 | DXCHANGE==4 | DXCHANGE==8 | DXCURREN == 2 ~ 'MCI',
                                         DIAGNOSIS==3 | DXCHANGE==3 | DXCHANGE==5 | DXCHANGE==6 | DXCURREN == 3 ~ 'AD'),
                RID = as.character(RID)) %>%
  dplyr::rename(DX.DATE = EXAMDATE)

diagnoses <- diagnoses_all %>%
  dplyr::filter(RID %in% centiloid_plot_data$RID)

diagnoses_bl <- diagnoses %>%
  dplyr::group_by(RID) %>%
  dplyr::arrange(as.Date(DX.DATE)) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

diagnoses_at_0 <- diagnoses %>%
  dplyr::left_join(lm_data) %>%
  dplyr::mutate(date_diff = abs(as.Date(predicted_date) - as.Date(DX.DATE))) %>%
  dplyr::filter(!is.na(DX.DATE))

diagnoses_at_0 <- diagnoses_at_0 %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(date_diff == min(date_diff)) %>%
  dplyr::ungroup()

#############################################################################
#continuing with Centiloid
#############################################################################

centiloid_plot_data <- as.data.frame(centiloid_plot_data) %>%
  dplyr::left_join(dem, by = "RID") %>%
  dplyr::mutate(age = lubridate::time_length(difftime(EXAMDATE_pet, birthdate), "years"))
centiloid_plot_data <- centiloid_plot_data %>%
  dplyr::mutate(age = round(age, digits = 1))
centiloid_plot_data$RID <- as.numeric(centiloid_plot_data$RID)

centiloid_splines <- lcmm(Centiloid ~ adjusted_new_time + age + age*adjusted_new_time + PTGENDER + PTEDUCAT + apoe, #link = c("5-quantsplines"), 
                          random = ~adjusted_new_time, subject="RID", maxiter = 300, data = centiloid_plot_data, link = "3-equi-splines")

summary(centiloid_splines) #3-spline: 2385.50

# numDeriv::grad(centiloid_splines)

centiloid_predict_splines <- predictlink(centiloid_splines)

centiloid_test <- predictY(centiloid_splines, datnew, var.time = "adjusted_new_time", draws = TRUE)

centiloid_bootstrapped_data <- lcmm_bootstrap_ci(new_data = datnew, n_iterations = 100, lcmm_data = centiloid_plot_data, name_of_biomarker = "Centiloid")
