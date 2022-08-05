#loading libraries
library(dplyr)
library(lubridate)
library(tidyr)
library(naniar)

########################################################### First prepping all the cross-sectional data

test_adni1 <- read.csv("C:\\Documents\\paper data longitudinal phases\\UCSFFSX_11_02_15.csv") %>% #ST41SV
  mutate(RID = as.factor(RID)) %>%
  mutate(EXAMDATE = as.POSIXct(EXAMDATE, format = "%Y-%m-%d"))

test_adni2 <- read.csv("C:\\Documents\\paper data longitudinal phases\\UCSFFSX6_06_07_22.csv") %>%
  mutate(RID = as.factor(RID)) %>%
  mutate(EXAMDATE = as.POSIXct(EXAMDATE, format = "%Y-%m-%d"))
         
test_adni3 <- read.csv("C:\\Documents\\paper data longitudinal phases\\UCSFFSX51_11_08_19.csv") %>%
  mutate(RID = as.factor(RID)) %>%
  mutate(EXAMDATE = as.POSIXct(EXAMDATE, format = "%Y-%m-%d"))
         
scans_3T_Duygu <- read.csv("C:\\Documents\\scanner information\\UCSF - ADNI-1 3T Cross-Sectional FreeSurfer (5.1) [ADNI1].csv") %>%
  mutate(FLDSTRENG = 3)%>%
  mutate(RID = as.factor(RID)) %>%
  mutate(EXAMDATE = as.POSIXct(EXAMDATE, format = "%Y-%m-%d"))

# creating new merged dataframe
ROI_data_list <- list(test_adni1, test_adni2, test_adni3, scans_3T_Duygu)

ROI_data_merged <- test_adni1
for ( df in ROI_data_list ) {
  ROI_data_merged <-merge(ROI_data_merged, df, all=T)
}

#removing the individual datasets
rm(test_adni1, test_adni2, test_adni3, scans_3T_Duygu, ROI_data_list, df)

# making sure that the date is in date format
ROI_data_merged$EXAMDATE <- as.POSIXct(ROI_data_merged$EXAMDATE, format = "%Y-%m-%d")

#some brain regions only have na values so getting rid of columns with na values only
ROI_data_merged <- ROI_data_merged %>%
  select_if(~sum(!is.na(.)) > 0)
  
########################################################### First prepping all the scanner meta data

scan_test <- read.csv("C:\\Documents\\scanner information\\1.5T scan information.csv") %>%
  select(PHASE, EXAMDATE, FIELD_STRENGTH, RID, SITEID, VISCODE)
# putting the dates in the same format
multidate <- function(data, formats){
  a<-list()
  for(i in 1:length(formats)){
    a[[i]]<- as.Date(data,format=formats[i])
    a[[1]][!is.na(a[[i]])]<-a[[i]][!is.na(a[[i]])]
  }
  a[[1]]
}
scan_test$EXAMDATE <- multidate(scan_test$EXAMDATE, 
                                       c("%m/%d/%Y","%Y-%m-%d %H:%M:%S","%Y-%m-%d"))

scan_test2 <- read.csv("C:\\Documents\\scanner information\\adni_scanner_meta_data.csv")  %>%
  dplyr::mutate(RID = case_when(grepl("ADNI2",ScanCode) | grepl("ADNI3", ScanCode) | grepl("ADNID", ScanCode) ~ sapply(strsplit(ScanCode, "_"), function(x) x[4]))) %>%
  dplyr::mutate(project = case_when(grepl("ADNI2", ScanCode) ~ "ADNI2",
                                    grepl("ADNI3", ScanCode) ~ "ADNI3")) %>%
  dplyr::mutate(RID = case_when(grepl("ADNI", ScanCode) ~ as.numeric(substr(RID, 1, 4)))) %>%
  dplyr::mutate(RID = as.character(RID)) %>%
  dplyr::mutate(RID = as.factor(RID)) %>%
  mutate(EXAMDATE = as.POSIXct(date, format = "%Y-%m-%d")) %>%
  select(-c(date, StudyInstanceUID, SeriesInstanceUID, DeviceSerialNumber))
# scan_test2[scan_test2 == "NULL"] = NA

scan_test3 <- read.csv("C:\\Documents\\paper data longitudinal phases\\MRI3META.csv") %>%
  select(RID, EXAMDATE, PHASE, FIELD_STRENGTH, SITEID, VISCODE) %>%
  mutate(EXAMDATE = as.POSIXct(EXAMDATE, format = "%Y-%m-%d"))

scan_test4 <- read.csv("C:\\Documents\\paper data longitudinal phases\\MRIMETA.csv") %>%
  select(PHASE, EXAMDATE, FIELD_STRENGTH, RID, SITEID, VISCODE) %>%
  mutate(EXAMDATE = as.POSIXct(EXAMDATE, format = "%Y-%m-%d"))

#now combining the scanner data
scan_data_list <- list(scan_test, scan_test2, scan_test3, scan_test4)

scan_data_merged <- scan_test
for ( df in scan_data_list ) {
  scan_data_merged <-merge(scan_data_merged, df, all=T)
} #42903

#making RID a factor
scan_data_merged$RID <- as.factor(scan_data_merged$RID)

#now getting rid of the individual scanner datasets
rm(scan_test, scan_test2, scan_test3, scan_test4, scan_data_list, df)

#filling in information about scanners if the ID, date, and field strength matches
scan_data_merged <- scan_data_merged %>%
  mutate(FIELD_STRENGTH = case_when(is.na(FIELD_STRENGTH) & MagneticFieldStrength == 3 ~ "3T",
                                    is.na(FIELD_STRENGTH) & MagneticFieldStrength == 1.5 ~ "1.5T",
                                    TRUE ~ FIELD_STRENGTH))

scan_data_merged <- scan_data_merged %>%
  mutate(Sex = case_when(PatientSex == "Female" ~ "F",
                         PatientSex == "Male" ~ "M",
                         TRUE ~ PatientSex),
         ADNI_Phase = case_when(is.na(PHASE) ~ project,
                                TRUE ~ PHASE)) %>%
  group_by(RID) %>% 
  fill(Sex, .direction = "downup") %>%
  ungroup()


scan_data_merged <- scan_data_merged %>% 
  group_by(RID, EXAMDATE, FIELD_STRENGTH) %>%
  fill(c(Manufacturer, ManufacturerModelName, SoftwareVersion), .direction = "downup") %>%
  ungroup() %>%
  unique()

#renaming some of thte columns in the scanning data (for sanity check that it is merging as it should) and dropping unnecessary columns
scan_data_merged <- rename(scan_data_merged, Field_Strength = FIELD_STRENGTH)
scan_data_merged <- rename(scan_data_merged, Code = ScanCode)

scan_data_merged <- scan_data_merged %>%
  select(-c(PHASE, project, PatientSex, MagneticFieldStrength))

scan_data_merged$EXAMDATE <- as.POSIXct(scan_data_merged$EXAMDATE, format = "Y-%m-%d")
scan_data_merged$EXAMDATE <- format(as.POSIXct(scan_data_merged$EXAMDATE,format='%-/%m-%d %H:%M:%S'),format='%Y-%m-%d') #getting rid of the timestamp

########################################################### Now filtering the participants to only participants who have been in the study at least 10 years

ROI_data_filtered <- ROI_data_merged %>%
  group_by(RID) %>%
  mutate(oldest_date = min(EXAMDATE),
         newest_date = max(EXAMDATE)) %>%
  ungroup()

ROI_data_filtered <- ROI_data_filtered %>%
  mutate(ADNI_years = as.numeric(difftime(newest_date, oldest_date, unit = "weeks"))/52.25)

ROI_data_filtered <- ROI_data_filtered %>%
  filter(ADNI_years >= 10) %>%
  unique()
length(unique(ROI_data_filtered$RID)) #86 people who have been in at least a decade

#looking at the scanner data with only the relevant ID's
RIDs_of_interest <- unique(ROI_data_filtered$RID)

scan_data_merged <- scan_data_merged %>%
  filter(RID %in% RIDs_of_interest)

########################################################### Now merging ROI volumes with scan data

#renaming some columns in ROI data for a sanity check after merging
ROI_data_filtered <- rename(ROI_data_filtered, Visit_code = VISCODE)
ROI_data_filtered <- rename(ROI_data_filtered, LONIID = LONIUID)
ROI_data_filtered <- rename(ROI_data_filtered, IMAGEID = IMAGEUID)
ROI_data_filtered <- rename(ROI_data_filtered, SID = LONISID)

ROI_data_filtered <- ROI_data_filtered %>%
  select(-c(RUNDATE, VERSION))

#getting the scan data ready for merging
scan_data_merged$EXAMDATE <- as.POSIXct(scan_data_merged$EXAMDATE, format = "%Y-%m-%d") + days(1) #dates seem to be off by one day in the scan data vs ROI data
scan_data_merged <- scan_data_merged %>%
  filter(!(is.na(EXAMDATE) | is.na(RID)))

#now merging the data a few different ways
merged_data <- merge(ROI_data_filtered, scan_data_merged, all.x = TRUE) #%>% #this will merge by RID and EXAMDATE
  #mutate(FIELD_STRENGTH = case_when())

merged_data <- merged_data[!duplicated(merged_data[,c('RID','EXAMDATE', "FLDSTRENG", "Field_Strength")]),] %>%
  filter((FLDSTRENG == 1.5 & Field_Strength == "1.5T") | (FLDSTRENG == 3 & Field_Strength == "3T") | is.na(FLDSTRENG) | is.na(Field_Strength)) %>%
  mutate(FIELD_STRENGTH = case_when(FLDSTRENG == 1.5 ~ "1.5T",
                                    FLDSTRENG == 3 ~ "3T",
                                    Field_Strength == "1.5T" ~ "1.5T",
                                    Field_Strength == "3T" ~ "3T")) %>%
  filter(!is.na(FIELD_STRENGTH))

#now merging the data in different ways to get the most accurate matching
merged_data <- merged_data[!duplicated(merged_data[,c('RID','EXAMDATE', "FIELD_STRENGTH")]),]


merged_data_by_viscode <- merge(ROI_data_filtered, scan_data_merged, all.x = TRUE, by.x = c("RID", "Visit_code"), by.y = c("RID", "VISCODE")) %>% #this will merge by RID and VISCODE - i noticed that sometimes the dates were off but the codes were the same - making the assumption these are the same scans
  filter((FLDSTRENG == 1.5 & Field_Strength == "1.5T") | (FLDSTRENG == 3 & Field_Strength == "3T") | is.na(FLDSTRENG) | is.na(Field_Strength)) %>%
  mutate(FIELD_STRENGTH = case_when(FLDSTRENG == 1.5 ~ "1.5T",
                                    FLDSTRENG == 3 ~ "3T",
                                    Field_Strength == "1.5T" ~ "1.5T",
                                    Field_Strength == "3T" ~ "3T")) %>%
  filter(!is.na(FIELD_STRENGTH))

merged_data_by_viscode <- merged_data_by_viscode[!duplicated(merged_data_by_viscode[,c('RID','EXAMDATE.x', "EXAMDATE.y", "FIELD_STRENGTH")]),]

merged_data_by_id <- merge(ROI_data_filtered, scan_data_merged, all.x = TRUE, by.x = c("RID", "LONIID"), by.y = c("RID", "LoniSeriesID")) %>% #this will merge by RID and SID - i noticed that sometimes the dates were off but the SID's were the same - making the assumption these are the same scans
  filter((FLDSTRENG == 1.5 & Field_Strength == "1.5T") | (FLDSTRENG == 3 & Field_Strength == "3T") | is.na(FLDSTRENG) | is.na(Field_Strength)) %>%
  mutate(FIELD_STRENGTH = case_when(FLDSTRENG == 1.5 ~ "1.5T",
                                    FLDSTRENG == 3 ~ "3T",
                                    Field_Strength == "1.5T" ~ "1.5T",
                                    Field_Strength == "3T" ~ "3T")) %>%
  filter(!is.na(FIELD_STRENGTH))

merged_data_by_id <- merged_data_by_id[!duplicated(merged_data_by_id[,c('RID','EXAMDATE.x', "EXAMDATE.y", "FIELD_STRENGTH")]),]

#making temporary datasets to filter to only the cases we want in the final dataset
merged_data_temp <- merged_data %>%
  select(RID, EXAMDATE, FLDSTRENG, Field_Strength, FIELD_STRENGTH, Visit_code, VISCODE, LoniSeriesID, LONIID, Manufacturer, ManufacturerModelName, SoftwareVersion, Sex)

merged_data_viscode_temp <- merged_data_by_viscode %>%
  select(RID, EXAMDATE.x, EXAMDATE.y, FLDSTRENG, Field_Strength, FIELD_STRENGTH, Visit_code, LoniSeriesID, LONIID, Manufacturer, ManufacturerModelName, SoftwareVersion, Sex) %>%
  rename(VISIT_CODE = Visit_code)

merged_data_id_temp <- merged_data_by_id %>%
  select(RID, EXAMDATE.x, EXAMDATE.y, FLDSTRENG, Field_Strength, FIELD_STRENGTH, Visit_code, LONIID, Manufacturer, ManufacturerModelName, SoftwareVersion, Sex) %>%
  rename(SID = LONIID)

# now merging those temporary datasets

all_merged_temp <- merge(merged_data_temp, merged_data_viscode_temp, all = TRUE) %>%
  merge(merged_data_id_temp, all = TRUE) %>%
  naniar::replace_with_na(replace = list(Visit_code = c("", " "))) %>%
  mutate(LoniSeriesID = as.character(LoniSeriesID),
         LONIID = as.character(LONIID),
         SID = as.character(SID))

all_merged_temp <- all_merged_temp %>%
  mutate(examDate = case_when(is.na(EXAMDATE) ~ EXAMDATE.x,
                              TRUE ~ EXAMDATE)) %>%
  mutate(Visit_Code = case_when(is.na(Visit_code) & is.na(VISCODE) & !is.na(VISIT_CODE) ~ VISIT_CODE,
                                is.na(Visit_code) & !is.na(VISCODE) & !is.na(VISIT_CODE) ~ VISIT_CODE,
                                !is.na(Visit_code) & is.na(VISCODE) & !is.na(VISIT_CODE) ~ VISIT_CODE,
                                !is.na(Visit_code) & !is.na(VISCODE) & !is.na(VISIT_CODE) ~ VISIT_CODE,
                                is.na(VISIT_CODE) & is.na(VISCODE) & !is.na(Visit_code) ~ Visit_code,
                                is.na(VISIT_CODE) & !is.na(VISCODE) & !is.na(Visit_code) ~ Visit_code,
                                TRUE ~ VISCODE)) %>%
  mutate(LONISID = case_when(is.na(LONIID) & is.na(LoniSeriesID) & !is.na(SID) ~ SID,
                             is.na(LONIID) & !is.na(LoniSeriesID) & !is.na(SID) ~ SID,
                             !is.na(LONIID) & is.na(LoniSeriesID) & !is.na(SID) ~ SID,
                             !is.na(LONIID) & !is.na(LoniSeriesID) & !is.na(SID) ~ SID,
                             is.na(SID) & is.na(LoniSeriesID) & !is.na(LONIID) ~ LONIID,
                             is.na(SID) & !is.na(LoniSeriesID) & !is.na(LONIID) ~ LONIID,
                             TRUE ~ LoniSeriesID))
all_merged_temp <- all_merged_temp %>%
  rename(LONISERIESID = LONISID)

all_merged_temp_test <- all_merged_temp %>%
   select(-c(FLDSTRENG, Field_Strength, EXAMDATE, EXAMDATE.x, EXAMDATE.y)) %>% #, Visit_code, VISIT_CODE, VISCODE, LONIID, LoniSeriesID, SID)) %>%
   unique()

########################################################### now merging this data back in with the filtered data
filtered_merge <- merge(ROI_data_filtered, all_merged_temp_test, by.x = c("RID", "EXAMDATE"), by.y = c("RID", "examDate"), all.x = TRUE) %>%
  filter((FIELD_STRENGTH == "1.5T" & FLDSTRENG == 1.5) | (FIELD_STRENGTH == "3T" & FLDSTRENG == 3) | is.na(FIELD_STRENGTH) | is.na(FLDSTRENG)) %>%
  select(-c(Visit_code.y, Visit_Code, VISCODE, VISIT_CODE, LONIID.y, LoniSeriesID, LONISERIESID, SID.x, SID.y, FLDSTRENG))%>%
  naniar::replace_with_na(replace = list(Visit_code.x = c("", " "))) %>%
  unique()
