# Example
data$batch <- case_when(data$ADNI_Phase == "ADNI1" & kdata$FIELD_STRENGTH == "1.5T" ~ "Protocol 1",
                            data$ADNI_Phase == "ADNI2" & data$FIELD_STRENGTH == "1.5T" ~ "Protocol 1",
                            data$ADNI_Phase == "ADNIGO" & data$FIELD_STRENGTH == "1.5T" ~ "Protocol 1",
                            data$ADNI_Phase == "ADNI2" & data$FIELD_STRENGTH == "3T" ~ "Protocol 2",
                            data$ADNI_Phase == "ADNIGO" & data$FIELD_STRENGTH == "3T" ~ "Protocol 2",
                            data$ADNI_Phase == "ADNI3" & data$FIELD_STRENGTH == "3T" ~ "Protocol 3")
