############################################################################################
############################################################################################
################# SCRIPT 8 - Danni Gadd - Blood vs brain DNAm project ######################
############################################################################################
############################################################################################

setwd("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts")
library(readxl)
library(plyr)
library(tibble)
library(tidyverse)
library(dplyr)

############################################################################################

### LOAD DATA

############################################################################################

## LBC 14 individuals subject IDs
IDs <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/Data/post_mortem_DNAm_extraction_Anna_14April2020/LBC_IDs.csv")
IDs$LBC_ID <- as.character(IDs$LBC_ID)

## Brain and blood methylation datasets 
brain_meth <- readRDS("/Volumes/marioni-lab/Danni/Blood_brain_project/Data/post_mortem_DNAm_extraction_Anna_14April2020/brain_methylation_data.rds")
blood_meth <- readRDS("/Volumes/marioni-lab/Danni/Blood_brain_project/Data/pre_mortem_DNAm_extraction_REM_14April2020/LBC1936_DNAm_blood_14longitudinal_Danni_14April2020.rds")

## Reordering brain dataset so relevant columns with info occupy the first few column positions 
brain_meth <- brain_meth[c(1,2391:2405,2:2390)]

############################################################################################

### FORMAT METHYLATION DATA SO THAT IT CAN BE USED IN PREDICTOR CALCULATIONS

############################################################################################

### BRAIN - across regions as a mean DNAm value 

# Check to see the columns I need and the length of the dataset in the CPG direction
colnames(brain_meth)
ncol(brain_meth)

# Include only methylation data, with the LBC IDs to identify people 
brain <- brain_meth[c(15,17:2405)]
brain$LBC_ID <- as.character(brain$LBC_ID)

# Subset to 14 LBC IDs
brain_LBC <- left_join(IDs, brain, by = "LBC_ID")

## Transposing the data 

# keep the first column 
names <-  brain_LBC[,1]

# Transpose everything other than the first column
t_brain <- as.data.frame(as.matrix(t(brain_LBC[,-1])))

# Assign first column as the column names of the transposed dataframe
colnames(t_brain) <- names

# Find mean brain DNAm for each person across regions
t_brain$LBC360291_m <- apply(t_brain[,c(1:5)],1,mean)
t_brain$LBC360232_m <- apply(t_brain[,c(6:10)],1,mean)
t_brain$LBC360999_m <- apply(t_brain[,c(11:15)],1,mean)
t_brain$LBC360202_m <- apply(t_brain[,c(16:20)],1,mean)
t_brain$LBC361079_m <- apply(t_brain[,c(21:25)],1,mean)
t_brain$LBC360849_m <- apply(t_brain[,c(26:30)],1,mean)
t_brain$LBC360021_m <- apply(t_brain[,c(31:35)],1,mean)
t_brain$LBC360873_m <- apply(t_brain[,c(36:40)],1,mean)
t_brain$LBC360666_m <- apply(t_brain[,c(41:44)],1,mean)
t_brain$LBC360705_m <- apply(t_brain[,c(45:49)],1,mean)
t_brain$LBC361021_m <- apply(t_brain[,c(50:54)],1,mean)
t_brain$LBC360620_m <- apply(t_brain[,c(55:59)],1,mean)
t_brain$LBC360993_m <- apply(t_brain[,c(60:64)],1,mean)
t_brain$LBC360838_m <- apply(t_brain[,c(65:69)],1,mean)

# Select mean brain DNAm only
brain_mean <- t_brain[c(70:83)]

# Take the _m indicator off as we know they are all means now 
names(brain_mean) <- substring(names(brain_mean),1,9)


###########################################################################################

### BRAIN - for each region specifically (split into 5 datasets)

###########################################################################################

# Subsetting, but this time for each region (no phenotypes)

# Include only methylation data, with the LBC IDs and regions
regions <- brain_meth[c(15,16,17:2405)]
regions$LBC_ID <- as.character(regions$LBC_ID)

# Subset to my 14 LBC peoples IDs
regions <- left_join(IDs, regions, by = "LBC_ID")

# Create versions for each of the regions in the brain

BA46 <- subset(regions, region == "BA46")
BA46 <- BA46[-2]

BA17 <- subset(regions, region == "BA17")
BA17 <- BA17[-2]

BA2021 <- subset(regions, region == "BA20/21")
BA2021 <- BA2021[-2]

HC <- subset(regions, region == "HC")
HC <- HC[-2]

BA24 <- subset(regions, region == "BA24")
BA24 <- BA24[-2]

## Try transposing these regional datasets (as above)

names <-  BA46[,1]
t_BA46 <- as.data.frame(as.matrix(t(BA46[,-1])))
colnames(t_BA46) <- names

names1 <-  BA17[,1]
t_BA17 <- as.data.frame(as.matrix(t(BA17[,-1])))
colnames(t_BA17) <- names1

names2 <-  BA2021[,1]
t_BA2021 <- as.data.frame(as.matrix(t(BA2021[,-1])))
colnames(t_BA2021) <- names2

names3 <-  HC[,1]
t_HC <- as.data.frame(as.matrix(t(HC[,-1])))
colnames(t_HC) <- names3 # LBC360666 is not present in HC (only missing sample)

names4 <-  BA24[,1]
t_BA24 <- as.data.frame(as.matrix(t(BA24[,-1])))
colnames(t_BA24) <- names4


###########################################################################################

### BLOOD - for the most recent wave measurement taken prior to death 

###########################################################################################

# Check to see the columns I need and the length of the dataset in the CPG direction
colnames(blood_meth)
ncol(blood_meth)

# Include only methylation data, with the LBC IDs to identify people, with wave included
blood <- blood_meth[c(5,4,20:2413)]
names(blood)[1] <- "LBC_ID"
blood$LBC_ID <- as.character(blood$LBC_ID)

# Write this file out as a sense check 
write.csv(blood, "/Volumes/marioni-lab/Danni/Blood_brain_project/Oct_2020_brain_comms/Additional_checks_collated/blood_DNAm_data_with_all_waves_available.csv", row.names = F)


## Get only the most recent wave for each person
blood_recent <- blood %>%
	group_by(LBC_ID) %>%
	slice(which.max(WAVE))

# Make table into a dataframe format 
blood_recent <- as.data.frame(blood_recent)

# Check to see if I've got any NA values, or whether I have the most recent wave and its complete 
table(is.na(blood_recent))

# write this as a record too 
write.csv(blood_recent, "/Volumes/marioni-lab/Danni/Blood_brain_project/Oct_2020_brain_comms/Additional_checks_collated/blood_DNAm_data_with_all_waves_available_filtered_to_most_recent_wave.csv", row.names = F)


# Remove the wave column, as i have just taken the most recent measurements for each person only 
blood_recent <- blood_recent[-2]

# Transpose so that LBC_IDs are the column names 

# keep the first column 
names <-  blood_recent[,1]

# Transpose everything other than the first column
t_blood <- as.data.frame(as.matrix(t(blood_recent[,-1])))

# Assign first column as the column names of the transposed dataframe
colnames(t_blood) <- names




###########################################################################################

### SAVING FILES  

###########################################################################################

# Add NA values for missing persons data in HC region
t_HC$LBC360666 <- NA

# Make sure LBC_ID columns are in the same order 
t_HC <- t_HC[match(colnames(t_BA17), colnames(t_HC))]
t_blood <- t_blood[match(colnames(t_BA17), colnames(t_blood))]
brain_mean <- brain_mean[match(colnames(t_BA17), colnames(brain_mean))]
t_BA46 <- t_BA46[match(colnames(t_BA17), colnames(t_BA46))]
t_BA24 <- t_BA24[match(colnames(t_BA17), colnames(t_BA24))]
t_BA2021 <- t_BA2021[match(colnames(t_BA17), colnames(t_BA2021))]


# Save files
write.csv(t_blood, "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_8_epigenetic_signature_data_processing_outputs/List_regions/blood_methylation_most_recent_only_LBC_IDS_as_colnames.csv")

write.csv(brain_mean, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_8_epigenetic_signature_data_processing_outputs/List_regions/brain_mean_methylation_LBC_IDS_as_colnames.csv")

write.csv(t_BA46, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_8_epigenetic_signature_data_processing_outputs/List_regions/BA46_methylation_LBC_IDS_as_colnames.csv")

write.csv(t_BA17, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_8_epigenetic_signature_data_processing_outputs/List_regions/BA17_methylation_LBC_IDS_as_colnames.csv")

write.csv(t_BA2021, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_8_epigenetic_signature_data_processing_outputs/List_regions/BA2021_methylation_LBC_IDS_as_colnames.csv")

write.csv(t_HC, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_8_epigenetic_signature_data_processing_outputs/List_regions/HC_methylation_LBC_IDS_as_colnames.csv")

write.csv(t_BA24, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_8_epigenetic_signature_data_processing_outputs/List_regions/BA24_methylation_LBC_IDS_as_colnames.csv")


############################################################################################

### SESSION INFO

############################################################################################

# > sessionInfo()
# R version 3.6.3 (2020-02-29)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Catalina 10.15.4

# Matrix products: default
# BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

# locale:
# [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
#  [1] forcats_0.5.0   stringr_1.4.0   dplyr_0.8.5     purrr_0.3.3    
#  [5] readr_1.3.1     tidyr_1.0.3     ggplot2_3.3.0   tidyverse_1.3.0
#  [9] tibble_3.0.1    plyr_1.8.4      readxl_1.3.1   

# loaded via a namespace (and not attached):
#  [1] Rcpp_1.0.3       cellranger_1.1.0 pillar_1.4.4     compiler_3.6.3  
#  [5] dbplyr_1.4.2     tools_3.6.3      jsonlite_1.6     lubridate_1.7.4 
#  [9] lifecycle_0.2.0  nlme_3.1-144     gtable_0.3.0     lattice_0.20-38 
# [13] pkgconfig_2.0.3  rlang_0.4.6      reprex_0.3.0     cli_1.1.0       
# [17] rstudioapi_0.10  DBI_1.0.0        haven_2.2.0      withr_2.1.2     
# [21] xml2_1.2.2       httr_1.4.1       fs_1.4.1         generics_0.0.2  
# [25] vctrs_0.2.4      hms_0.5.2        grid_3.6.3       tidyselect_1.0.0
# [29] glue_1.4.0       R6_2.4.1         fansi_0.4.0      modelr_0.1.5    
# [33] magrittr_1.5     backports_1.1.5  scales_1.1.0     ellipsis_0.3.0  
# [37] rvest_0.3.5      assertthat_0.2.1 colorspace_1.4-1 utf8_1.1.4      
# [41] stringi_1.4.3    munsell_0.5.0    broom_0.5.5      crayon_1.3.4    



############################################################################################


















