############################################################################################
############################################################################################
################# SCRIPT 3 - Danni Gadd - Blood vs brain DNAm project ######################
############################################################################################
############################################################################################

setwd("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts")
library(tidyverse)
library(dplyr)

############################################################################################

### LOAD DATA

############################################################################################

# Read in the wide format with DNAm at blood waves generated in the previous script 
bloodDNAm <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_2_smoking_CPG_data_processing_part_one_for_correlations_outputs/smoking_CPG_data_wide_format_by_wave_and_region.csv")

# LBC phenotypes file prepared in script 1
LBC_phenos <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_1_extracting_14_LBC_individuals_from_the_LBC_phenotype_dataset_outputs/LBC_phenotypes_dataset_14_individuals.csv")

############################################################################################

### Get the most recent blood DNAm and the corresponding age which it was taken

############################################################################################

# Assign a wave identifier column
bloodDNAm$Wave_id<-ifelse(bloodDNAm$X4 == "NA",1,
		ifelse(bloodDNAm$X4 >0,4,
		))

bloodDNAm$Wave_id[is.na(bloodDNAm$Wave_id)] <- 3

# For any values at wave 4 which were NA, replace with values at wave 3 so that each persons most recent blood DNAm is used 
bloodDNAm$X4[is.na(bloodDNAm$X4)] <- bloodDNAm$X3[is.na(bloodDNAm$X4)]

# Grab the X4 column which now contains the most recent blood DNAm data 
blood_DNAm_recent <- bloodDNAm[c(2,6,12)]

# Get age from phenotype data 
agedays_DNAm <- LBC_phenos[c(2,4:7)]

# For any values which are NA in wave 4, replace with age at wave 3 so that the age at most recent DNAm sample is reflected 
agedays_DNAm$agedays_w4[is.na(agedays_DNAm$agedays_w4)] <- agedays_DNAm$agedays_w3[is.na(agedays_DNAm$agedays_w4)]

# Grab the agedaysW4 column which now represents the age at which the most recent blood DNAm was taken prior to death
age_blood_DNAm_recent <- agedays_DNAm[c(1,5)]

# Join these together so I have DNAm at last time point before death and age at last time point before death 
blood_DNAm_recent <- merge(blood_DNAm_recent, age_blood_DNAm_recent, by = "LBC_ID")

# Rename so that the variables reflect what they actually are now
names(blood_DNAm_recent)[2] <- "DNAm_blood_recent"
names(blood_DNAm_recent)[4] <- "age_blood_recent"

# Now convert so that age values are in years not days 
blood_DNAm_recent$age_blood_recent <- blood_DNAm_recent$age_blood_recent / 365.25

# Save this as its quite a handy thing to have reference to without doing these steps 
write.csv(blood_DNAm_recent, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_3_extract_most_recent_blood_wave_and_age_at_which_blood_was_taken_outputs/blood_DNAm_recent_with_age_at_recent.csv")


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
# [1] forcats_0.5.0   stringr_1.4.0   dplyr_0.8.5     purrr_0.3.3    
# [5] readr_1.3.1     tidyr_1.0.3     tibble_3.0.1    ggplot2_3.3.0  
# [9] tidyverse_1.3.0

# loaded via a namespace (and not attached):
#  [1] Rcpp_1.0.3       cellranger_1.1.0 pillar_1.4.4     compiler_3.6.3  
#  [5] dbplyr_1.4.2     tools_3.6.3      jsonlite_1.6     lubridate_1.7.4 
#  [9] lifecycle_0.2.0  nlme_3.1-144     gtable_0.3.0     lattice_0.20-41 
# [13] pkgconfig_2.0.3  rlang_0.4.6      reprex_0.3.0     cli_1.1.0       
# [17] rstudioapi_0.10  DBI_1.0.0        haven_2.2.0      withr_2.1.2     
# [21] xml2_1.2.2       httr_1.4.1       fs_1.4.1         generics_0.0.2  
# [25] vctrs_0.2.4      hms_0.5.2        grid_3.6.3       tidyselect_1.0.0
# [29] glue_1.4.0       R6_2.4.1         readxl_1.3.1     modelr_0.1.5    
# [33] magrittr_1.5     backports_1.1.5  scales_1.1.0     ellipsis_0.3.0  
# [37] rvest_0.3.5      assertthat_0.2.1 colorspace_1.4-1 stringi_1.4.3   
# [41] munsell_0.5.0    broom_0.5.5      crayon_1.3.4     



############################################################################################


















