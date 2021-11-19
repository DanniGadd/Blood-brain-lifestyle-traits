############################################################################################
############################################################################################
################# SCRIPT 4 - Danni Gadd - Blood vs brain DNAm project ######################
############################################################################################
############################################################################################

setwd("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts")
library(tidyverse)
library(dplyr)
library(haven)

############################################################################################

### LOAD DATA

############################################################################################

# Read in the .sav file which has the age at death data for the LBC total cohort 
path = file.path("/Volumes/marioni-lab/Danni/Blood_brain_project/Data/LBC_age_at_death_data_complete_27_April_2020/LBC1936_Blood_DNAm_And_Brain_DNAm_RM_26APR2020.sav")
agedeath = read_sav(path)
agedeath <- as.data.frame(agedeath)

## LBC 14 individuals subject IDs
IDs <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/Data/post_mortem_DNAm_extraction_Anna_14April2020/LBC_IDs.csv")

############################################################################################

### Sorting out age at death data for inclusion 

############################################################################################

# Get the IDs I want to extract for 
LBC_ID <- IDs[1]
names(LBC_ID)[1] <- "ID"

# now try to merge based on the 14 IDs
brain_data_LBC <- merge(LBC_ID, agedeath, by.x = "ID", by.y= "lbc36no")

# Get age related data 
brain_data_LBC <- brain_data_LBC[c(1:8)]

# Grab just the info I need 
df <- brain_data_LBC[c(1,8)]
names(df)[1] <- "LBC_ID"
names(df)[2] <- "agedays_death_new"

# Save as a reference csv file 
write.csv(df, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_4_extract_age_at_death_for_inclusion_using_updated_death_info_outputs/age_days_at_death_data_LBC_updated_deaths.csv")


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
#  [1] haven_2.2.0     forcats_0.5.0   stringr_1.4.0   dplyr_0.8.5    
#  [5] purrr_0.3.3     readr_1.3.1     tidyr_1.0.3     tibble_3.0.1   
#  [9] ggplot2_3.3.0   tidyverse_1.3.0

# loaded via a namespace (and not attached):
#  [1] Rcpp_1.0.3       cellranger_1.1.0 pillar_1.4.4     compiler_3.6.3  
#  [5] dbplyr_1.4.2     tools_3.6.3      jsonlite_1.6     lubridate_1.7.4 
#  [9] lifecycle_0.2.0  nlme_3.1-144     gtable_0.3.0     lattice_0.20-41 
# [13] pkgconfig_2.0.3  rlang_0.4.6      reprex_0.3.0     cli_1.1.0       
# [17] rstudioapi_0.10  DBI_1.0.0        withr_2.1.2      xml2_1.2.2      
# [21] httr_1.4.1       fs_1.4.1         generics_0.0.2   vctrs_0.2.4     
# [25] hms_0.5.2        grid_3.6.3       tidyselect_1.0.0 glue_1.4.0      
# [29] R6_2.4.1         readxl_1.3.1     modelr_0.1.5     magrittr_1.5    
# [33] backports_1.1.5  scales_1.1.0     ellipsis_0.3.0   rvest_0.3.5     
# [37] assertthat_0.2.1 colorspace_1.4-1 stringi_1.4.3    munsell_0.5.0   
# [41] broom_0.5.5      crayon_1.3.4 



############################################################################################


















