############################################################################################
############################################################################################
################# SCRIPT 7 - Danni Gadd - Blood vs brain DNAm project ######################
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

# Assign predictor weights file path

xl_data <- "/Volumes/marioni-lab/Danni/Blood_brain_project/Data/GS_complex_trait_predictor_weights_from_genome_biology_4_May_2020/13059_2018_1514_MOESM1_ESM.xlsx"

# View sheet names

excel_sheets(path = xl_data)

#  [1] "Table S1 - BMI"                "Table S2 - Smoking"           
#  [3] "Table S3 - Alcohol"            "Table S4 - Education"         
#  [5] "Table S5  - Total cholesterol" "Table S6 - HDL cholesterol"   
#  [7] "Table S7 - LDL cholesterol"    "Table S8 - Total-HDL ratio"   
#  [9] "Table S9 - Waist-to-Hip ratio" "Table S10 - % Body fat" 

# Read in all tabs as a list of dataframes using lapply()

tab_name <- excel_sheets(path = xl_data)

list_all <- lapply(tab_name, function(x) read_excel(path = xl_data, sheet = x))

str(list_all)

############################################################################################

### SAVE WEIGHTS

############################################################################################

# Save weights for traits of interest as dataframes
BMI <- list_all[1] %>% as.data.frame
SMOK <- list_all[2] %>% as.data.frame
ALC <- list_all[3] %>% as.data.frame
EDU <- list_all[4] %>% as.data.frame
CHOL <- list_all[5] %>% as.data.frame
HDL <- list_all[6] %>% as.data.frame
LDL <- list_all[7] %>% as.data.frame
CHOL_RAT <- list_all[8] %>% as.data.frame
WAIST_HIP <- list_all[9] %>% as.data.frame
BODY_FAT <- list_all[10] %>% as.data.frame

# Output the files so that they can be easily read in and used in the predictor generation in the next script
write.csv(BMI, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_7_extract_predictor_weights_for_traits_outputs/BMI_weights.csv")
write.csv(SMOK, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_7_extract_predictor_weights_for_traits_outputs/SMOK_weights.csv")
write.csv(ALC, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_7_extract_predictor_weights_for_traits_outputs/ALC_weights.csv")
write.csv(EDU, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_7_extract_predictor_weights_for_traits_outputs/EDU_weights.csv")
write.csv(CHOL, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_7_extract_predictor_weights_for_traits_outputs/CHOL_weights.csv")
write.csv(HDL, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_7_extract_predictor_weights_for_traits_outputs/HDL_weights.csv")
write.csv(LDL, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_7_extract_predictor_weights_for_traits_outputs/LDL_weights.csv")
write.csv(CHOL_RAT, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_7_extract_predictor_weights_for_traits_outputs/CHOL_RAT_weights.csv")
write.csv(WAIST_HIP, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_7_extract_predictor_weights_for_traits_outputs/WAIST_HIP_weights.csv")
write.csv(BODY_FAT, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_7_extract_predictor_weights_for_traits_outputs/BODY_FAT_weights.csv")

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
# [29] glue_1.4.0       R6_2.4.1         modelr_0.1.5     magrittr_1.5    
# [33] backports_1.1.5  scales_1.1.0     ellipsis_0.3.0   rvest_0.3.5     
# [37] assertthat_0.2.1 colorspace_1.4-1 stringi_1.4.3    munsell_0.5.0   
# [41] broom_0.5.5      crayon_1.3.4    



############################################################################################


















