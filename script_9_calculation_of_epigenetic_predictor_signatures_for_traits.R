############################################################################################
############################################################################################
################# SCRIPT 9 - Danni Gadd - Blood vs brain DNAm project ######################
############################################################################################
############################################################################################

setwd("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts")
library(tidyverse)
library(lubridate)
library(dplyr)

############################################################################################

### LOAD DATA

############################################################################################

# Save sample files as a list to read into functions below
files = list.files("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_8_epigenetic_signature_data_processing_outputs/List_regions", pattern="*.csv", full.names = TRUE)

# Print files list for ordering
# > files
# [1] "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_8_epigenetic_signature_data_processing_outputs/List_regions/BA17_methylation_LBC_IDS_as_colnames.csv"                  
# [2] "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_8_epigenetic_signature_data_processing_outputs/List_regions/BA2021_methylation_LBC_IDS_as_colnames.csv"                
# [3] "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_8_epigenetic_signature_data_processing_outputs/List_regions/BA24_methylation_LBC_IDS_as_colnames.csv"                  
# [4] "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_8_epigenetic_signature_data_processing_outputs/List_regions/BA46_methylation_LBC_IDS_as_colnames.csv"                  
# [5] "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_8_epigenetic_signature_data_processing_outputs/List_regions/blood_methylation_most_recent_only_LBC_IDS_as_colnames.csv"
# [6] "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_8_epigenetic_signature_data_processing_outputs/List_regions/brain_mean_methylation_LBC_IDS_as_colnames.csv"            
# [7] "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_8_epigenetic_signature_data_processing_outputs/List_regions/HC_methylation_LBC_IDS_as_colnames.csv"

# Make sure this order is matched in the folder with the listed files
file_order <- c("BA17","BA2021", "BA24", "BA46", "blood", "brain", "HC")

###########################################################################################

## SMOKING PREDICTOR SCORES

###########################################################################################

## Smoking predictor calculation for each sample (mean brain DNAm, blood DNAm and 5 brain regions DNAm separately)
## The check for this (the way i first did it manually) are at the very end of this script. I wanted to keep them
## as they are a useful sense check to make sure the lapply function is doing what it is supposed to.

# Use lapply to do the same steps for predictor score generation for each file 
dfList <- lapply(files, function(i) {
     df <- read.csv(i) # read in sample DNAm data file
     weight <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_7_extract_predictor_weights_for_traits_outputs/SMOK_weights.csv") # read in smoking weights
     df <- df %>% filter(X %in% weight$CpG) # filter the sample dataset by the CpG sites in the predictor weight file
     not_names <- setdiff(weight$CpG, df$X) # find out if there are any CpGs in the weights file that arent in the sample dataset
     weight_re <- weight[weight$CpG %in% df$X,] # remove the CpGs which dont match up from the predictor weights 
     re_sm_match <- df[match(weight_re$CpG,df$X),] # make sure the order of CpGs now matches between the files
     re_weights <- re_sm_match[,c(2:15)]*weight_re[,3] # multiply the methylation values with the predicto weights values 
     re_sum <- colSums(re_weights) # sum the columns to get a predictor score for each individual
     return(re_sum) # output scores
})

smoking_predictors <- do.call(rbind, dfList) # bind scores together so that each sample dataset has one row in the resulting table 
smoking_predictors <- as.data.frame(smoking_predictors)
smoking_predictors$files <- file_order
smoking_predictors <- smoking_predictors[c(15,1:14)]


###########################################################################################

## ALCOHOL PREDICTOR SCORES

###########################################################################################

# Use lapply to do the same steps for predictor score generation for each file 
dfList <- lapply(files, function(i) {
     df <- read.csv(i) # read in sample DNAm data file
     weight <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_7_extract_predictor_weights_for_traits_outputs/ALC_weights.csv") # read in alcohol weights
     df <- df %>% filter(X %in% weight$CpG) # filter the sample dataset by the CpG sites in the predictor weight file
     not_names <- setdiff(weight$CpG, df$X) # find out if there are any CpGs in the weights file that arent in the sample dataset
     weight_re <- weight[weight$CpG %in% df$X,] # remove the CpGs which dont match up from the predictor weights 
     re_sm_match <- df[match(weight_re$CpG,df$X),] # make sure the order of CpGs now matches between the files
     re_weights <- re_sm_match[,c(2:15)]*weight_re[,3] # multiply the methylation values with the predicto weights values 
     re_sum <- colSums(re_weights) # sum the columns to get a predictor score for each individual
     return(re_sum) # output scores
})

alcohol_predictors <- do.call(rbind, dfList) # bind scores together so that each sample dataset has one row in the resulting table 
alcohol_predictors <- as.data.frame(alcohol_predictors)
alcohol_predictors$files <- file_order
alcohol_predictors <- alcohol_predictors[c(15,1:14)]

###########################################################################################

## BMI PREDICTOR SCORES

###########################################################################################

# Use lapply to do the same steps for predictor score generation for each file 
dfList <- lapply(files, function(i) {
     df <- read.csv(i) # read in sample DNAm data file
     weight <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_7_extract_predictor_weights_for_traits_outputs/BMI_weights.csv") # read in alcohol weights
     df <- df %>% filter(X %in% weight$CpG) # filter the sample dataset by the CpG sites in the predictor weight file
     not_names <- setdiff(weight$CpG, df$X) # find out if there are any CpGs in the weights file that arent in the sample dataset
     weight_re <- weight[weight$CpG %in% df$X,] # remove the CpGs which dont match up from the predictor weights 
     re_sm_match <- df[match(weight_re$CpG,df$X),] # make sure the order of CpGs now matches between the files
     re_weights <- re_sm_match[,c(2:15)]*weight_re[,3] # multiply the methylation values with the predicto weights values 
     re_sum <- colSums(re_weights) # sum the columns to get a predictor score for each individual
     return(re_sum) # output scores
})

BMI_predictors <- do.call(rbind, dfList) # bind scores together so that each sample dataset has one row in the resulting table 
BMI_predictors <- as.data.frame(BMI_predictors)
BMI_predictors$files <- file_order
BMI_predictors <- BMI_predictors[c(15,1:14)]

###########################################################################################

## HDL PREDICTOR SCORES

###########################################################################################

# Use lapply to do the same steps for predictor score generation for each file 
dfList <- lapply(files, function(i) {
     df <- read.csv(i) # read in sample DNAm data file
     weight <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_7_extract_predictor_weights_for_traits_outputs/HDL_weights.csv") # read in alcohol weights
     df <- df %>% filter(X %in% weight$CpG) # filter the sample dataset by the CpG sites in the predictor weight file
     not_names <- setdiff(weight$CpG, df$X) # find out if there are any CpGs in the weights file that arent in the sample dataset
     weight_re <- weight[weight$CpG %in% df$X,] # remove the CpGs which dont match up from the predictor weights 
     re_sm_match <- df[match(weight_re$CpG,df$X),] # make sure the order of CpGs now matches between the files
     re_weights <- re_sm_match[,c(2:15)]*weight_re[,3] # multiply the methylation values with the predicto weights values 
     re_sum <- colSums(re_weights) # sum the columns to get a predictor score for each individual
     return(re_sum) # output scores
})

HDL_predictors <- do.call(rbind, dfList) # bind scores together so that each sample dataset has one row in the resulting table 
HDL_predictors <- as.data.frame(HDL_predictors)
HDL_predictors$files <- file_order
HDL_predictors <- HDL_predictors[c(15,1:14)]



###########################################################################################

### SAVE PREDICTOR SCORE FILES 

###########################################################################################


write.csv(smoking_predictors, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_9_calculation_of_epigenetic_predictor_signatures_for_traits_outputs/smoking_predictors.csv")
write.csv(alcohol_predictors, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_9_calculation_of_epigenetic_predictor_signatures_for_traits_outputs/alcohol_predictors.csv")
write.csv(BMI_predictors, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_9_calculation_of_epigenetic_predictor_signatures_for_traits_outputs/BMI_predictors.csv")
write.csv(HDL_predictors, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_9_calculation_of_epigenetic_predictor_signatures_for_traits_outputs/HDL_predictors.csv")


###########################################################################################

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
#  [1] lubridate_1.7.4 forcats_0.5.0   stringr_1.4.0   dplyr_0.8.5    
#  [5] purrr_0.3.3     readr_1.3.1     tidyr_1.0.3     tibble_3.0.1   
#  [9] ggplot2_3.3.0   tidyverse_1.3.0

# loaded via a namespace (and not attached):
#  [1] Rcpp_1.0.3       cellranger_1.1.0 pillar_1.4.4     compiler_3.6.3  
#  [5] dbplyr_1.4.2     tools_3.6.3      jsonlite_1.6     lifecycle_0.2.0 
#  [9] nlme_3.1-144     gtable_0.3.0     lattice_0.20-38  pkgconfig_2.0.3 
# [13] rlang_0.4.6      reprex_0.3.0     cli_1.1.0        rstudioapi_0.10 
# [17] DBI_1.0.0        haven_2.2.0      withr_2.1.2      xml2_1.2.2      
# [21] httr_1.4.1       fs_1.4.1         generics_0.0.2   vctrs_0.2.4     
# [25] hms_0.5.2        grid_3.6.3       tidyselect_1.0.0 glue_1.4.0      
# [29] R6_2.4.1         readxl_1.3.1     modelr_0.1.5     magrittr_1.5    
# [33] backports_1.1.5  scales_1.1.0     ellipsis_0.3.0   rvest_0.3.5     
# [37] assertthat_0.2.1 colorspace_1.4-1 stringi_1.4.3    munsell_0.5.0   
# [41] broom_0.5.5      crayon_1.3.4  


############################################################################################

### END OF MAIN SCRIPT

###########################################################################################

### MANUAL CHECK TO MAKE SURE LAPPLY FUNCTIONS ABOVE ARE WORKING CORRECTLY

# Samples 
blood <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_8_epigenetic_signature_data_processing_outputs/List_regions/blood_methylation_most_recent_only_LBC_IDS_as_colnames.csv")
brain <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_8_epigenetic_signature_data_processing_outputs/List_regions/brain_mean_methylation_LBC_IDS_as_colnames.csv")
BA17 <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_8_epigenetic_signature_data_processing_outputs/List_regions/BA17_methylation_LBC_IDS_as_colnames.csv")
BA24 <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_8_epigenetic_signature_data_processing_outputs/List_regions/BA24_methylation_LBC_IDS_as_colnames.csv")
BA46 <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_8_epigenetic_signature_data_processing_outputs/List_regions/BA46_methylation_LBC_IDS_as_colnames.csv")
BA2021 <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_8_epigenetic_signature_data_processing_outputs/List_regions/BA2021_methylation_LBC_IDS_as_colnames.csv")
HC <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_8_epigenetic_signature_data_processing_outputs/List_regions/HC_methylation_LBC_IDS_as_colnames.csv")

# Predictor weights 
SMOK <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_7_extract_predictor_weights_for_traits_outputs/SMOK_weights.csv")
ALC <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_7_extract_predictor_weights_for_traits_outputs/ALC_weights.csv")
BMI <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_7_extract_predictor_weights_for_traits_outputs/BMI_weights.csv")


###########################################################################################

## BRAIN - SMOKING PREDICTOR LAPPLY FUNCTION MANUAL CHECK 

###########################################################################################

## Subset to CPGs for specific trait 
brain_sm <- brain %>%
	filter(X %in% SMOK$CpG)

## Find CPG sites which are in predictor list but not in brain list  

# Check to see where the values are that dont match up (8 values)
not <- !(SMOK$CpG %in% brain_sm$X) 

# Get names of the non-matching culprits
not_names <- setdiff(SMOK$CpG, brain_sm$X)

# Remove the non-matching names of the CPGs from the smoking weights dataset 
SMOK_br <- SMOK[SMOK$CpG %in% brain_sm$X,]

# Make sure the order of the CPGs in my file is matched to the order in the predictor weights file
brain_sm_match <- brain_sm[match(SMOK_br$CpG,brain_sm$X),]

# Convert all variables to characters/numerics 
brain_sm_match <- brain_sm_match %>% mutate_if(is.integer, as.numeric) -> brain_sm_match
brain_sm_match <- brain_sm_match %>% mutate_if(is.factor, as.character) -> brain_sm_match
SMOK_br <- SMOK_br %>% mutate_if(is.integer, as.numeric) -> SMOK_br
SMOK_br <- SMOK_br %>% mutate_if(is.factor, as.character) -> SMOK_br

# Check these now are identically matched 
identical(brain_sm_match$X, SMOK_br$CpG)

# Multiply the predictor weights by the CPGs
brain_weights <- brain_sm_match[,c(2:15)]*SMOK_br[,3]

# Sum the CPG values to get a predictor value for each person for this trait as a new column
brain_sum <- colSums(brain_weights)


###########################################################################################

## BLOOD - SMOKING PREDICTOR LAPPLY FUNCTION MANUAL CHECK 

###########################################################################################

# Follow the process as above for brain but with blood 
blood_sm <- blood %>%
	filter(X %in% SMOK$CpG)

not_names <- setdiff(SMOK$CpG, blood_sm$X)

SMOK_bl <- SMOK[SMOK$CpG %in% blood_sm$X,]

blood_sm_match <- blood_sm[match(SMOK_bl$CpG,blood_sm$X),]

blood_sm_match <- blood_sm_match %>% mutate_if(is.integer, as.numeric) -> blood_sm_match
blood_sm_match <- blood_sm_match %>% mutate_if(is.factor, as.character) -> blood_sm_match
SMOK_bl <- SMOK_bl %>% mutate_if(is.integer, as.numeric) -> SMOK_bl
SMOK_bl <- SMOK_bl %>% mutate_if(is.factor, as.character) -> SMOK_bl

identical(blood_sm_match$X, SMOK_bl$CpG)

blood_weights <- blood_sm_match[,c(2:15)]*SMOK_bl[,3]

blood_sum <- colSums(blood_weights)


############################################################################################
















