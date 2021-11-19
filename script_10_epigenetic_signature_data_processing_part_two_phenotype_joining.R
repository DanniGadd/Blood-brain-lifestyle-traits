############################################################################################
############################################################################################
################# SCRIPT 10 - Danni Gadd - Blood vs brain DNAm project ######################
############################################################################################
############################################################################################

setwd("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts")
library(tidyverse)

############################################################################################

### LOAD DATA

############################################################################################

# Predictor scores 
smoking <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_9_calculation_of_epigenetic_predictor_signatures_for_traits_outputs/smoking_predictors.csv")
alcohol <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_9_calculation_of_epigenetic_predictor_signatures_for_traits_outputs/alcohol_predictors.csv")
BMI <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_9_calculation_of_epigenetic_predictor_signatures_for_traits_outputs/BMI_predictors.csv")
HDL <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_9_calculation_of_epigenetic_predictor_signatures_for_traits_outputs/HDL_predictors.csv")

# Save these predictor scores as a list too, so they can be read into the functions below
files = list.files("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_9_calculation_of_epigenetic_predictor_signatures_for_traits_outputs", pattern="*.csv", full.names = TRUE)

# Data used in script 5 for previous regressions on the smoking CPG 
data_recent <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_6_smoking_CPG_data_processing_part_two_phenotype_joining_for_regressions_outputs/smoking_CPG_data_for_statistical_testing_most_recent_blood_data_only.csv")

# Load in trait info for the epiegentic predictor relevant ones generated earlier 
alcohol_recent <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_5_extract_phenotype_most_recent_measure_before_death_outputs/alcohol_info_most_recent.csv")
alcohol_recent <- alcohol_recent[-1]

BMI_recent <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_5_extract_phenotype_most_recent_measure_before_death_outputs/BMI_info_most_recent.csv")
BMI_recent <- BMI_recent[-1]

HDL_recent <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_5_extract_phenotype_most_recent_measure_before_death_outputs/HDL_info_most_recent.csv")
HDL_recent <- HDL_recent[-1]


###########################################################################################

### Join new predictors up to the main dataframe for testing 

###########################################################################################

## Transposing the predictor data 

# Use lapply to do the same steps for predictor score transposing for each file 
dfList <- lapply(files, function(i) {
     df <- read.csv(i) # read in predictor score file 
     df <- df[-1] # get rid of first column as its just junk 
     names <- df[,1] # get the names of the samples before transposing so that these match up 
     t_df <- as.data.frame(as.matrix(t(df[,-1]))) # transpose everything other than the first column
     colnames(t_df) <- names # assign first column as the names of the transposed dataframe 
     t_df$LBC_ID <- rownames(t_df) # Create LBC_ID column for merging
     t_df <- t_df[c(8,1:7)] # Reorder so the LBC_ID column is first 
     return(t_df) # output the transformed datasets for the predictor scores
})

# Extract the ones I want as separate dataframes (double checked based on order of files list above)
alcohol <- as.data.frame(dfList[1])
BMI <- as.data.frame(dfList[2])
HDL <- as.data.frame(dfList[3])
smoking <- as.data.frame(dfList[4])

# Add identifiers for each trait on datasets 
names(alcohol) <- paste0(names(alcohol), "_ALC")
names(BMI) <- paste0(names(BMI), "_BMI")
names(smoking) <- paste0(names(smoking), "_SMO")
names(HDL) <- paste0(names(HDL), "_HDL")

# Change LBC back for merging
names(alcohol)[1] <- "LBC_ID"
names(BMI)[1] <- "LBC_ID"
names(smoking)[1] <- "LBC_ID"
names(HDL)[1] <- "LBC_ID"

#### Create a version with gathered data for each trait in brain regions

### HDL 
regions <- HDL[c(-6,-7)]
names(regions) <- c("LBC_ID", "BA17", "BA20/21", "BA24", "BA46", "HC")

gather <- gather(regions, key= LBC_ID, value=region)
names(gather)[1] <- "region"
names(gather)[2] <- "HDL_score_all_brain"

# Get the IDs again
ID <- HDL$LBC_ID 

# Add the IDs abck in again in the right order 
gather$LBC_ID <- rep(ID, 5)

# Remove the missing value at position 65
gather <- gather[-65,]
HDL_collated <- gather

### BMI
regions <- BMI[c(-6,-7)]
names(regions) <- c("LBC_ID", "BA17", "BA20/21", "BA24", "BA46", "HC")

gather <- gather(regions, key= LBC_ID, value=region)
names(gather)[1] <- "region"
names(gather)[2] <- "BMI_score_all_brain"
ID <- BMI$LBC_ID 
gather$LBC_ID <- rep(ID, 5)

# Remove the missing value at position 65
gather <- gather[-65,]
BMI_collated <- gather

### ALC
regions <- alcohol[c(-6,-7)]
names(regions) <- c("LBC_ID", "BA17", "BA20/21", "BA24", "BA46", "HC")

gather <- gather(regions, key= LBC_ID, value=region)
names(gather)[1] <- "region"
names(gather)[2] <- "ALC_score_all_brain"
ID <- alcohol$LBC_ID 
gather$LBC_ID <- rep(ID, 5)

# Remove the missing value at position 65
gather <- gather[-65,]
ALC_collated <- gather

### SMOK
regions <- smoking[c(-6,-7)]
names(regions) <- c("LBC_ID", "BA17", "BA20/21", "BA24", "BA46", "HC")

gather <- gather(regions, key= LBC_ID, value=region)
names(gather)[1] <- "region"
names(gather)[2] <- "SMOK_score_all_brain"
ID <- smoking$LBC_ID 
gather$LBC_ID <- rep(ID, 5)

# Remove the missing value at position 65
gather <- gather[-65,]
SMOK_collated <- gather


#### ADD PREDICTOR SCORE VARIABLES INTO MAIN DATASET

# See if you can merge together with main dataset 
# the individual regions as columns and the combined cross-region column for each trait

# Individual scores by-sample as variables 
data <- merge(data_recent, smoking, by = "LBC_ID")
data <- merge(data, alcohol, by = "LBC_ID")
data <- merge(data, BMI, by = "LBC_ID")
data <- merge(data, HDL, by = "LBC_ID")

# Change factors to characters for joining 
data$LBC_ID <- as.character(data$LBC_ID)
data$region <- as.character(data$region)

# Collated scores that have each region of the brain under a trait variable 
data <- left_join(data, SMOK_collated, by = c("LBC_ID" = "LBC_ID", "region" = "region"))
data <- left_join(data, ALC_collated, by = c("LBC_ID" = "LBC_ID", "region" = "region"))
data <- left_join(data, HDL_collated, by = c("LBC_ID" = "LBC_ID", "region" = "region"))
data <- left_join(data, BMI_collated, by = c("LBC_ID" = "LBC_ID", "region" = "region"))


### ADD THE EXTRA PHENOTYPES TO DATA FRAME FOR HDL, BMI AND ALCOHOL
data <- merge(data, alcohol_recent, by = "LBC_ID")
data <- merge(data, BMI_recent, by = "LBC_ID")
data <- merge(data, HDL_recent, by = "LBC_ID")

# Tidy up dataset as a final version to use in regression analyses

# Remove junk columns
data <- data[-2]
data <- data[-5]

# Rename smoking CpG site brain and blood beta DNAm values 
names(data)[8] <- "cg05575921_brain" 
names(data)[11] <- "cg05575921_blood" 

# Rename age as its years not days
names(data)[19] <- "age_yrs_death"

# Save data file for summary stats table production in the next script
write.csv(data, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_10_epigenetic_signature_data_processing_part_two_phenotype_joining_for_regressions_outputs/data_withmost_recent_clinical_phenotypes_added.csv")


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
# [1] forcats_0.5.0   stringr_1.4.0   dplyr_0.8.5     purrr_0.3.3    
# [5] readr_1.3.1     tidyr_1.0.3     tibble_3.0.1    ggplot2_3.3.0  
# [9] tidyverse_1.3.0

# loaded via a namespace (and not attached):
#  [1] Rcpp_1.0.3       cellranger_1.1.0 pillar_1.4.4     compiler_3.6.3  
#  [5] dbplyr_1.4.2     tools_3.6.3      jsonlite_1.6     lubridate_1.7.4 
#  [9] lifecycle_0.2.0  nlme_3.1-144     gtable_0.3.0     lattice_0.20-38 
# [13] pkgconfig_2.0.3  rlang_0.4.6      reprex_0.3.0     cli_1.1.0       
# [17] rstudioapi_0.10  DBI_1.0.0        haven_2.2.0      withr_2.1.2     
# [21] xml2_1.2.2       httr_1.4.1       fs_1.4.1         generics_0.0.2  
# [25] vctrs_0.2.4      hms_0.5.2        grid_3.6.3       tidyselect_1.0.0
# [29] glue_1.4.0       R6_2.4.1         readxl_1.3.1     modelr_0.1.5    
# [33] magrittr_1.5     backports_1.1.5  scales_1.1.0     ellipsis_0.3.0  
# [37] rvest_0.3.5      assertthat_0.2.1 colorspace_1.4-1 stringi_1.4.3   
# [41] munsell_0.5.0    broom_0.5.5      crayon_1.3.4    


############################################################################################












