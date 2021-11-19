############################################################################################
############################################################################################
################# SCRIPT 2 - Danni Gadd - Blood vs brain DNAm project ######################
############################################################################################
############################################################################################

setwd("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts")
library(tidyverse)
library(dplyr)

############################################################################################

### LOAD DATA 

############################################################################################

## LBC 14 individuals subject IDs
IDs <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/Data/post_mortem_DNAm_extraction_Anna_14April2020/LBC_IDs.csv")

## Brain and blood methylation datasets 
brain_meth <- readRDS("/Volumes/marioni-lab/Danni/Blood_brain_project/Data/post_mortem_DNAm_extraction_Anna_14April2020/brain_methylation_data.rds")
blood_meth <- readRDS("/Volumes/marioni-lab/Danni/Blood_brain_project/Data/pre_mortem_DNAm_extraction_REM_14April2020/LBC1936_DNAm_blood_14longitudinal_Danni_14April2020.rds")

## Reordering brain dataset so relevant columns with info occupy the first few column positions 
brain_meth <- brain_meth[c(1,2391:2405,2:2390)]

############################################################################################

### PREP DATA FOR CORRELATIONS 

############################################################################################

# Required data 
brain <- brain_meth[c(-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12,-13,-14)]
blood <- blood_meth[c(-2,-3,-6,-7,-8,-9,-10,-11,-12,-13,-14,-15,-16,-17,-18,-19)]

# Reorder 
names(blood)[3] <- "LBC_ID"
blood <- blood[c(3,1,2,4:2392)]
brain <- brain[c(2,1,3,4:2392)]

# Name samples and identifiers 
names(blood)[2] <- "Sample_ID"
names(brain)[2] <- "Sample_ID"
names(blood)[3] <- "Identifier"
names(brain)[3] <- "Identifier"

# Change to character variables
blood$LBC_ID <- as.character(blood$LBC_ID)
blood$Sample_ID <- as.character(blood$Sample_ID)
blood$Identifier <- as.character(blood$Identifier)
brain$LBC_ID <- as.character(brain$LBC_ID)
brain$Sample_ID <- as.character(brain$Sample_ID)
brain$Identifier <- as.character(brain$Identifier)

# Index position of CPG of interest 
D <- which(colnames(blood) == "cg05575921") # 353
E <- which(colnames(brain) == "cg05575921") # 195

# Subsetting to the smoking CpG of interest 
blood_test <- blood[c(1,2,3,353)]
brain_test <- brain[c(1,2,3,195)]

# Dont need sample_ID
blood_s <- blood_test[-2]
brain_s <- brain_test[-2]

# Spreading so that IDs get brought up as columns
spread <- blood_s %>%
	spread(key = Identifier, value = cg05575921)

spread2 <- brain_s %>%
	spread(key = Identifier, value = cg05575921)

# Check to see if they do actually match up 
identical(rownames(spread), rownames(spread2)) 

# Binding them together (the first wide format)
merge4 <- bind_cols(spread, spread2)

# Remove column duplicate
merge4 <- merge4[-6]

# Save
write.csv(merge4, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_2_smoking_CPG_data_processing_part_one_for_correlations_outputs/smoking_CPG_data_wide_format_by_wave_and_region.csv")


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











