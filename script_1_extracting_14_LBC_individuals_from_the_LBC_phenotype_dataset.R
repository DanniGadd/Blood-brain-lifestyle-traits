############################################################################################
############################################################################################
################# SCRIPT 1 - Danni Gadd - Blood vs brain DNAm project ######################
############################################################################################
############################################################################################

############################################################################################

setwd("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts")
library(haven)

############################################################################################

### LOAD DATA 

############################################################################################

## LBC 14 individuals subject IDs
IDs <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/Data/post_mortem_DNAm_extraction_Anna_14April2020/LBC_IDs.csv")

## Data request from the LBC team
LBC_phenotypes <- read_sav("/Volumes/marioni-lab/Danni/Blood_brain_project/Data/LBC_data_request_dataset_22_April_2020/LBC1936_Blood_DNAm_And_Brain_DNAm_RM_22APR2020.sav")

############################################################################################

### EXRACT DATA

############################################################################################

test <- as.data.frame(LBC_phenotypes)

merge <- merge(IDs, test, by.x = "LBC_ID", by.y = "lbc36no")

# Save this as an output 
write.csv(merge, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_1_extracting_14_LBC_individuals_from_the_LBC_phenotype_dataset_outputs/LBC_phenotypes_dataset_14_individuals_haven.csv")

############################################################################################

### SENSE CHECK PHENOTYPE READ FILE METHODS 

############################################################################################

# An additional sense check to make sure the read.spss (foreign package) and the new read_sav (haven) are producing the same dataset 
test <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_1_extracting_14_LBC_individuals_from_the_LBC_phenotype_dataset_outputs/LBC_phenotypes_dataset_14_individuals_haven.csv")
test2 <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_1_extracting_14_LBC_individuals_from_the_LBC_phenotype_dataset_outputs/LBC_phenotypes_dataset_14_individuals.csv")

identical(test,test2)

diff <- setdiff(test,test2)

write.csv(merge, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_1_extracting_14_LBC_individuals_from_the_LBC_phenotype_dataset_outputs/LBC_phenotypes_dataset_14_individuals.csv")


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
# [1] haven_2.2.0

# loaded via a namespace (and not attached):
#  [1] readr_1.3.1     compiler_3.6.3  R6_2.4.1        ellipsis_0.3.0 
#  [5] magrittr_1.5    hms_0.5.2       tools_3.6.3     pillar_1.4.4   
#  [9] tibble_3.0.1    crayon_1.3.4    Rcpp_1.0.3      vctrs_0.2.4    
# [13] forcats_0.5.0   lifecycle_0.2.0 pkgconfig_2.0.3 rlang_0.4.6 


############################################################################################

































