############################################################################################
############################################################################################
################# SCRIPT 6 - Danni Gadd - Blood vs brain DNAm project ######################
############################################################################################
############################################################################################

setwd("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts")
library(tidyverse)
library(dplyr)

############################################################################################

### LOAD DATA

############################################################################################

# Load in extra data generated that will need to be joined up to main DNAm datasets for brain and blood 
smoking_cat <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_5_extract_phenotype_most_recent_measure_before_death_outputs/smoking_cat_most_recent.csv")
smoking_cat <- smoking_cat[-1]
pack_years <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_5_extract_phenotype_most_recent_measure_before_death_outputs/pack_years_recent.csv")
pack_years <- pack_years[-1]
recent_blood <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_3_extract_most_recent_blood_wave_and_age_at_which_blood_was_taken_outputs/blood_DNAm_recent_with_age_at_recent.csv")
recent_blood <- recent_blood[-1]
age_days_death <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_4_extract_age_at_death_for_inclusion_using_updated_death_info_outputs/age_days_at_death_data_LBC_updated_deaths.csv")
age_days_death <- age_days_death[-1]

# Load in proportion of Neurons variable for brain samples to add in too
prop_neuron <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/Data/LBC_propneuron_from_Anna_19_May_2020/LBC_propneuron.csv")

## Brain methylation dataset 
brain <- readRDS("/Volumes/marioni-lab/Danni/Blood_brain_project/Data/post_mortem_DNAm_extraction_Anna_14April2020/brain_methylation_data.rds")
brain <- as.data.frame(brain)

## Blood methylation dataset 
blood <- readRDS("/Volumes/marioni-lab/Danni/Blood_brain_project/Data/pre_mortem_DNAm_extraction_REM_14April2020/LBC1936_DNAm_blood_14longitudinal_Danni_14April2020.rds")
blood <- as.data.frame(blood)

############################################################################################

## CREATE A PHNEOTYPE DATASET THAT CAN BE USED FOR REGRESSIONS

############################################################################################

# Change classess in methylation datasets
brain <- brain %>% mutate_if(is.factor, as.character) -> brain
brain <- brain %>% mutate_if(is.integer, as.numeric) -> brain
blood <- blood %>% mutate_if(is.factor, as.character) -> blood
blood <- blood %>% mutate_if(is.integer, as.numeric) -> blood

# Select variables of interest 
select <- brain %>% 
	select(LBC_ID, region, PMI, sex, age, brain_pH, brain_weight_g, ApoE_genotype, cg05575921)


# Do the same for blood
select2 <- blood %>% 
	select(ID, WAVE, sex, age, neut, lymph, mono, eosin, baso, cg05575921)

# Merge smoking info with dataset
brain_s <- merge(select, smoking_cat, by = "LBC_ID")
brain_s <- merge(brain_s, pack_years, by = "LBC_ID")

# Adding most recent blood DNAm and age
brain_sa <- merge(brain_s, recent_blood, by = "LBC_ID")

# Merge blood data based on LBC_ID and Wave identifier in order to get the cell counts in the main dataset 
brain_sa <- brain_sa %>% mutate_if(is.integer, as.numeric) -> brain_sa
join <- left_join(brain_sa, select2, by = c("LBC_ID" = "ID", "Wave_id" = "WAVE"))
join <- join[c(-15,-16,-22)]
names(join)[4] <- "sex"
names(join)[5] <- "age"
names(join)[9] <- "DNAm_brain_death"

# Add age at death in 
join2 <- merge(join, age_days_death, by = "LBC_ID")

# Time to death calculated by age at death minus age at blood CpG measurement
join2 <- join2 %>% mutate_if(is.integer, as.numeric) -> join2
join2$time_to_death <- join2$agedays_death_new - join2$age_blood_recent

# Change prop_neuron variable to chracter vs factor for joining 
prop_neuron <- prop_neuron %>% mutate_if(is.factor, as.character) -> prop_neuron

# Add prop_neuron variable in for regression analysis
join3 <- left_join(join2, prop_neuron, by = c("LBC_ID" = "LBC_ID", "region" = "region"))

# Save this version output 
write.csv(join3 , file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_6_smoking_CPG_data_processing_part_two_phenotype_joining_for_regressions_outputs/smoking_CPG_data_for_statistical_testing_most_recent_blood_data_only.csv")


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


















