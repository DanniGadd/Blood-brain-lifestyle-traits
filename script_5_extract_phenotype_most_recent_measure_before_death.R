############################################################################################
############################################################################################
################# SCRIPT 5 - Danni Gadd - Blood vs brain DNAm project ######################
############################################################################################
############################################################################################

### SET WORKING DIRECTORY 
setwd("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts")

### LOAD ANY PACKAGES USED UP HERE 
library(tidyverse)
library(dplyr)

############################################################################################

### LOAD DATA

############################################################################################

# Now get the agedaysblood at each time point 
LBC_phenos <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_1_extracting_14_LBC_individuals_from_the_LBC_phenotype_dataset_outputs/LBC_phenotypes_dataset_14_individuals.csv")

###########################################################################################

### GET TRAIT PHENOTYPE DATA THAT IS MOST RECENT PRIOR TO DEATH FOR ALL TRAITS USED

###########################################################################################

# Get phenotype info from all available waves
BMI <- LBC_phenos[c(2,25,26,27)]
ALC <- LBC_phenos[c(2,15,18,23)]
HDL <- LBC_phenos[c(2,36,45,52)]
SMOK <- LBC_phenos[c(2,68,75)]

### GET MOST RECENT SMOKING INFO FIRST 

# Get most recent phenotype info prior to death  
SMOK$smokcurr_w4[is.na(SMOK$smokcurr_w4)] <- SMOK$smokcurr_w3[is.na(SMOK$smokcurr_w4)]

# Write this out so you can check which waves data was from
write.csv(SMOK, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_5_extract_phenotype_most_recent_measure_before_death_outputs/smoking_cat_wave_split.csv")


# Grab stuff I need 
SMOK <- SMOK[c(1,3)]

# Rename
names(SMOK)[2] <- "smokcat_recent"

# Write this to csv so i can use it for the smoking CPG analysis (as i was using wave 1 before)
write.csv(SMOK, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_5_extract_phenotype_most_recent_measure_before_death_outputs/smoking_cat_most_recent.csv")


### GET MOST RECENT ALCOHOL 

# Sort out comments and replace with NA
ALC$alcunitwk_w4[c(7,11)] <- NA
ALC$alcunitwk_w3[c(7,4,8)] <- NA
ALC$alcunitwk_w2[c(7)] <- NA

# Write this out so you can check which waves data was from
write.csv(ALC, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_5_extract_phenotype_most_recent_measure_before_death_outputs/alcohol_info_wave_split.csv")


ALC$alcunitwk_w4[is.na(ALC$alcunitwk_w4)] <- ALC$alcunitwk_w3[is.na(ALC$alcunitwk_w4)]

ALC$alcunitwk_w4[is.na(ALC$alcunitwk_w4)] <- ALC$alcunitwk_w2[is.na(ALC$alcunitwk_w4)]

# Individual at row 7 should be 0 as it was unrecordably low levels so I think this constitutes 0 as they clearly did drink from the comment just very little
# In the original file it says "alcohol levels too low to clasify"

ALC$alcunitwk_w4[c(7)] <- 0

# Grab alcohol stuff I need that is now most recent 
ALC <- ALC[c(1,4)]

# name properly 
names(ALC)[2] <- "allunitwk_recent"

# Write this to csv for reference
write.csv(ALC, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_5_extract_phenotype_most_recent_measure_before_death_outputs/alcohol_info_most_recent.csv")


### GET MOST RECENT HDL

# Write this out so you can check which waves data was from
write.csv(HDL, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_5_extract_phenotype_most_recent_measure_before_death_outputs/HDL_info_wave_split.csv")


# replace wave 4 with anythign from wave 3 
HDL$bld_hdlchol_w4[is.na(HDL$bld_hdlchol_w4)] <- HDL$bld_hdlchol_w3[is.na(HDL$bld_hdlchol_w4)]

# replace remaining NA values in wave 4 with wave 2 info
HDL$bld_hdlchol_w4[is.na(HDL$bld_hdlchol_w4)] <- HDL$bld_hdlchol_w2[is.na(HDL$bld_hdlchol_w4)]

# Grab what I need 
HDL <- HDL[c(1,4)]

# Rename propelry 
names(HDL)[2] <- "HDL_recent"

# Save file as per others for reference
write.csv(HDL, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_5_extract_phenotype_most_recent_measure_before_death_outputs/HDL_info_most_recent.csv")

### GET MOST RECENT BMI

# Replace the NA values in wave 4 with wave 3 values 
BMI$bmi_w4[is.na(BMI$bmi_w4)] <- BMI$bmi_w3[is.na(BMI$bmi_w4)]

# Grab what I need
BMI <- BMI[c(1,4)]

# Name properly
names(BMI)[2] <- "BMI_recent"

# Save that as reference file 
write.csv(BMI, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_5_extract_phenotype_most_recent_measure_before_death_outputs/BMI_info_wave_split.csv")

# Save that as reference file 
write.csv(BMI, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_5_extract_phenotype_most_recent_measure_before_death_outputs/BMI_info_most_recent.csv")


### PACK YEARS 

## number of packs of cigarettes smoked per day (20 per pack) multiplied by the number of years the person has smoked

# Grab packs relevant smoking info
packs <- LBC_phenos[c(2,55:76)]

# Narrow down so its easier to read 

packs <- packs[c(1,3,4,5,12,13,16,17,20,23)]

# Get most recent smoking status 

# Join up to take a look properly 
packs <- merge(packs, SMOK, by = "LBC_ID")

# Okay person 6 and 8 have never smoked - so we dont need to worry about them but will save them for now 
never <- packs %>% filter(smokcat_recent == "0")

# save so i have wave identifier 
write.csv(packs, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_5_extract_phenotype_most_recent_measure_before_death_outputs/packs_info_wave_split.csv")

# Now lets try to get the years these people have smoked for 
# To do this ill split them off into current or former smokers 
current <- packs %>% filter(smokcat_recent == "2")
former <- packs %>% filter(smokcat_recent == "1")

## FORMER - lets do this first 

# Populate W4 with any from wave 3 
former$smokagestop_w4[is.na(former$smokagestop_w4)] <- former$smokagestop_w3[is.na(former$smokagestop_w4)]

# Populate W4 with any from wave 2
former$smokagestop_w4[is.na(former$smokagestop_w4)] <- former$smokagestop_w2[is.na(former$smokagestop_w4)]

# Populate W4 with any from wave 1
former$smokagestop_w4[is.na(former$smokagestop_w4)] <- former$smokagestop_w1[is.na(former$smokagestop_w4)]

# As one person has an end age of 30 but no start age, ill impute the mean start age for this individual
mean_start <- mean(as.numeric(former$smokagestart_w1), na.rm = TRUE)

# mean age for start is 16.5 so ill impute this 
former$smokagestart_w1[is.na(former$smokagestart_w1)] <- mean_start

# Now calculate the time between starting the stopping smoking for each person
former$years_smoked <- former$smokagestop_w4 - former$smokagestart_w1 

# Now add number of cigarettes to W3 (no W4) column to get most recent data for this too 
former$smoknumcigs_w3[is.na(former$smoknumcigs_w3)] <- former$smoknumcigs_w2[is.na(former$smoknumcigs_w3)]

# And W1 to W3 now 
former$smoknumcigs_w3[is.na(former$smoknumcigs_w3)] <- former$smoknumcigs_w1[is.na(former$smoknumcigs_w3)]

# Work out packs per day for those that it is possible for 
former$packs_day <- former$smoknumcigs_w3 / 20

# Work out pack years for those that it is possibele for 
former$pack_years <- former$years_smoked * former$packs_day

# Select what I will need 
former_pack <- former %>% select("LBC_ID", "pack_years")



## CURRENT - now do this one 

# Get age to figure out the age difference from starting smoking to most recent
ages <- LBC_phenos[c(2,4:7)]

# Get most recent age measure in W4 column
ages$agedays_w4[is.na(ages$agedays_w4)] <- ages$agedays_w3[is.na(ages$agedays_w4)]
ages <- ages[c(1,5)]

# Join up with current 
current <- left_join(current, ages, by = "LBC_ID")

# Convert age days recent to age years recent 
current$agedays_w4 <- current$agedays_w4 / 365.25

# Work out the number of years smokers had been smoking for 
current$years_smoking <- current$agedays_w4 - current$smokagestart_w1

# Now get the most recent measure for the number of cigarettes smoked
current$smoknumcigs_w4[is.na(current$smoknumcigs_w4)] <- current$smoknumcigs_w3[is.na(current$smoknumcigs_w4)]

# Now work out the packs smoked per day (/20)
current$packs_smoked <- current$smoknumcigs_w4 / 20 

# Now calculate pack years for these current smokers 
current$pack_years <- current$years_smoking * current$packs_smoked

current_pack <- current %>% select("LBC_ID", "pack_years")

# Add the 2 individuals back in who never smoked (i.e. have a pack_years of 0)
never$pack_years <- 0

never_pack <- never %>% select("LBC_ID", "pack_years")

pack_years <- rbind(current_pack, never_pack)

# Now add whatever is possible in from the former group 
pack_years <- rbind(pack_years, former_pack)

# Save that as reference file 
write.csv(pack_years, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_5_extract_phenotype_most_recent_measure_before_death_outputs/pack_years_recent.csv")


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


















