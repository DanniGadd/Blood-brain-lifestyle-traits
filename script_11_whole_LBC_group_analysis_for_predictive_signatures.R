############################################################################################
############################################################################################
################# SCRIPT 17 - Danni Gadd - Blood vs brain DNAm project ######################
############################################################################################
############################################################################################

setwd("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts")
library(tidyverse)
library(dplyr)
library(haven)
library(readxl)
library(plyr)
library(tibble)
library(lubridate)

############################################################################################

### LOAD DATA

############################################################################################

# Extracted dataset for all LBC individuals available, with CpGs of interest for the traits included 
blood <- readRDS("/Volumes/marioni-lab/Danni/Blood_brain_project/Data/Extraction_of_LBC_total_group_19_May_2020/LBC1936_DNAm_blood_all_longitudinal_Danni_19May2020.rds")

## LBC 14 individuals subject IDs
IDs <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/Data/post_mortem_DNAm_extraction_Anna_14April2020/LBC_IDs.csv")

## Data request from the LBC team
LBC_phenotypes <- read_sav("/Volumes/marioni-lab/Danni/Blood_brain_project/Data/LBC_data_request_dataset_22_April_2020/LBC1936_Blood_DNAm_And_Brain_DNAm_RM_22APR2020.sav")

# list of 14 LBC1936 brain donor participants
brain <- c("LBC360291", "LBC360232", "LBC360999", "LBC360202", "LBC361079", "LBC360849", "LBC360021", "LBC360873", "LBC360666", "LBC360705", "LBC361021", "LBC360620", "LBC360993", "LBC360838")


###########################################################################################

### PREP FILES FOR ANALYSIS IN WAVE 4

###########################################################################################

### BLOOD DNAm 

# Check dimensions of the original group first 
dim(blood) # [1] 3525 2413 

# Get wave 4 only 
LBC_w4 <- subset(blood, WAVE == 4)

LBC_W4 <- blood[blood$WAVE==4 & blood$cohort=="LBC36",]

# Check dimentions
dim(LBC_w4) # [1]  589 2413 (589 people in wave 4)

# Remove 14 people which we have DNAm in brain
blood_meth <- LBC_w4[!LBC_w4$ID %in% brain, ]

# Check if they were removed
dim(blood_meth) # [1]  583 2413

### PHENOTYPES 

# Save down phenotype data 
phenos <- as.data.frame(LBC_phenotypes)
dim(phenos) # [1] 1091   83

# Subset to W4 group 
matched_phenos <- phenos$lbc36no %in% blood_meth$ID
phenos <- phenos[phenos$lbc36no %in% blood_meth$ID,]
dim(phenos)
# [1] 500  83

# Find overlap
overlap <- blood_meth$ID %in% phenos$lbc36no
blood_subset <- blood_meth[overlap,]

dim(blood_subset) 
# [1]  500 2413
# 500 people in my final blood subset! for W4

overlap <- which(blood_subset$ID %in% "LBC360420")
blood_subset <- blood_subset[-overlap,]


# Save W4 n 500 methylation and phenotype files for reference 

write.csv(phenos, file="/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/wave_4_phenotype_n_500_dataset.csv")

write.csv(blood_subset, file="/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/wave_4_DNAm_subset_n_500_dataset.csv")



### FORMATTING BLOOD FOR PREDICTOR CALCULATIONS

# Include only methylation data, with the LBC IDs to identify people, with wave included
blood <- blood_subset[c(5,20:2413)]
names(blood)[1] <- "LBC_ID"
blood$LBC_ID <- as.character(blood$LBC_ID)


# Transpose so that LBC_IDs are the column names 

# keep the first column 
names <-  blood[,1]

# Transpose everything other than the first column
t_blood <- as.data.frame(as.matrix(t(blood[,-1])))

# Assign first column as the column names of the transposed dataframe
colnames(t_blood) <- names

# Save this file 
write.csv(t_blood, "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/List_files/blood_methylation_most_recent_only_LBC_WHOLE_GROUP_IDs_as_colnames.csv")

###########################################################################################

### CALCULATE PREDICTOR SCORES FOR BLOOD - READ IN FILES NEEDED

###########################################################################################

# Save this file 
blood_csv <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/List_files/blood_methylation_most_recent_only_LBC_WHOLE_GROUP_IDs_as_colnames.csv")


# Save sample files as a list to read into functions below
files = list.files("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/List_files", pattern="*.csv", full.names = TRUE)

# Make sure this order is matched in the folder with the listed files
file_order <- c("blood")


###########################################################################################

## BLOOD - SMOKING PREDICTOR SCORES

###########################################################################################


# Use lapply to do the same steps for predictor score generation for each file 
dfList <- lapply(files, function(i) {
     df <- read.csv(i) # read in sample DNAm data file
     weight <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_7_extract_predictor_weights_for_traits_outputs/SMOK_weights.csv") # read in smoking weights
     df <- df %>% filter(X %in% weight$CpG) # filter the sample dataset by the CpG sites in the predictor weight file
     not_names <- setdiff(weight$CpG, df$X) # find out if there are any CpGs in the weights file that arent in the sample dataset
     weight_re <- weight[weight$CpG %in% df$X,] # remove the CpGs which dont match up from the predictor weights 
     re_sm_match <- df[match(weight_re$CpG,df$X),] # make sure the order of CpGs now matches between the files
     re_weights <- re_sm_match[,c(2:500)]*weight_re[,3] # multiply the methylation values with the predicto weights values 
     re_sum <- colSums(re_weights) # sum the columns to get a predictor score for each individual
     return(re_sum) # output scores
})

smoking_predictors <- do.call(rbind, dfList) # bind scores together so that each sample dataset has one row in the resulting table 
smoking_predictors <- as.data.frame(smoking_predictors)
smoking_predictors$files <- file_order
smoking_predictors <- smoking_predictors[c(500,1:499)]


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
     re_weights <- re_sm_match[,c(2:500)]*weight_re[,3] # multiply the methylation values with the predicto weights values 
     re_sum <- colSums(re_weights) # sum the columns to get a predictor score for each individual
     return(re_sum) # output scores
})

alcohol_predictors <- do.call(rbind, dfList) # bind scores together so that each sample dataset has one row in the resulting table 
alcohol_predictors <- as.data.frame(alcohol_predictors)
alcohol_predictors$files <- file_order
alcohol_predictors <- alcohol_predictors[c(500,1:499)]

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
     re_weights <- re_sm_match[,c(2:500)]*weight_re[,3] # multiply the methylation values with the predicto weights values 
     re_sum <- colSums(re_weights) # sum the columns to get a predictor score for each individual
     return(re_sum) # output scores
})

BMI_predictors <- do.call(rbind, dfList) # bind scores together so that each sample dataset has one row in the resulting table 
BMI_predictors <- as.data.frame(BMI_predictors)
BMI_predictors$files <- file_order
BMI_predictors <- BMI_predictors[c(500,1:499)]

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
     re_weights <- re_sm_match[,c(2:500)]*weight_re[,3] # multiply the methylation values with the predicto weights values 
     re_sum <- colSums(re_weights) # sum the columns to get a predictor score for each individual
     return(re_sum) # output scores
})

HDL_predictors <- do.call(rbind, dfList) # bind scores together so that each sample dataset has one row in the resulting table 
HDL_predictors <- as.data.frame(HDL_predictors)
HDL_predictors$files <- file_order
HDL_predictors <- HDL_predictors[c(500,1:499)]



###########################################################################################

### SAVE PREDICTOR SCORE FILES 

###########################################################################################


write.csv(smoking_predictors, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/Predictor_scores/smoking_predictors.csv")
write.csv(alcohol_predictors, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/Predictor_scores/alcohol_predictors.csv")
write.csv(BMI_predictors, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/Predictor_scores/BMI_predictors.csv")
write.csv(HDL_predictors, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/Predictor_scores/HDL_predictors.csv")



###########################################################################################

### PHENOTYPE ANALYSIS - GETTING THE MOST RECENT PHENOTYPES FOR INCLUSION

###########################################################################################

### PHENOTYPES 

# Save down phenotype data 
phenos <- as.data.frame(LBC_phenotypes)
dim(phenos) # [1] 1091   83

# Subset to W4 group 
matched_phenos <- phenos$lbc36no %in% blood_meth$ID
phenos <- phenos[phenos$lbc36no %in% blood_meth$ID,]
dim(phenos)
# [1] 500  83


overlap <- which(phenos$lbc36no %in% "LBC360420")
phenos <- phenos[-overlap,] # 499 update 


### Get the most recent phenotype info for the individuals remining 

# Extract phenotypic-relevant variables for each trait 
BMI <- phenos[c(1,23:26)]
ALC <- phenos[c(1,11,14,17,22)]
HDL <- phenos[c(1,35,42,44,51)]
SMOK <- phenos[c(1,55,63,67,74)]

# Sort out weird haven labelling
SMOK$smokever_w1 <- as.numeric(SMOK$smokever_w1)
SMOK$smokcurr_w2 <- as.numeric(SMOK$smokcurr_w2)
SMOK$smokcurr_w3 <- as.numeric(SMOK$smokcurr_w3)
SMOK$smokcurr_w4 <- as.numeric(SMOK$smokcurr_w4)

# Sort out weird haven labelling
ALC$alcunitwk_w1 <- as.numeric(ALC$alcunitwk_w1)
ALC$alcunitwk_w2 <- as.numeric(ALC$alcunitwk_w2)
ALC$alcunitwk_w3 <- as.numeric(ALC$alcunitwk_w3)
ALC$alcunitwk_w4 <- as.numeric(ALC$alcunitwk_w4)

####################################################################################################################################

# w4 smoking only 
SMOK_w4 <- SMOK[c(1,5)]

# name properly 
names(SMOK_w4)[2] <- "smoking_w4"

# Write this to csv so i can use it for the smoking CPG analysis (as i was using wave 1 before)
write.csv(SMOK_w4, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/smoking_cat_w4.csv")


####################################################################################################################################

### MOST RECENT SMOKING 

# Anything not present in W4 gets W3 assigned 
SMOK$smokcurr_w4[is.na(SMOK$smokcurr_w4)] <- SMOK$smokcurr_w3[is.na(SMOK$smokcurr_w4)]

# Now anything that is not present in W4 still gets W2 assigned value 
SMOK$smokcurr_w4[is.na(SMOK$smokcurr_w4)] <- SMOK$smokcurr_w2[is.na(SMOK$smokcurr_w4)]

# Now anything that is not present in W4 still gets W1 assigned value 
SMOK$smokcurr_w4[is.na(SMOK$smokcurr_w4)] <- SMOK$smokever_w1[is.na(SMOK$smokcurr_w4)]

# Grab stuff I need 
SMOK <- SMOK[c(1,5)]

# Write this to csv so i can use it for the smoking CPG analysis (as i was using wave 1 before)
write.csv(SMOK, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/smoking_cat_most_recent.csv")

####################################################################################################################################

### ALCOHOL 

# Replace -99.000 with NA
ALC[ALC == -99.000] <- NA

####################################################################################################################################

# w4 ALC only 
ALC_w4 <- ALC[c(1,5)]

# name properly 
names(ALC_w4)[2] <- "allunitwk_w4"

# Write this to csv for reference
write.csv(ALC_w4, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/alcohol_info_w4.csv")


####################################################################################################################################

# Anything not present in W4 gets W3 assigned 
ALC$alcunitwk_w4[is.na(ALC$alcunitwk_w4)] <- ALC$alcunitwk_w3[is.na(ALC$alcunitwk_w4)]

# Anything not present in W4 gets W2 now assigned 
ALC$alcunitwk_w4[is.na(ALC$alcunitwk_w4)] <- ALC$alcunitwk_w2[is.na(ALC$alcunitwk_w4)]

# Anything not present in W4 gets W1 now assigned 
ALC$alcunitwk_w4[is.na(ALC$alcunitwk_w4)] <- ALC$alcunitwk_w1[is.na(ALC$alcunitwk_w4)]


# Grab alcohol stuff I need that is now most recent 
ALC <- ALC[c(1,5)]

# change character info to numeric for ALC file 
ALC$alcunitwk_w4 <- as.numeric(ALC$alcunitwk_w4)

# name properly 
names(ALC)[2] <- "allunitwk_recent"

# Write this to csv for reference
write.csv(ALC, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/alcohol_info_most_recent.csv")

####################################################################################################################################

### HDL

# Replace -999.0 with NA as its clealry invalid
HDL[HDL == -999.0] <- NA

####################################################################################################################################

# w4 HDL
HDL_w4 <- HDL[c(1,5)]

# Rename propelry 
names(HDL_w4)[2] <- "HDL_w4"

# Save file as per others for reference
write.csv(HDL_w4, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/HDL_info_w4.csv")


####################################################################################################################################


# replace wave 4 with anythign form wave 3 
HDL$bld_hdlchol_w4[is.na(HDL$bld_hdlchol_w4)] <- HDL$bld_hdlchol_w3[is.na(HDL$bld_hdlchol_w4)]

# replace remaining NA values in wave 4 with wave 2 info
HDL$bld_hdlchol_w4[is.na(HDL$bld_hdlchol_w4)] <- HDL$bld_hdlchol_w2[is.na(HDL$bld_hdlchol_w4)]

# replace remaining NA values in wave 4 with wave 1 info
HDL$bld_hdlchol_w4[is.na(HDL$bld_hdlchol_w4)] <- HDL$bld_hdlchol_w1[is.na(HDL$bld_hdlchol_w4)]

# Grab what I need 
HDL <- HDL[c(1,5)]

# Rename propelry 
names(HDL)[2] <- "HDL_recent"

# Save file as per others for reference
write.csv(HDL, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/HDL_info_most_recent.csv")


####################################################################################################################################

### BMI

# w4 BMI
BMI_w4 <- BMI[c(1,5)]

# Name properly
names(BMI_w4)[2] <- "BMI_w4"

# Save that as reference file 
write.csv(BMI_w4, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/BMI_info_w4.csv")


####################################################################################################################################


# Replace the NA values in wave 4 with wave 3 values 
BMI$bmi_w4[is.na(BMI$bmi_w4)] <- BMI$bmi_w3[is.na(BMI$bmi_w4)]

# Replace the NA values in wave 4 with wave 2 values 
BMI$bmi_w4[is.na(BMI$bmi_w4)] <- BMI$bmi_w2[is.na(BMI$bmi_w4)]

# Replace the NA values in wave 4 with wave 1 values 
BMI$bmi_w4[is.na(BMI$bmi_w4)] <- BMI$bmi_w1[is.na(BMI$bmi_w4)]

# Grab what I need
BMI <- BMI[c(1,5)]

# Name properly
names(BMI)[2] <- "BMI_recent"

# Save that as reference file 
write.csv(BMI, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/BMI_info_most_recent.csv")

####################################################################################################################################

### PACK YEARS 

# w4 pack years 
# Grab packs relevant smoking info
packs <- phenos[c(1,54:75)]

# Narrow down so its easier to read 
packs <- packs[c(1,4,5,6,12,13,16,17,20,23)]

# Get most recent smoking status added in  
# Join up to take a look properly 
packs_w4 <- merge(packs, SMOK_w4, by = "lbc36no")

# Now calculate the time between starting the stopping smoking for each person
packs_w4$years_smoked <- packs_w4$smokagestop_w4 - packs_w4$smokagestart_w1 

# Work out packs per day for those that it is possible for 
packs_w4$packs_day <- packs_w4$smoknumcigs_w4 / 20

# Work out pack years for those that it is possibele for 
packs_w4$pack_years <- packs_w4$years_smoked * packs_w4$packs_day

# Grab relevant stuff 
packs_w4 <- packs_w4[c(1,11,14)]

# name
names(packs_w4)[3] <- "pack_years_w4"

# If someone has never smoked (0) then we should assign 0 pack years 

# Categorise a never smoked information variable only 
packs_w4 <- packs_w4%>% mutate(never = case_when(
     smoking_w4 == "0" ~ "0",
     smoking_w4 == "1" ~ "NA",
     smoking_w4 == "2" ~ "NA"))

# Convert never 
packs_w4$pack_years_w4 <- as.numeric(packs_w4$pack_years_w4)

packs_w4$never <- as.numeric(packs_w4$never)

# Merge never smoked values into pack years as intended
packs_w4$pack_years_w4[is.na(packs_w4$pack_years_w4)] <- packs_w4$never[is.na(packs_w4$pack_years_w4)]


# Save that as reference file 
write.csv(packs_w4, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/pack_years_info_w4.csv")

####################################################################################################################################

# Most recent pack years 

## number of packs of cigarettes smoked per day (20 per pack) multiplied by the number of years the person has smoked

# Grab packs relevant smoking info
packs <- phenos[c(1,54:75)]

# Narrow down so its easier to read 
packs <- packs[c(1,4,5,6,12,13,16,17,20,23)]

# Get most recent smoking status added in  
# Join up to take a look properly 
packs <- merge(packs, SMOK, by = "lbc36no")

# Populate W4 with any from wave 3 
packs$smokagestop_w4[is.na(packs$smokagestop_w4)] <- packs$smokagestop_w3[is.na(packs$smokagestop_w4)]

# Populate W4 with any from wave 2
packs$smokagestop_w4[is.na(packs$smokagestop_w4)] <- packs$smokagestop_w2[is.na(packs$smokagestop_w4)]

# Populate W4 with any from wave 1
packs$smokagestop_w4[is.na(packs$smokagestop_w4)] <- packs$smokagestop_w1[is.na(packs$smokagestop_w4)]

# Now calculate the time between starting the stopping smoking for each person
packs$years_smoked <- packs$smokagestop_w4 - packs$smokagestart_w1 

# Now add number of cigarettes from W3 to W4
packs$smoknumcigs_w4[is.na(packs$smoknumcigs_w4)] <- packs$smoknumcigs_w3[is.na(packs$smoknumcigs_w4)]

# Now add number of cigarettes from W2 to W4 
packs$smoknumcigs_w4[is.na(packs$smoknumcigs_w4)] <- packs$smoknumcigs_w2[is.na(packs$smoknumcigs_w4)]

# And W1 to W4 
packs$smoknumcigs_w4[is.na(packs$smoknumcigs_w4)] <- packs$smoknumcigs_w1[is.na(packs$smoknumcigs_w4)]

# Work out packs per day for those that it is possible for 
packs$packs_day <- packs$smoknumcigs_w4 / 20

# Work out pack years for those that it is possibele for 
packs$pack_years <- packs$years_smoked * packs$packs_day

# Grab what I need 
packs <- packs[c(1,11,14)]

# Change label class
packs$pack_years <- as.numeric(packs$pack_years)

# If someone has never smoked (0) then we should assign 0 pack years 

# Categorise a never smoked information variable only 
packs <- packs%>% mutate(never = case_when(
     smokcurr_w4 == "0" ~ "0",
     smokcurr_w4 == "1" ~ "NA",
     smokcurr_w4 == "2" ~ "NA"))

# Merge never smoked values into pack years as intended
packs$pack_years[is.na(packs$pack_years)] <- packs$never[is.na(packs$pack_years)]

# Grab relevant stuff 
packs <- packs[c(1:3)]

# Save that as reference file 
write.csv(packs, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/pack_years_info_most_recent.csv")

####################################################################################################################################




###########################################################################################

### JOIN PREDICTORS UP TO MAIN DATASET NOW 

############################################################################################

# Predictor scores 
smoking <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/Predictor_scores/smoking_predictors.csv")
alcohol <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/Predictor_scores/alcohol_predictors.csv")
BMI <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/Predictor_scores/BMI_predictors.csv")
HDL <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/Predictor_scores/HDL_predictors.csv")

# Save these predictor scores as a list too, so they can be read into the functions below
files = list.files("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/Predictor_scores", pattern="*.csv", full.names = TRUE)

# Get order of files 
# [1] "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/Predictor_scores/alcohol_predictors.csv"
# [2] "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/Predictor_scores/BMI_predictors.csv"    
# [3] "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/Predictor_scores/HDL_predictors.csv"    
# [4] "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/Predictor_scores/smoking_predictors.csv"


## Transposing the predictor data 

# Use lapply to do the same steps for predictor score transposing for each file 
dfList <- lapply(files, function(i) {
     df <- read.csv(i) # read in predictor score file 
     df <- df[-1] # get rid of first column as its just junk 
     names <- df[,1] # get the names of the samples before transposing so that these match up 
     t_df <- as.data.frame(as.matrix(t(df[,-1]))) # transpose everything other than the first column
     colnames(t_df) <- names # assign first column as the names of the transposed dataframe 
     t_df$LBC_ID <- rownames(t_df) # Create LBC_ID column for merging
     return(t_df) # output the transformed datasets for the predictor scores
})

# Extract the ones I want as separate dataframes
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
names(alcohol)[2] <- "LBC_ID"
names(BMI)[2] <- "LBC_ID"
names(smoking)[2] <- "LBC_ID"
names(HDL)[2] <- "LBC_ID"

# Save predictor scores for merging 
write.csv(alcohol, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/alcohol_predictor_scores_for_merging.csv")
write.csv(smoking, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/smoking_predictor_scores_for_merging.csv")
write.csv(BMI, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/BMI_predictor_scores_for_merging.csv")
write.csv(HDL, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/HDL_predictor_scores_for_merging.csv")



###########################################################################################

### Get the main dataset sorted out with variables of interest here 

# Read in the main blood subset of 500 individuals 
blood_main <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/wave_4_DNAm_subset_n_500_dataset.csv")

# Convert dodgy classess
blood_main$ID <- as.character(blood_main$ID)
blood_main$sex <- as.character(blood_main$sex)

# Grab just the smoking CpGs only 
blood_main_selected <- blood_main %>% select("ID", "WAVE", "age", "sex", "cg05575921", "neut", "lymph", "mono", "baso", "eosin")
names(blood_main_selected)[1] <- "LBC_ID"



### Now merge in all predictor scores and trait info (which was generated above and saved)
# I need to make sure we keep this to the 466 with trait info still when merging for now 

# Get the predictor scores loaded back for ease 
alcohol <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/alcohol_predictor_scores_for_merging.csv")
smoking <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/smoking_predictor_scores_for_merging.csv")
BMI <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/BMI_predictor_scores_for_merging.csv")
HDL <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/HDL_predictor_scores_for_merging.csv")

# get rid of junk column
alcohol <- alcohol[-1]
smoking <- smoking[-1]
BMI <- BMI[-1]
HDL <- HDL[-1]

# make factor LBC_ID into character 
alcohol$LBC_ID <- as.character(alcohol$LBC_ID)
smoking$LBC_ID <- as.character(smoking$LBC_ID)
BMI$LBC_ID <- as.character(BMI$LBC_ID)
HDL$LBC_ID <- as.character(HDL$LBC_ID)

# Now load in the w4 trait info only group for sensitvity analysis 
alcohol_trait_w4 <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/alcohol_info_w4.csv")
smoking_trait_w4 <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/smoking_cat_w4.csv")
BMI_trait_w4 <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/BMI_info_w4.csv")
HDL_trait_w4 <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/HDL_info_w4.csv")
packs_w4 <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/pack_years_info_w4.csv")

# Get rid of junk first column (annoying)
alcohol_trait_w4 <- alcohol_trait_w4[-1]
smoking_trait_w4 <- smoking_trait_w4[-1]
BMI_trait_w4 <- BMI_trait_w4[-1]
HDL_trait_w4 <- HDL_trait_w4[-1]
packs_w4 <- packs_w4[c(2,4)]

# Rename LBC column
names(alcohol_trait_w4)[1] <- "LBC_ID"
names(smoking_trait_w4)[1] <- "LBC_ID"
names(BMI_trait_w4)[1] <- "LBC_ID"
names(HDL_trait_w4)[1] <- "LBC_ID"
names(packs_w4)[1] <- "LBC_ID"

# Join up everything to one dataset for predictor scores 
blood <- left_join(blood_main_selected, smoking, by = "LBC_ID")
blood <- left_join(blood, alcohol, by = "LBC_ID")
blood <- left_join(blood, BMI, by = "LBC_ID")
blood <- left_join(blood, HDL, by = "LBC_ID")

# Convert LBC_ID in the phenotype files to a character 
alcohol_trait_w4$LBC_ID <- as.character(alcohol_trait_w4$LBC_ID)
smoking_trait_w4$LBC_ID <- as.character(smoking_trait_w4$LBC_ID)
BMI_trait_w4$LBC_ID <- as.character(BMI_trait_w4$LBC_ID)
HDL_trait_w4$LBC_ID <- as.character(HDL_trait_w4$LBC_ID)
packs_w4$LBC_ID <- as.character(packs_w4$LBC_ID)

# Now add in phenotypes (most recent, then also w4 only - so i can do sensitivity analysis on them both)
blood <- left_join(blood, packs_w4, by = "LBC_ID")
blood <- left_join(blood, alcohol_trait_w4, by = "LBC_ID")
blood <- left_join(blood, smoking_trait_w4, by = "LBC_ID")
blood <- left_join(blood, BMI_trait_w4, by = "LBC_ID")
blood <- left_join(blood, HDL_trait_w4, by = "LBC_ID")

# Save this dataset!!! 
write.csv(blood, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/whole_group_LBC_blood_analysis_data_ready_for_plots_499.csv")


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
#  [1] lubridate_1.7.4 plyr_1.8.4      readxl_1.3.1    haven_2.2.0    
#  [5] forcats_0.5.0   stringr_1.4.0   dplyr_0.8.5     purrr_0.3.3    
#  [9] readr_1.3.1     tidyr_1.0.3     tibble_3.0.1    ggplot2_3.3.0  
# [13] tidyverse_1.3.0

# loaded via a namespace (and not attached):
#  [1] Rcpp_1.0.3       cellranger_1.1.0 pillar_1.4.4     compiler_3.6.3  
#  [5] dbplyr_1.4.2     tools_3.6.3      jsonlite_1.6     lifecycle_0.2.0 
#  [9] nlme_3.1-144     gtable_0.3.0     lattice_0.20-41  pkgconfig_2.0.3 
# [13] rlang_0.4.6      reprex_0.3.0     cli_1.1.0        rstudioapi_0.10 
# [17] DBI_1.0.0        withr_2.1.2      xml2_1.2.2       httr_1.4.1      
# [21] fs_1.4.1         generics_0.0.2   vctrs_0.2.4      hms_0.5.2       
# [25] grid_3.6.3       tidyselect_1.0.0 glue_1.4.0       R6_2.4.1        
# [29] modelr_0.1.5     magrittr_1.5     backports_1.1.5  scales_1.1.0    
# [33] ellipsis_0.3.0   rvest_0.3.5      assertthat_0.2.1 colorspace_1.4-1
# [37] stringi_1.4.3    munsell_0.5.0    broom_0.5.5      crayon_1.3.4    



############################################################################################

###########################################################################################
###########################################################################################
###########################################################################################

### END OF SCRIPT - MANUAL CHECK BELOW FOR LAPPLY FUNCTIONS ONLY AS REFERENCE

###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################

## BLOOD SMOKING PREDICTOR LAPPLY FUNCTION MANUAL CHECK 

###########################################################################################

# Read in smoking weights
SMOK <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_7_extract_predictor_weights_for_traits_outputs/SMOK_weights.csv")

# filter to smoking CpG sites only
blood_sm <- blood %>%
     filter(X %in% SMOK$CpG)

# Check previous
dim(blood)
# [1] 2394  895

# Check filtered 
dim(blood_sm)
# [1] 230 895

# find out if there are any CpGs in the weights file that arent in the sample dataset
not_names <- setdiff(SMOK$CpG, blood_sm$X) # 3 sites 

# remove the CpGs which dont match up from the predictor weights file
SMOK_bl <- SMOK[SMOK$CpG %in% blood_sm$X,]

# make sure the order of CpGs now matches between the files
blood_sm_match <- blood_sm[match(SMOK_bl$CpG,blood_sm$X),]

# Are there NA values
table(is.na(blood_sm_match))
#  FALSE 
# 205850 

# test whether matches up 
blood_sm_match <- blood_sm_match %>% mutate_if(is.integer, as.numeric) -> blood_sm_match
blood_sm_match <- blood_sm_match %>% mutate_if(is.factor, as.character) -> blood_sm_match
SMOK_bl <- SMOK_bl %>% mutate_if(is.integer, as.numeric) -> SMOK_bl
SMOK_bl <- SMOK_bl %>% mutate_if(is.factor, as.character) -> SMOK_bl

identical(blood_sm_match$X, SMOK_bl$CpG)

# Find the number of columns im multiplying by this time (i.e. people)
colnames(blood_sm_match)
ncol(blood_sm_match)
# [1] 895

# Multiply the methylation values with the predictor weights values 
blood_weights <- blood_sm_match[,c(2:895)]*SMOK_bl[,3]

blood_sum <- colSums(blood_weights)

blood_n <- as.data.frame(blood_sum)



###########################################################################################
























