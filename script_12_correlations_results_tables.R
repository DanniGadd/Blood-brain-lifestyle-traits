############################################################################################
############################################################################################
################# SCRIPT 12 - Danni Gadd - Blood vs brain DNAm project ######################
############################################################################################
############################################################################################

setwd("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts") 
library(tidyverse)
library(corrplot)
library(dplyr)
library(Hmisc)
library(ggplot2)
library(pwr)
library(ggpubr)
library(RcmdrMisc)

############################################################################################

### LOAD DATA

############################################################################################

# Load in main dataset as generated in script 10 with everything in it for anlayses
data <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_10_epigenetic_signature_data_processing_part_two_phenotype_joining_for_regressions_outputs/data_withmost_recent_clinical_phenotypes_added.csv")

# Load in the data I went back and generated in script 17 for the whole LBC group 
whole <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/whole_group_LBC_blood_analysis_data_ready_for_plots_499.csv")



############################################################################################

### CORRELATIONS FOR SMOKING CPG

############################################################################################

### For the 14 dataset 

# calualte mean brain measure for each person 
means <- data %>% select(c("LBC_ID", "cg05575921_brain")) 
mean <- tapply(means$cg05575921_brain, means$LBC_ID, mean)
mean <- as.data.frame(mean)
names(mean)[1] <- "mean_cg05575921_brain"

# get IDs as list and rejoin
IDs <- data %>% select("LBC_ID") %>% unique() 
mean <- cbind(mean, IDs)

# Join back to original dataset
data <- merge(data, mean, by = "LBC_ID")

# Get the data I need
corrdata <- data %>% select("LBC_ID", "region", "mean_cg05575921_brain", "cg05575921_brain", "cg05575921_blood", "smokcat_recent", "pack_years") %>%
unique()

# Save blood measures 
corrblood <- corrdata[c(1,5)] %>% unique()

# Get brain measures
corrbrain <- corrdata[c(1,2,4)] %>% unique()

# Spread the data for brain so I can correlate by brain region
spread <- corrbrain %>%
	spread(key = region, value = cg05575921_brain)

# Add mean brain measure to the spread brains 
spread <- merge(spread, mean, by = "LBC_ID")

# Join the blood most recent values back in now that we have variables for each brain region
corrdata <- merge(spread, corrblood, by = "LBC_ID")

# Extract smoking info
corrsmok <- data %>% select("LBC_ID", "smokcat_recent", "pack_years") %>% unique()

# Add to main correlation dataset 
corrdata <- merge(corrdata, corrsmok, by = "LBC_ID")

# make smoking category continuous 
corrdata$smokcat_recent <- as.numeric(corrdata$smokcat_recent)

# Get rid of LBC identifier to run 
corrdata <- corrdata[-1]

# Rename it so it looks nice in the plot 
names(corrdata)[1] <- "BA17 cg05575921"
names(corrdata)[2] <- "BA20/21 cg05575921"
names(corrdata)[3] <- "BA24 cg05575921"
names(corrdata)[4] <- "BA46 cg05575921"
names(corrdata)[5] <- "BA35 cg05575921"
names(corrdata)[6] <- "Mean brain cg05575921"
names(corrdata)[7] <- "Blood cg05575921"
names(corrdata)[8] <- "Smoking status"
names(corrdata)[9] <- "Pack years smoked"

# Remove mean brain - no longer needed 
corrdata <- corrdata[-6]

###########################################################################################

# Use method to get p values as well as same r values as above methods
test <- rcorr.adjust(as.matrix(corrdata), type = c("spearman"), use=c("pairwise.complete.obs"))

# Extract R table
str(test)
tableR <- test$R

str(tableR)
tabler <- tableR$r
tablep <- tableR$P

# Load function to tabulate r and p values
flat_cor_mat <- function(cor_r, cor_p){
  #This function provides a simple formatting of a correlation matrix
  #into a table with 4 columns containing :
    # Column 1 : row names (variable 1 for the correlation test)
    # Column 2 : column names (variable 2 for the correlation test)
    # Column 3 : the correlation coefficients
    # Column 4 : the p-values of the correlations
  library(tidyr)
  library(tibble)
  cor_r <- rownames_to_column(as.data.frame(cor_r), var = "row")
  cor_r <- gather(cor_r, column, cor, -1)
  cor_p <- rownames_to_column(as.data.frame(cor_p), var = "row")
  cor_p <- gather(cor_p, column, p, -1)
  cor_p_matrix <- left_join(cor_r, cor_p, by = c("row", "column"))
  cor_p_matrix
}

# Tabulate 
cor <- flat_cor_mat(tableR$r, tableR$P)

# Data frame it 
cor <- as.data.frame(cor)

# Rename
names(cor) <- c("Tissue_A", "Tissue_B", "r", "p")

# Make values rounded to 2 digits 
cor$r <- round(cor$r, digits = 2)
cor$p <- round(cor$p, digits = 3)

# Save results table 
write.csv(cor, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/Oct_2020_brain_comms/Correlation_results/smoking_cpg_spearmans.csv", row.names = F)

###########################################################################################

### For the whole group

# the 'w4' phenotype only n=499 group

# Get the data I need
corrdata <- whole %>% select("LBC_ID", "cg05575921", "smoking_w4", "pack_years_w4") %>%
unique()

# Make sure class correct 
corrdata$smoking_w4 <- as.numeric(corrdata$smoking_w4)

# Rename for results 
names(corrdata) <- c("LBC_ID", "Blood cg05575921", "Smoking status", "Pack years smoked")

# remove things not correlating
corrdata <- corrdata[-1]

# Use method to get p values as well as same r values as above methods
test <- rcorr.adjust(as.matrix(corrdata), type = c("spearman"), use=c("pairwise.complete.obs"))

# Extract R table
str(test)
tableR <- test$R

str(tableR)
tabler <- tableR$r
tablep <- tableR$P

# Load function to tabulate r and p values
flat_cor_mat <- function(cor_r, cor_p){
  #This function provides a simple formatting of a correlation matrix
  #into a table with 4 columns containing :
    # Column 1 : row names (variable 1 for the correlation test)
    # Column 2 : column names (variable 2 for the correlation test)
    # Column 3 : the correlation coefficients
    # Column 4 : the p-values of the correlations
  library(tidyr)
  library(tibble)
  cor_r <- rownames_to_column(as.data.frame(cor_r), var = "row")
  cor_r <- gather(cor_r, column, cor, -1)
  cor_p <- rownames_to_column(as.data.frame(cor_p), var = "row")
  cor_p <- gather(cor_p, column, p, -1)
  cor_p_matrix <- left_join(cor_r, cor_p, by = c("row", "column"))
  cor_p_matrix
}

# Tabulate 
cor <- flat_cor_mat(tableR$r, tableR$P)

# Data frame it 
cor <- as.data.frame(cor)

# Rename
names(cor) <- c("Tissue_A", "Tissue_B", "r", "p")

# Make values rounded to 2 digits 
cor$r <- round(cor$r, digits = 2)
cor$p <- format(cor$p, digits = 3, scientific = TRUE)


# Save results table 
write.csv(cor, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/Oct_2020_brain_comms/Correlation_results/smoking_cpg_spearmans_whole_group.csv", row.names = F)


############################################################################################

### CORRELATIONS FOR SCORES 

############################################################################################

# Reload datasets 

# Load in main dataset as generated in script 10 with everything in it for anlayses
data <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_10_epigenetic_signature_data_processing_part_two_phenotype_joining_for_regressions_outputs/data_withmost_recent_clinical_phenotypes_added.csv")

# Load in the data I went back and generated in script 17 for the whole LBC group 
whole <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/whole_group_LBC_blood_analysis_data_ready_for_plots.csv")

###########################################################################################

### CORRELATION PLOT

###########################################################################################
 
# Get all data in one big doc across sample regions  
corrdata <- data %>% select("LBC_ID", 
	"smokcat_recent", "pack_years", "allunitwk_recent", "BMI_recent", "HDL_recent", 
	"blood_SMO", "blood_ALC", "blood_BMI", "blood_HDL", 
	"brain_SMO", "brain_ALC", "brain_BMI", "brain_HDL", 
	"HC_SMO", "HC_ALC", "HC_BMI", "HC_HDL", 
     "BA24_SMO", "BA24_ALC", "BA24_BMI", "BA24_HDL",
      "BA46_SMO", "BA46_ALC", "BA46_BMI", "BA46_HDL", 
      "BA17_SMO", "BA17_ALC", "BA17_BMI", "BA17_HDL", 
      "BA2021_SMO", "BA2021_ALC", "BA2021_BMI", "BA2021_HDL") %>%
unique()

# Get rid of LBC identifier to run 
corrdata <- corrdata[-1]

# Sort out class of variables to run 
corrdata <- corrdata %>% mutate_if(is.factor, as.numeric) -> corrdata
corrdata <- corrdata %>% mutate_if(is.integer, as.numeric) -> corrdata
corrdata <- corrdata %>% mutate_if(is.character, as.numeric) -> corrdata

# Rename variables so it looks good in the plot and is nice and clear 
names(corrdata) <- c("smoking status trait", "pack years smoked trait", "alcohol units per week trait", "BMI trait", "HDL cholesterol trait",
 "blood smoking predictor", "blood alcohol predictor", "blood BMI predictor", "blood HDL predictor", 
 "brain smoking predictor", "brain alcohol predictor", "brain BMI predictor", "brain HDL predictor", 
 "BA35 smoking predictor", "BA35 alcohol predictor", "BA35 BMI predictor", "BA35 HDL predictor", 
     "BA24 smoking predictor", "BA24 alcohol predictor", "BA24 BMI predictor", "BA24 HDL predictor", 
     "BA46 smoking predictor", "BA46 alcohol predictor", "BA46 BMI predictor", "BA46 HDL predictor", 
     "BA17 smoking predictor", "BA17 alcohol predictor", "BA17 BMI predictor", "BA17 HDL predictor", 
     "BA2021 smoking predictor", "BA2021 alcohol predictor", "BA2021 BMI predictor", "BA2021 HDL predictor")


# Remove mean brain (no longer needed)
corrdata <- corrdata[c(-10,-11,-12,-13)]

# Use method to get p values as well as same r values as above methods
test <- rcorr.adjust(as.matrix(corrdata), type = c("spearman"), use=c("pairwise.complete.obs"))

# Extract R table
str(test)
tableR <- test$R

str(tableR)
tabler <- tableR$r
tablep <- tableR$P

# Load function to tabulate r and p values
flat_cor_mat <- function(cor_r, cor_p){
  #This function provides a simple formatting of a correlation matrix
  #into a table with 4 columns containing :
    # Column 1 : row names (variable 1 for the correlation test)
    # Column 2 : column names (variable 2 for the correlation test)
    # Column 3 : the correlation coefficients
    # Column 4 : the p-values of the correlations
  library(tidyr)
  library(tibble)
  cor_r <- rownames_to_column(as.data.frame(cor_r), var = "row")
  cor_r <- gather(cor_r, column, cor, -1)
  cor_p <- rownames_to_column(as.data.frame(cor_p), var = "row")
  cor_p <- gather(cor_p, column, p, -1)
  cor_p_matrix <- left_join(cor_r, cor_p, by = c("row", "column"))
  cor_p_matrix
}

# Tabulate 
cor <- flat_cor_mat(tableR$r, tableR$P)

# Data frame it 
cor <- as.data.frame(cor)

# Rename
names(cor) <- c("Tissue_A", "Tissue_B", "r", "p")

# Make values rounded to 2 digits 
cor$r <- round(cor$r, digits = 2)
cor$p <- round(cor$p, digits =3)

# Save results table 
write.csv(cor, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/Oct_2020_brain_comms/Correlation_results/predictor_correlations_spearmans.csv", row.names = F)

###########################################################################################

### Now in the whole LBC group 

corrdata <- whole %>% select("LBC_ID", "blood_SMO", "blood_ALC", "blood_BMI", "blood_HDL", "pack_years_w4", "allunitwk_w4", "smoking_w4", "BMI_w4", "HDL_w4") %>% unique()

corrdata <- corrdata[-1]

# Sort out class of variables to run 
corrdata <- corrdata %>% mutate_if(is.factor, as.numeric) -> corrdata
corrdata <- corrdata %>% mutate_if(is.integer, as.numeric) -> corrdata
corrdata <- corrdata %>% mutate_if(is.character, as.numeric) -> corrdata

names(corrdata) <- c("blood smoking predictor", "blood alcohol predictor", "blood BMI predictor", "blood HDL predictor", 
  "pack years smoked trait", "alcohol units per week trait", "smoking status trait", "BMI trait", "HDL cholesterol trait")

# Use method to get p values as well as same r values as above methods
test <- rcorr.adjust(as.matrix(corrdata), type = c("spearman"), use=c("pairwise.complete.obs"))

# Extract R table
str(test)
tableR <- test$R

str(tableR)
tabler <- tableR$r
tablep <- tableR$P

# Load function to tabulate r and p values
flat_cor_mat <- function(cor_r, cor_p){
  #This function provides a simple formatting of a correlation matrix
  #into a table with 4 columns containing :
    # Column 1 : row names (variable 1 for the correlation test)
    # Column 2 : column names (variable 2 for the correlation test)
    # Column 3 : the correlation coefficients
    # Column 4 : the p-values of the correlations
  library(tidyr)
  library(tibble)
  cor_r <- rownames_to_column(as.data.frame(cor_r), var = "row")
  cor_r <- gather(cor_r, column, cor, -1)
  cor_p <- rownames_to_column(as.data.frame(cor_p), var = "row")
  cor_p <- gather(cor_p, column, p, -1)
  cor_p_matrix <- left_join(cor_r, cor_p, by = c("row", "column"))
  cor_p_matrix
}

# Tabulate 
cor <- flat_cor_mat(tableR$r, tableR$P)

# Data frame it 
cor <- as.data.frame(cor)

# Rename
names(cor) <- c("Tissue_A", "Tissue_B", "r", "p")

# Make values rounded to 2 digits 
cor$r <- round(cor$r, digits = 2)
cor$p <- format(cor$p, digits =3, scientific = TRUE)
# Save results table 
write.csv(cor, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/Oct_2020_brain_comms/Correlation_results/predictor_correlations_spearmans_whole_group.csv", row.names = F)






############################################################################################


################# PEARSON VERSION 


############################################################################################

### LOAD DATA

############################################################################################

# Load in main dataset as generated in script 10 with everything in it for anlayses
data <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_10_epigenetic_signature_data_processing_part_two_phenotype_joining_for_regressions_outputs/data_withmost_recent_clinical_phenotypes_added.csv")

# Load in the data I went back and generated in script 17 for the whole LBC group 
whole <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/whole_group_LBC_blood_analysis_data_ready_for_plots.csv")


############################################################################################

### CORRELATIONS FOR SMOKING CPG

############################################################################################

### For the 14 dataset 

# calualte mean brain measure for each person 
means <- data %>% select(c("LBC_ID", "cg05575921_brain")) 
mean <- tapply(means$cg05575921_brain, means$LBC_ID, mean)
mean <- as.data.frame(mean)
names(mean)[1] <- "mean_cg05575921_brain"

# get IDs as list and rejoin
IDs <- data %>% select("LBC_ID") %>% unique() 
mean <- cbind(mean, IDs)

# Join back to original dataset
data <- merge(data, mean, by = "LBC_ID")

# Get the data I need
corrdata <- data %>% select("LBC_ID", "region", "mean_cg05575921_brain", "cg05575921_brain", "cg05575921_blood", "smokcat_recent", "pack_years") %>%
unique()

# Save blood measures 
corrblood <- corrdata[c(1,5)] %>% unique()

# Get brain measures
corrbrain <- corrdata[c(1,2,4)] %>% unique()

# Spread the data for brain so I can correlate by brain region
spread <- corrbrain %>%
  spread(key = region, value = cg05575921_brain)

# Add mean brain measure to the spread brains 
spread <- merge(spread, mean, by = "LBC_ID")

# Join the blood most recent values back in now that we have variables for each brain region
corrdata <- merge(spread, corrblood, by = "LBC_ID")

# Extract smoking info
corrsmok <- data %>% select("LBC_ID", "smokcat_recent", "pack_years") %>% unique()

# Add to main correlation dataset 
corrdata <- merge(corrdata, corrsmok, by = "LBC_ID")

# make smoking category continuous 
corrdata$smokcat_recent <- as.numeric(corrdata$smokcat_recent)

# Get rid of LBC identifier to run 
corrdata <- corrdata[-1]

# Rename it so it looks nice in the plot 
names(corrdata)[1] <- "BA17 cg05575921"
names(corrdata)[2] <- "BA20/21 cg05575921"
names(corrdata)[3] <- "BA24 cg05575921"
names(corrdata)[4] <- "BA46 cg05575921"
names(corrdata)[5] <- "BA35 cg05575921"
names(corrdata)[6] <- "Mean brain cg05575921"
names(corrdata)[7] <- "Blood cg05575921"
names(corrdata)[8] <- "Smoking status"
names(corrdata)[9] <- "Pack years smoked"

# Remove mean brain - no longer needed 
corrdata <- corrdata[-6]

###########################################################################################

# Use method to get p values as well as same r values as above methods
test <- rcorr.adjust(as.matrix(corrdata), type = c("pearson"), use=c("pairwise.complete.obs"))

# Extract R table
str(test)
tableR <- test$R

str(tableR)
tabler <- tableR$r
tablep <- tableR$P

# Load function to tabulate r and p values
flat_cor_mat <- function(cor_r, cor_p){
  #This function provides a simple formatting of a correlation matrix
  #into a table with 4 columns containing :
    # Column 1 : row names (variable 1 for the correlation test)
    # Column 2 : column names (variable 2 for the correlation test)
    # Column 3 : the correlation coefficients
    # Column 4 : the p-values of the correlations
  library(tidyr)
  library(tibble)
  cor_r <- rownames_to_column(as.data.frame(cor_r), var = "row")
  cor_r <- gather(cor_r, column, cor, -1)
  cor_p <- rownames_to_column(as.data.frame(cor_p), var = "row")
  cor_p <- gather(cor_p, column, p, -1)
  cor_p_matrix <- left_join(cor_r, cor_p, by = c("row", "column"))
  cor_p_matrix
}

# Tabulate 
cor <- flat_cor_mat(tableR$r, tableR$P)

# Data frame it 
cor <- as.data.frame(cor)

# Rename
names(cor) <- c("Tissue_A", "Tissue_B", "r", "p")

# Make values rounded to 2 digits 
cor$r <- round(cor$r, digits = 2)
cor$p <- round(cor$p, digits = 3)

# Save results table 
write.csv(cor, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/Oct_2020_brain_comms/Correlation_results/smoking_cpg_pearson.csv", row.names = F)


###########################################################################################

### For the whole group

# the 'w4' phenotype only n=500 group

# Get the data I need
corrdata <- whole %>% select("LBC_ID", "cg05575921", "smoking_w4", "pack_years_w4") %>%
unique()

# Make sure class correct 
corrdata$smoking_w4 <- as.numeric(corrdata$smoking_w4)

# Rename for results 
names(corrdata) <- c("LBC_ID", "Blood cg05575921", "Smoking status", "Pack years smoked")

# remove things not correlating
corrdata <- corrdata[-1]

# Use method to get p values as well as same r values as above methods
test <- rcorr.adjust(as.matrix(corrdata), type = c("pearson"), use=c("pairwise.complete.obs"))

# Extract R table
str(test)
tableR <- test$R

str(tableR)
tabler <- tableR$r
tablep <- tableR$P

# Load function to tabulate r and p values
flat_cor_mat <- function(cor_r, cor_p){
  #This function provides a simple formatting of a correlation matrix
  #into a table with 4 columns containing :
    # Column 1 : row names (variable 1 for the correlation test)
    # Column 2 : column names (variable 2 for the correlation test)
    # Column 3 : the correlation coefficients
    # Column 4 : the p-values of the correlations
  library(tidyr)
  library(tibble)
  cor_r <- rownames_to_column(as.data.frame(cor_r), var = "row")
  cor_r <- gather(cor_r, column, cor, -1)
  cor_p <- rownames_to_column(as.data.frame(cor_p), var = "row")
  cor_p <- gather(cor_p, column, p, -1)
  cor_p_matrix <- left_join(cor_r, cor_p, by = c("row", "column"))
  cor_p_matrix
}

# Tabulate 
cor <- flat_cor_mat(tableR$r, tableR$P)

# Data frame it 
cor <- as.data.frame(cor)

# Rename
names(cor) <- c("Tissue_A", "Tissue_B", "r", "p")

# Make values rounded to 2 digits 
cor$r <- round(cor$r, digits = 2)
cor$p <- format(cor$p, digits = 3, scientific = TRUE)


# Save results table 
write.csv(cor, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/Oct_2020_brain_comms/Correlation_results/smoking_cpg_pearson_whole_group.csv", row.names = F)


############################################################################################

### CORRELATIONS FOR SCORES 

############################################################################################

# Reload datasets 

# Load in main dataset as generated in script 10 with everything in it for anlayses
data <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_10_epigenetic_signature_data_processing_part_two_phenotype_joining_for_regressions_outputs/data_withmost_recent_clinical_phenotypes_added.csv")

# Load in the data I went back and generated in script 17 for the whole LBC group 
whole <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/whole_group_LBC_blood_analysis_data_ready_for_plots.csv")

###########################################################################################

### CORRELATION PLOT

###########################################################################################
 
# Get all data in one big doc across sample regions  
corrdata <- data %>% select("LBC_ID", 
  "smokcat_recent", "pack_years", "allunitwk_recent", "BMI_recent", "HDL_recent", 
  "blood_SMO", "blood_ALC", "blood_BMI", "blood_HDL", 
  "brain_SMO", "brain_ALC", "brain_BMI", "brain_HDL", 
  "HC_SMO", "HC_ALC", "HC_BMI", "HC_HDL", 
     "BA24_SMO", "BA24_ALC", "BA24_BMI", "BA24_HDL",
      "BA46_SMO", "BA46_ALC", "BA46_BMI", "BA46_HDL", 
      "BA17_SMO", "BA17_ALC", "BA17_BMI", "BA17_HDL", 
      "BA2021_SMO", "BA2021_ALC", "BA2021_BMI", "BA2021_HDL") %>%
unique()

# Get rid of LBC identifier to run 
corrdata <- corrdata[-1]

# Sort out class of variables to run 
corrdata <- corrdata %>% mutate_if(is.factor, as.numeric) -> corrdata
corrdata <- corrdata %>% mutate_if(is.integer, as.numeric) -> corrdata
corrdata <- corrdata %>% mutate_if(is.character, as.numeric) -> corrdata

# Rename variables so it looks good in the plot and is nice and clear 
names(corrdata) <- c("smoking status trait", "pack years smoked trait", "alcohol units per week trait", "BMI trait", "HDL cholesterol trait",
 "blood smoking predictor", "blood alcohol predictor", "blood BMI predictor", "blood HDL predictor", 
 "brain smoking predictor", "brain alcohol predictor", "brain BMI predictor", "brain HDL predictor", 
 "BA35 smoking predictor", "BA35 alcohol predictor", "BA35 BMI predictor", "BA35 HDL predictor", 
     "BA24 smoking predictor", "BA24 alcohol predictor", "BA24 BMI predictor", "BA24 HDL predictor", 
     "BA46 smoking predictor", "BA46 alcohol predictor", "BA46 BMI predictor", "BA46 HDL predictor", 
     "BA17 smoking predictor", "BA17 alcohol predictor", "BA17 BMI predictor", "BA17 HDL predictor", 
     "BA2021 smoking predictor", "BA2021 alcohol predictor", "BA2021 BMI predictor", "BA2021 HDL predictor")


# Remove mean brain (no longer needed)
corrdata <- corrdata[c(-10,-11,-12,-13)]


# Use method to get p values as well as same r values as above methods
test <- rcorr.adjust(as.matrix(corrdata), type = c("pearson"), use=c("pairwise.complete.obs"))

# Extract R table
str(test)
tableR <- test$R

str(tableR)
tabler <- tableR$r
tablep <- tableR$P

# Load function to tabulate r and p values
flat_cor_mat <- function(cor_r, cor_p){
  #This function provides a simple formatting of a correlation matrix
  #into a table with 4 columns containing :
    # Column 1 : row names (variable 1 for the correlation test)
    # Column 2 : column names (variable 2 for the correlation test)
    # Column 3 : the correlation coefficients
    # Column 4 : the p-values of the correlations
  library(tidyr)
  library(tibble)
  cor_r <- rownames_to_column(as.data.frame(cor_r), var = "row")
  cor_r <- gather(cor_r, column, cor, -1)
  cor_p <- rownames_to_column(as.data.frame(cor_p), var = "row")
  cor_p <- gather(cor_p, column, p, -1)
  cor_p_matrix <- left_join(cor_r, cor_p, by = c("row", "column"))
  cor_p_matrix
}

# Tabulate 
cor <- flat_cor_mat(tableR$r, tableR$P)

# Data frame it 
cor <- as.data.frame(cor)

# Rename
names(cor) <- c("Tissue_A", "Tissue_B", "r", "p")

# Make values rounded to 2 digits 
cor$r <- round(cor$r, digits = 2)
cor$p <- round(cor$p, digits =3)

# Save results table 
write.csv(cor, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/Oct_2020_brain_comms/Correlation_results/predictor_correlations_pearson.csv", row.names = F)

###########################################################################################

### Now in the whole LBC group 

corrdata <- whole %>% select("LBC_ID", "blood_SMO", "blood_ALC", "blood_BMI", "blood_HDL", "pack_years_w4", "allunitwk_w4", "smoking_w4", "BMI_w4", "HDL_w4") %>% unique()

corrdata <- corrdata[-1]

# Sort out class of variables to run 
corrdata <- corrdata %>% mutate_if(is.factor, as.numeric) -> corrdata
corrdata <- corrdata %>% mutate_if(is.integer, as.numeric) -> corrdata
corrdata <- corrdata %>% mutate_if(is.character, as.numeric) -> corrdata

names(corrdata) <- c("blood smoking predictor", "blood alcohol predictor", "blood BMI predictor", "blood HDL predictor", 
  "pack years smoked trait", "alcohol units per week trait", "smoking status trait", "BMI trait", "HDL cholesterol trait")

# Use method to get p values as well as same r values as above methods
test <- rcorr.adjust(as.matrix(corrdata), type = c("pearson"), use=c("pairwise.complete.obs"))

# Extract R table
str(test)
tableR <- test$R

str(tableR)
tabler <- tableR$r
tablep <- tableR$P

# Load function to tabulate r and p values
flat_cor_mat <- function(cor_r, cor_p){
  #This function provides a simple formatting of a correlation matrix
  #into a table with 4 columns containing :
    # Column 1 : row names (variable 1 for the correlation test)
    # Column 2 : column names (variable 2 for the correlation test)
    # Column 3 : the correlation coefficients
    # Column 4 : the p-values of the correlations
  library(tidyr)
  library(tibble)
  cor_r <- rownames_to_column(as.data.frame(cor_r), var = "row")
  cor_r <- gather(cor_r, column, cor, -1)
  cor_p <- rownames_to_column(as.data.frame(cor_p), var = "row")
  cor_p <- gather(cor_p, column, p, -1)
  cor_p_matrix <- left_join(cor_r, cor_p, by = c("row", "column"))
  cor_p_matrix
}

# Tabulate 
cor <- flat_cor_mat(tableR$r, tableR$P)

# Data frame it 
cor <- as.data.frame(cor)

# Rename
names(cor) <- c("Tissue_A", "Tissue_B", "r", "p")

# Make values rounded to 2 digits 
cor$r <- round(cor$r, digits = 2)
cor$p <- format(cor$p, digits =3, scientific = TRUE)
# Save results table 
write.csv(cor, file = "/Volumes/marioni-lab/Danni/Blood_brain_project/Oct_2020_brain_comms/Correlation_results/predictor_correlations_pearson_whole_group.csv", row.names = F)







###########################################################################################

### SESSION INFO

############################################################################################

# > sessionInfo()
# R version 3.6.3 (2020-02-29)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Catalina 10.15.5

# Matrix products: default
# BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

# locale:
# [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
#  [1] RcmdrMisc_2.7-0 sandwich_2.5-1  car_3.0-8       carData_3.0-4  
#  [5] ggpubr_0.4.0    pwr_1.3-0       Hmisc_4.4-0     Formula_1.2-3  
#  [9] survival_3.2-3  lattice_0.20-41 corrplot_0.84   forcats_0.5.0  
# [13] stringr_1.4.0   dplyr_1.0.0     purrr_0.3.4     readr_1.3.1    
# [17] tidyr_1.1.0     tibble_3.0.1    ggplot2_3.3.2   tidyverse_1.3.0

# loaded via a namespace (and not attached):
#  [1] nlme_3.1-148        fs_1.4.2            lubridate_1.7.9    
#  [4] RColorBrewer_1.1-2  httr_1.4.1          tools_3.6.3        
#  [7] backports_1.1.8     R6_2.4.1            rpart_4.1-15       
# [10] nortest_1.0-4       DBI_1.1.0           colorspace_1.4-1   
# [13] nnet_7.3-14         withr_2.2.0         tidyselect_1.1.0   
# [16] gridExtra_2.3       curl_4.3            compiler_3.6.3     
# [19] cli_2.0.2           rvest_0.3.5         htmlTable_2.0.0    
# [22] xml2_1.3.2          scales_1.1.1        checkmate_2.0.0    
# [25] digest_0.6.25       foreign_0.8-76      rio_0.5.16         
# [28] base64enc_0.1-3     jpeg_0.1-8.1        pkgconfig_2.0.3    
# [31] htmltools_0.5.0     dbplyr_1.4.4        htmlwidgets_1.5.1  
# [34] rlang_0.4.6         readxl_1.3.1        rstudioapi_0.11    
# [37] generics_0.0.2      zoo_1.8-8           jsonlite_1.7.0     
# [40] acepack_1.4.1       zip_2.0.4           magrittr_1.5       
# [43] Matrix_1.2-18       Rcpp_1.0.3          munsell_0.5.0      
# [46] fansi_0.4.1         abind_1.4-5         lifecycle_0.2.0    
# [49] stringi_1.4.6       MASS_7.3-51.6       grid_3.6.3         
# [52] blob_1.2.1          crayon_1.3.4        haven_2.3.1        
# [55] splines_3.6.3       hms_0.5.3           knitr_1.29         
# [58] pillar_1.4.4        ggsignif_0.6.0      reprex_0.3.0       
# [61] glue_1.4.1          latticeExtra_0.6-29 data.table_1.12.8  
# [64] modelr_0.1.8        png_0.1-7           vctrs_0.3.1        
# [67] cellranger_1.1.0    gtable_0.3.0        assertthat_0.2.1   
# [70] xfun_0.15           openxlsx_4.1.5      broom_0.5.6        
# [73] e1071_1.7-3         rstatix_0.6.0       class_7.3-17       
# [76] cluster_2.1.0       ellipsis_0.3.1     

############################################################################################












