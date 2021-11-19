############################################################################################
############################################################################################
################# SCRIPT 13 - Danni Gadd - Blood vs brain DNAm project ######################
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
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr)
library(ggplot2)

############################################################################################

### LOAD DATA

############################################################################################

# Load in main dataset as generated in script 10 with everything in it for anlayses
data <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_10_epigenetic_signature_data_processing_part_two_phenotype_joining_for_regressions_outputs/data_withmost_recent_clinical_phenotypes_added.csv")

# Load in the data I went back and generated in script 17 for the whole LBC group 
whole <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/whole_group_LBC_blood_analysis_data_ready_for_plots_499.csv")

names(whole)[12] <- "blood_SMO"
############################################################################################

### Plot 1 - blood vs brain relationships 

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


# Plots for cg05575921 pearson


P1 <- ggplot(corrdata, aes(x=BA17, y=cg05575921_blood)) +
geom_point(colour = "blue4", size = 2) +
geom_smooth(method='lm', colour = "blue4") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = 0.79)


P2 <- ggplot(corrdata, aes(x=BA24, y=cg05575921_blood)) +
geom_point(colour = "blue4", size = 2) +
geom_smooth(method='lm', colour = "blue4") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = 0.80)


P3 <- ggplot(corrdata, aes(x=BA46, y=cg05575921_blood)) +
geom_point(colour = "blue4", size = 2) +
geom_smooth(method='lm', colour = "blue4") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = 0.81, label.y = 0.95)


P4 <- ggplot(corrdata, aes(x=HC, y=cg05575921_blood)) +
geom_point(colour = "blue4", size = 2) +
geom_smooth(method='lm', colour = "blue4") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = 0.76)


names(corrdata)[3] <- "BA20"

P5 <- ggplot(corrdata, aes(x=BA20, y=cg05575921_blood)) +
geom_point(colour = "blue4", size = 2) +
geom_smooth(method='lm', colour = "blue4") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = 0.83)


ggarrange(P3, P2, P4, P5, P1,
          ncol = 5, nrow = 1)


# Now load in the smoking & HDL signature plots 


# BA35
BA35 <- data %>% select("blood_SMO","HC_SMO") %>% unique()

BA35 <- ggplot(BA35, aes(x=blood_SMO, y=HC_SMO)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = 4.5, label.y = 4.4)


# BA46
BA46 <- data %>% select("blood_SMO","BA46_SMO") %>% unique()

BA46 <- ggplot(BA46, aes(x=blood_SMO, y=BA46_SMO)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = 5, label.y = 4.2)


# BA24 
BA24 <- data %>% select("blood_SMO","BA24_SMO") %>% unique()

BA24 <- ggplot(BA24, aes(x=blood_SMO, y=BA24_SMO)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = 5, label.y = 4.45)



# BA20/21
BA2021 <- data %>% select("blood_SMO","BA2021_SMO") %>% unique()

BA2021 <- ggplot(BA2021, aes(x=blood_SMO, y=BA2021_SMO)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = 5, label.y = 4.30)



# BA17
BA17 <- data %>% select("blood_SMO","BA17_SMO") %>% unique()

BA17 <- ggplot(BA17, aes(x=blood_SMO, y=BA17_SMO)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = 5, label.y = 4.15)




### Do a plot between blood score & brain scores

# BA35
HDL_BA35 <- data %>% select("blood_HDL","HC_HDL") %>% unique()

HDL_BA35 <- ggplot(HDL_BA35, aes(x=blood_HDL, y=HC_HDL)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = 2.6, label.y =1)


# BA46
HDL_BA46 <- data %>% select("blood_HDL","BA46_HDL") %>% unique()

HDL_BA46 <- ggplot(HDL_BA46, aes(x=blood_HDL, y=BA46_HDL)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = 2.6, label.y =0.9)


# BA24 
HDL_BA24 <- data %>% select("blood_HDL","BA24_HDL") %>% unique()

HDL_BA24 <- ggplot(HDL_BA24, aes(x=blood_HDL, y=BA24_HDL))+
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = 2.3, label.y =1.06)


# BA20/21
HDL_BA2021 <- data %>% select("blood_HDL","BA2021_HDL") %>% unique()

HDL_BA2021 <- ggplot(HDL_BA2021, aes(x=blood_HDL, y=BA2021_HDL))+
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = 2.6, label.y =1.1)


# BA17
HDL_BA17 <- data %>% select("blood_HDL","BA17_HDL") %>% unique()

HDL_BA17 <- ggplot(HDL_BA17, aes(x=blood_HDL, y=BA17_HDL))+
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = 2.6, label.y =0.95)


# So for the smoking cpg 
ggarrange(P3, P2, P4, P5, P1,
          ncol = 5, nrow = 1)

# And for the signatures 
ggarrange(BA46, BA24, BA35, BA2021, BA17,
          ncol = 5, nrow = 1)

# And for the signatures 
ggarrange(HDL_BA46, BA24, BA35, BA2021, BA17,
          ncol = 5, nrow = 1)

# Trying them altogether 
figure <- ggarrange(P3, P2, P4, P5, P1,
  BA46, BA24, BA35, BA2021, BA17,
  HDL_BA46, HDL_BA24, HDL_BA35, HDL_BA2021, HDL_BA17,
          ncol = 5, nrow = 3)


# Save as PDF
#pdf("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts/Chosen_code_for_submission_to_GitHub/Figure1_blood_brain_correlations.pdf",width=15, height =7.5) 
#figure <- ggarrange(P3, P2, P4, P5, P1,
#  BA46, BA24, BA35, BA2021, BA17,
#  HDL_BA46, HDL_BA24, HDL_BA35, HDL_BA2021, HDL_BA17,
#          ncol = 5, nrow = 3)
#ggsave(filename = paste0("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts/Chosen_code_for_submission_to_GitHub/Figure1_blood_brain_correlations.pdf"), plot=figure)




# Do alcohol and BMI ones for supplementary section

# BA35
ALC_BA35 <- data %>% select("blood_ALC","HC_ALC") %>% unique()

ALC_BA35 <- ggplot(ALC_BA35, aes(x=blood_ALC, y=HC_ALC)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7,  label.x = -11.8, label.y =-8)


# BA46
ALC_BA46 <- data %>% select("blood_ALC","BA46_ALC") %>% unique()

ALC_BA46 <- ggplot(ALC_BA46, aes(x=blood_ALC, y=BA46_ALC)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = -11.8, label.y =-8)


# BA24 
ALC_BA24 <- data %>% select("blood_ALC","BA24_ALC") %>% unique()

ALC_BA24 <- ggplot(ALC_BA24, aes(x=blood_ALC, y=BA24_ALC))+
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = -11.8, label.y =-8.1)


# BA20/21
ALC_BA2021 <- data %>% select("blood_ALC","BA2021_ALC") %>% unique()

ALC_BA2021 <- ggplot(ALC_BA2021, aes(x=blood_ALC, y=BA2021_ALC))+
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = -11, label.y =-8.1)


# BA17
ALC_BA17 <- data %>% select("blood_ALC","BA17_ALC") %>% unique()

ALC_BA17 <- ggplot(ALC_BA17, aes(x=blood_ALC, y=BA17_ALC))+
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = -11.8, label.y =-8.25)


# BMI plot 

# BA35
BMI_BA35 <- data %>% select("blood_BMI","HC_BMI") %>% unique()

BMI_BA35 <- ggplot(BMI_BA35, aes(x=blood_BMI, y=HC_BMI)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7,  label.x = -0.69, label.y = )


# BA46
BMI_BA46 <- data %>% select("blood_BMI","BA46_BMI") %>% unique()

BMI_BA46 <- ggplot(BMI_BA46, aes(x=blood_BMI, y=BA46_BMI)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = -0.7, label.y =0.22)


# BA24 
BMI_BA24 <- data %>% select("blood_BMI","BA24_BMI") %>% unique()

BMI_BA24 <- ggplot(BMI_BA24, aes(x=blood_BMI, y=BA24_BMI))+
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = -0.80, label.y =)


# BA20/21
BMI_BA2021 <- data %>% select("blood_BMI","BA2021_BMI") %>% unique()

BMI_BA2021 <- ggplot(BMI_BA2021, aes(x=blood_BMI, y=BA2021_BMI))+
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = -0.86, label.y =0.19)


# BA17
BMI_BA17 <- data %>% select("blood_BMI","BA17_BMI") %>% unique()

BMI_BA17 <- ggplot(BMI_BA17, aes(x=blood_BMI, y=BA17_BMI))+
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = -0.80, label.y =)


# Trying them altogether 
figure <- ggarrange(P3, P2, P4, P5, P1,
  BA46, BA24, BA35, BA2021, BA17,
  HDL_BA46, HDL_BA24, HDL_BA35, HDL_BA2021, HDL_BA17,
  ALC_BA46, ALC_BA24, ALC_BA35, ALC_BA2021, ALC_BA17,
  BMI_BA46, BMI_BA24, BMI_BA35, BMI_BA2021, BMI_BA17,
          ncol = 5, nrow = 5)


# Save as PDF
pdf("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts/Chosen_code_for_submission_to_GitHub/Figure1_blood_brain_correlations_all_traits_PEARSON.pdf",width=16, height =13) 
figure <- ggarrange(P3, P2, P4, P5, P1,
  BA46, BA24, BA35, BA2021, BA17,
  HDL_BA46, HDL_BA24, HDL_BA35, HDL_BA2021, HDL_BA17,
  ALC_BA46, ALC_BA24, ALC_BA35, ALC_BA2021, ALC_BA17,
  BMI_BA46, BMI_BA24, BMI_BA35, BMI_BA2021, BMI_BA17,
          ncol = 5, nrow = 5)
ggsave(filename = paste0("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts/Chosen_code_for_submission_to_GitHub/Figure1_blood_brain_correlations_all_traits_PEARSON.pdf"), plot=figure)







############################################################################################

### Plot 2 - blood & brain vs trait relationships 

############################################################################################

### PACK YEARS

# Work out how many have pack years for plots 
test <- whole$pack_years_w4
test2 <- table(is.na(test))
test2 # 306

# Filter so pack years 0 removed 
#test <- whole %>% select("LBC_ID", "pack_years_w4", "blood_SMO", "cg05575921")
#test <- test %>% filter(!pack_years_w4 == 0)
# r = 0.27 with never smoke removed, wheras r = 0.31 with them in - not much different 

# Whole LBC population Blood 
LBC_a <- ggplot(whole, aes(x=blood_SMO, y=pack_years_w4)) +
geom_point(colour = "red", size = 1) +
geom_smooth(method='lm', colour = "red") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = , label.y =)


# Blood 
blood_a <- ggplot(data, aes(y=pack_years, x=blood_SMO)) +
geom_point(colour = "red", size = 2) +
geom_smooth(method='lm', colour = "red") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = 3.0, label.y =105)


# Brain regions
HC <- subset(data, region == "HC")
BA24 <- subset(data, region == "BA24")
BA17 <- subset(data, region == "BA17")
BA46 <- subset(data, region == "BA46")
BA2021 <- subset(data, region == "BA20/21")

HC_a <- ggplot(HC, aes(y=pack_years, x=HC_SMO)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = 4.0, label.y =120)


BA24_a <- ggplot(BA24, aes(y=pack_years, x=BA24_SMO)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = 4.0, label.y =125)


BA17_a <- ggplot(BA17, aes(y=pack_years, x=BA17_SMO)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = 3.5, label.y =130)


BA46_a <- ggplot(BA46, aes(y=pack_years, x=BA46_SMO)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x =3.5, label.y =115)


BA2021_a <- ggplot(BA2021, aes(y=pack_years, x=BA2021_SMO)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = 3.6, label.y =115)




###########################################################################################


### HDL

# Work out how many have pack years for plots 
test <- whole$HDL_w4
test2 <- table(is.na(test))
test2 # 483

# Whole LBC population Blood 
LBC_p <- ggplot(whole, aes(y=HDL_w4, x= blood_HDL)) +
geom_point(colour = "red", size = 1) +
geom_smooth(method='lm', colour = "red") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = , label.y =)


# Blood 
blood_p <- ggplot(data, aes(y=HDL_recent, x=blood_HDL)) +
geom_point(colour = "red", size = 2) +
geom_smooth(method='lm', colour = "red") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = , label.y =)

# Brain regions
HC <- subset(data, region == "HC")
BA24 <- subset(data, region == "BA24")
BA17 <- subset(data, region == "BA17")
BA46 <- subset(data, region == "BA46")
BA2021 <- subset(data, region == "BA20/21")

HC_p <- ggplot(HC, aes(y=HDL_recent, x=HC_HDL)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = , label.y =)


BA24_p <- ggplot(BA24, aes(y=HDL_recent, x=BA24_HDL)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x =0.77 , label.y =2)


BA17_p <- ggplot(BA17, aes(y=HDL_recent, x=BA17_HDL)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = , label.y =2.1)


BA46_p <- ggplot(BA46, aes(y=HDL_recent, x=BA46_HDL)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x =0.73 , label.y =2.1)


BA2021_p <- ggplot(BA2021, aes(y=HDL_recent, x=BA2021_HDL)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = , label.y =)




###########################################################################################

### ALCOHOL

# Work out how many have pack years for plots 
test <- whole$allunitwk_w4
test2 <- table(is.na(test))
test2 # 358

# Whole LBC population Blood 
LBC_q <- ggplot(whole, aes(y=allunitwk_w4, x= blood_ALC)) +
geom_point(colour = "red", size = 1) +
geom_smooth(method='lm', colour = "red") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = , label.y =)


# Blood 
blood_q <- ggplot(data, aes(y=allunitwk_recent, x=blood_ALC)) +
geom_point(colour = "red", size = 2) +
geom_smooth(method='lm', colour = "red") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = , label.y =)


# Brain regions
HC <- subset(data, region == "HC")
BA24 <- subset(data, region == "BA24")
BA17 <- subset(data, region == "BA17")
BA46 <- subset(data, region == "BA46")
BA2021 <- subset(data, region == "BA20/21")

HC_q <- ggplot(HC, aes(y=allunitwk_recent, x=HC_ALC)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x =-9, label.y =35)


BA24_q <- ggplot(BA24, aes(y=allunitwk_recent, x=BA24_ALC)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x =-9.3 , label.y =35)


BA17_q <- ggplot(BA17, aes(y=allunitwk_recent, x=BA17_ALC)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x =-9.1 , label.y =35)


BA46_q <- ggplot(BA46, aes(y=allunitwk_recent, x=BA46_ALC)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x =-9.5 , label.y =45)


BA2021_q <- ggplot(BA2021, aes(y=allunitwk_recent, x=BA2021_ALC)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x =-9.0 , label.y =52)



### BMI

# Work out how many have pack years for plots 
test <- whole$BMI_w4
test2 <- table(is.na(test))
test2 # 498

# Whole LBC population Blood 
LBC_r <- ggplot(whole, aes(y=BMI_w4, x= blood_BMI)) +
geom_point(colour = "red", size = 1) +
geom_smooth(method='lm', colour = "red") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = , label.y =)

# Blood 
blood_r <- ggplot(data, aes(y=BMI_recent, x=blood_BMI)) +
geom_point(colour = "red", size = 2) +
geom_smooth(method='lm', colour = "red") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = , label.y =)


# Brain regions
HC <- subset(data, region == "HC")
BA24 <- subset(data, region == "BA24")
BA17 <- subset(data, region == "BA17")
BA46 <- subset(data, region == "BA46")
BA2021 <- subset(data, region == "BA20/21")

HC_r <- ggplot(HC, aes(y=BMI_recent, x=HC_BMI)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = , label.y =)


BA24_r <- ggplot(BA24, aes(y=BMI_recent, x=BA24_BMI)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = , label.y =)


BA17_r <- ggplot(BA17, aes(y=BMI_recent, x=BA17_BMI)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x =0.08 , label.y =32)


BA46_r <- ggplot(BA46, aes(y=BMI_recent, x=BA46_BMI)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = , label.y =)


BA2021_r <- ggplot(BA2021, aes(y=BMI_recent, x=BA2021_BMI)) +
geom_point(colour = "dodgerblue", size = 2) +
geom_smooth(method='lm', colour = "dodgerblue") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = 0.04, label.y =35)





# CPG FOR SMOKING AND PACK YEARS 

### PACK YEARS SECOND 

# Load in main dataset as generated in script 10 with everything in it for anlayses
data <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_10_epigenetic_signature_data_processing_part_two_phenotype_joining_for_regressions_outputs/data_withmost_recent_clinical_phenotypes_added.csv")

# Load in the data I went back and generated in script 17 for the whole LBC group 
whole <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/whole_group_LBC_blood_analysis_data_ready_for_plots_499.csv")


names(whole)[12] <- "blood_SMO"

# Get number of people in LBC main group with w4 info 
test <- whole$pack_years_w4
table <- table(is.na(test))
#(N=306) with never smokers inlcuded, but n=37 without 
# when pack years = 0 removed r = -0.22, p = 0.19, which is similar to -0.28 when 0 kept in for never smokers - however p value becomes significant when 0's are included 

# Blood  plot
LBC_b <- ggplot(whole, aes(y=pack_years_w4, x=cg05575921)) +
geom_point(colour = "red", size = 1) +
geom_smooth(method='lm', colour = "red") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = , label.y =)

# Blood  plot
blood_b <- ggplot(data, aes(y=pack_years, x=cg05575921_blood)) +
geom_point(colour = "red", size = 2) +
geom_smooth(method='lm', colour = "red") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x = , label.y =)

# Brain plots
HC <- subset(data, region == "HC")
BA24 <- subset(data, region == "BA24")
BA17 <- subset(data, region == "BA17")
BA46 <- subset(data, region == "BA46")
BA2021 <- subset(data, region == "BA20/21")

HC_b <- ggplot(HC, aes(y=pack_years, x=cg05575921_brain)) +
geom_point(colour = "blue4", size = 2) +
geom_smooth(method='lm', colour = "blue4") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x =0.78, label.y =120)

BA24_b <- ggplot(BA24, aes(y=pack_years, x=cg05575921_brain)) +
geom_point(colour = "blue4", size = 2) +
geom_smooth(method='lm', colour = "blue4") + xlim(0.78,0.83) +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x =0.80, label.y =132)

BA17_b <- ggplot(BA17, aes(y=pack_years, x=cg05575921_brain)) +
geom_point(colour = "blue4", size = 2) +
geom_smooth(method='lm', colour = "blue4") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x =0.79, label.y =132)

BA46_b <- ggplot(BA46, aes(y=pack_years, x=cg05575921_brain)) +
geom_point(colour = "blue4", size = 2) +
geom_smooth(method='lm', colour = "blue4") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + xlim(0.8,0.86) +
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x =0.81, label.y =133)

BA2021_b <- ggplot(BA2021, aes(y=pack_years, x=cg05575921_brain)) +
geom_point(colour = "blue4", size = 2) +
geom_smooth(method='lm', colour = "blue4") +
theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12),     
      axis.title.y=element_blank(), axis.text.y=element_text(size=12)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7, label.x =0.81, label.y =133)


# Save as PDF
pdf("/Volumes/marioni-lab/Danni/Blood_brain_project/Oct_2020_brain_comms/Additional_checks_collated/Figure2_correlation_with_traits_PEARSON.pdf",width=19, height =13) 
figure <- ggarrange(LBC_b, blood_b, BA46_b, BA24_b, HC_b, BA2021_b, BA17_b,
LBC_a, blood_a, BA46_a, BA24_a, HC_a, BA2021_a, BA17_a,
LBC_p, blood_p, BA46_p, BA24_p, HC_p, BA2021_p, BA17_p,
LBC_q, blood_q, BA46_q, BA24_q, HC_q, BA2021_q, BA17_q,
LBC_r, blood_r, BA46_r, BA24_r, HC_r,BA2021_r, BA17_r,
ncol = 7, nrow = 5)
ggsave(filename = paste0("/Volumes/marioni-lab/Danni/Blood_brain_project/Oct_2020_brain_comms/Additional_checks_collated/Figure2_correlation_with_traits_PEARSON.pdf"), plot=figure)


















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
#  [1] Hmisc_4.4-0     Formula_1.2-3   survival_3.1-8  lattice_0.20-41
#  [5] corrplot_0.84   ggpubr_0.3.0    forcats_0.5.0   stringr_1.4.0  
#  [9] dplyr_0.8.5     purrr_0.3.3     readr_1.3.1     tidyr_1.0.3    
# [13] tibble_3.0.1    ggplot2_3.3.0   tidyverse_1.3.0

# loaded via a namespace (and not attached):
#  [1] httr_1.4.1          jsonlite_1.6        splines_3.6.3      
#  [4] carData_3.0-3       modelr_0.1.5        assertthat_0.2.1   
#  [7] latticeExtra_0.6-28 cellranger_1.1.0    pillar_1.4.4       
# [10] backports_1.1.5     glue_1.4.0          digest_0.6.23      
# [13] RColorBrewer_1.1-2  ggsignif_0.6.0      checkmate_2.0.0    
# [16] rvest_0.3.5         colorspace_1.4-1    cowplot_1.0.0      
# [19] htmltools_0.4.0     Matrix_1.2-18       pkgconfig_2.0.3    
# [22] broom_0.5.5         haven_2.2.0         scales_1.1.0       
# [25] openxlsx_4.1.3      rio_0.5.16          htmlTable_1.13.2   
# [28] mgcv_1.8-31         farver_2.0.1        generics_0.0.2     
# [31] car_3.0-7           ellipsis_0.3.0      withr_2.1.2        
# [34] nnet_7.3-12         cli_1.1.0           magrittr_1.5       
# [37] crayon_1.3.4        readxl_1.3.1        fs_1.4.1           
# [40] nlme_3.1-144        rstatix_0.5.0       xml2_1.2.2         
# [43] foreign_0.8-76      tools_3.6.3         data.table_1.12.6  
# [46] hms_0.5.2           lifecycle_0.2.0     munsell_0.5.0      
# [49] reprex_0.3.0        cluster_2.1.0       zip_2.0.4          
# [52] compiler_3.6.3      rlang_0.4.6         grid_3.6.3         
# [55] rstudioapi_0.10     htmlwidgets_1.5.1   labeling_0.3       
# [58] base64enc_0.1-3     gtable_0.3.0        abind_1.4-5        
# [61] DBI_1.0.0           curl_4.2            R6_2.4.1           
# [64] gridExtra_2.3       lubridate_1.7.4     knitr_1.28         
# [67] stringi_1.4.3       Rcpp_1.0.3          vctrs_0.2.4        
# [70] rpart_4.1-15        acepack_1.4.1       dbplyr_1.4.2       
# [73] tidyselect_1.0.0    xfun_0.11     



############################################################################################































