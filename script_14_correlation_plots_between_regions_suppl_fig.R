####################################################################################
####################################################################################
############################## CORR PLOTS ##########################################
####################################################################################
#################################################################################### 

library(tidyverse)

####################################################################################

### Load data 

####################################################################################

# Load in main dataset as generated in script 10 with everything in it for anlayses
data <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_10_epigenetic_signature_data_processing_part_two_phenotype_joining_for_regressions_outputs/data_withmost_recent_clinical_phenotypes_added.csv")

# Load in the data I went back and generated in script 17 for the whole LBC group 
whole <- read.csv("/Volumes/marioni-lab/Danni/Blood_brain_project/R_scripts_outputs/script_17_whole_LBC_group_analysis_for_predictive_signatures_outputs/whole_group_LBC_blood_analysis_data_ready_for_plots.csv")


####################################################################################

### CORRELATE

####################################################################################

library(corrplot)

# install.packages("ggcorrplot")
library(ggcorrplot)

# Isolate the data you want to correlate 

pred <- data[c(23:50)] %>% unique()

pred <- pred %>% select("BA17_SMO",   "BA2021_SMO", "BA24_SMO",   "BA46_SMO",  "HC_SMO",     "BA17_ALC",   "BA2021_ALC", "BA24_ALC",  
"BA46_ALC",   "HC_ALC",     "BA17_BMI",  
"BA2021_BMI", "BA24_BMI",   "BA46_BMI",   
"HC_BMI" ,    "BA17_HDL" ,  "BA2021_HDL", "BA24_HDL" ,  "BA46_HDL"  
, "HC_HDL")

# cor <- cor(pred, use = "pairwise.complete.obs")
cor <- cor(pred, method = "spearman", use = "pairwise.complete.obs")

colnames(cor) <- c("BA17 SMO",   "BA2021 SMO", "BA24 SMO",   "BA46 SMO",  "BA35 SMO",     "BA17 ALC",   "BA2021 ALC", "BA24 ALC",  
"BA46 ALC",   "BA35 ALC",     "BA17 BMI",  
"BA2021 BMI", "BA24 BMI",   "BA46 BMI",   
"BA35 BMI" ,    "BA17 HDL" ,  "BA2021 HDL", "BA24 HDL" ,  "BA46 HDL"  
, "BA35 HDL")

rownames(cor) <-  c("BA17 SMO",   "BA2021 SMO", "BA24 SMO",   "BA46 SMO",  "BA35 SMO",     "BA17 ALC",   "BA2021 ALC", "BA24 ALC",  
"BA46 ALC",   "BA35 ALC",     "BA17 BMI",  
"BA2021 BMI", "BA24 BMI",   "BA46 BMI",   
"BA35 BMI" ,    "BA17 HDL" ,  "BA2021 HDL", "BA24 HDL" ,  "BA46 HDL"  
, "BA35 HDL")


pdf(paste0("/Volumes/marioni-lab/Danni/Blood_brain_project/Oct_2020_brain_comms/Correlation_results/ggcorrplot_for_predictors_against_one_another_spearman.pdf"))
ggcorrplot(cor, 
           hc.order = FALSE,
           type = "lower", 
           lab = TRUE,
           lab_size = 2.5,
           colors = c("blue", "white", "red")) + 
theme(legend.title = element_blank(), text = element_text(size = 12), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12))
dev.off()



# cor <- cor(pred, use = "pairwise.complete.obs")
cor <- cor(pred, method = "pearson", use = "pairwise.complete.obs")

colnames(cor) <- c("BA17 SMO",   "BA2021 SMO", "BA24 SMO",   "BA46 SMO",  "BA35 SMO",     "BA17 ALC",   "BA2021 ALC", "BA24 ALC",  
"BA46 ALC",   "BA35 ALC",     "BA17 BMI",  
"BA2021 BMI", "BA24 BMI",   "BA46 BMI",   
"BA35 BMI" ,    "BA17 HDL" ,  "BA2021 HDL", "BA24 HDL" ,  "BA46 HDL"  
, "BA35 HDL")

rownames(cor) <-  c("BA17 SMO",   "BA2021 SMO", "BA24 SMO",   "BA46 SMO",  "BA35 SMO",     "BA17 ALC",   "BA2021 ALC", "BA24 ALC",  
"BA46 ALC",   "BA35 ALC",     "BA17 BMI",  
"BA2021 BMI", "BA24 BMI",   "BA46 BMI",   
"BA35 BMI" ,    "BA17 HDL" ,  "BA2021 HDL", "BA24 HDL" ,  "BA46 HDL"  
, "BA35 HDL")


pdf(paste0("/Volumes/marioni-lab/Danni/Blood_brain_project/Oct_2020_brain_comms/Correlation_results/ggcorrplot_for_predictors_against_one_another_pearson.pdf"))
ggcorrplot(cor, 
           hc.order = FALSE,
           type = "lower", 
           lab = TRUE,
           lab_size = 2.5,
           colors = c("blue", "white", "red")) + 
theme(legend.title = element_blank(), text = element_text(size = 12), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12))
dev.off()





### CPG SMOKING 

pred <- data[c(3,9,12)] %>% unique()

spread <- pred %>% spread("region", "cg05575921_brain", "cg05575921_blood")

spread <- spread[-1] %>% as.data.frame()
colnames(spread)[2] <- "BA2021"
colnames(spread)[5] <- "BA35"

spread$BA17 <- as.numeric(spread$BA17)
spread$BA35 <- as.numeric(spread$BA35)
spread$BA24 <- as.numeric(spread$BA24)
spread$BA46 <- as.numeric(spread$BA46)
spread$BA2021 <- as.numeric(spread$BA2021)


cor <- cor(spread, method = "spearman", use = "pairwise.complete.obs")


pdf(paste0("/Volumes/marioni-lab/Danni/Blood_brain_project/Oct_2020_brain_comms/Correlation_results/ggcorrplot_for_smoking_cpg_across_regions_spearman.pdf"))
ggcorrplot(cor, 
           hc.order = FALSE,
           type = "lower", 
           lab = TRUE,
           lab_size = 10,
           colors = c("blue", "white", "red")) + 
theme(legend.title = element_blank(), text = element_text(size = 25), axis.text.x=element_text(size=25), axis.text.y=element_text(size=25))
dev.off()


cor <- cor(spread, method = "pearson", use = "pairwise.complete.obs")

pdf(paste0("/Volumes/marioni-lab/Danni/Blood_brain_project/Oct_2020_brain_comms/Correlation_results/ggcorrplot_for_smoking_cpg_across_regions_pearson.pdf"))
ggcorrplot(cor, 
           hc.order = FALSE,
           type = "lower", 
           lab = TRUE,
           lab_size = 10,
           colors = c("blue", "white", "red")) + 
theme(legend.title = element_blank(), text = element_text(size = 25), axis.text.x=element_text(size=25), axis.text.y=element_text(size=25))
dev.off()
























