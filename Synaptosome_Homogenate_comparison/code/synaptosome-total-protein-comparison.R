# Synaptosome total protein test code
# 8/22/24

# clear console and global environment
rm(list=ls()); gc()
options(stringsAsFactors = F)

#set working directory
setwd("C:/Users/Shelby_macdonaldlab/OneDrive - University of Pittsburgh/MacDonaldLab_Data/ASD_Development_March2023/Paper_files/Analysis_Code")

## load libraries
library(tidyverse)
library(data.table)
library(readxl)
library(dplyr)

##########################################################################
# load synaptosome prep file
syn=read_excel('Synaptosome_enrichment/input/Synaptosome_prep.xlsx')
# change Hom_ID to match meta file subject ids
names(syn)[1]="ID"

# load meta data
load("Code_Homogenate/WorkData/RAW-HOMOGENATE-PEPTIDE.RData")

# combine subject information with syn total protein data
df=merge(syn, meta.samples[,c(4:8)],by="ID") # 62 x 10

# retain only subjects that had exact known starting mg weights
# this excludes block prep 5
# excludes samples that were exactly 40 mg (these values were estimated during tissue collection)
df=df[df$Block_prep != 5, ] # 58 x 10
df=df[df$Hom_mg != 40.0, ] # 47 x 10

#############################
### create new column that is the ratio of Synaptosome ug to Homogenate mg
df$ratio=df$Syn_totalug/df$Hom_mg

# fit the linear model
#model <- lm(ratio ~ DX + AGEyr + PMIhr + Block_prep + SEX, data=df)


# calculate residuals
residuals=resid(model)
# store residuals in df
df$residuals=residuals

# add diagnosis and age back to the residuals
# estimate the effects of diagnosis and age and then add them back to the adjusted residuals

# Get the predicted values for Diagnosis and Age
diagnosis_age_effect <- predict(model, newdata = df, terms = c("DX", "AGEyr"))

# Add these effects back to the residuals
df$adjusted_residuals <- df$residuals + diagnosis_age_effect


# test for an age association

# create new levels
df$DX=factor(df$DX, levels=unique(df$DX))

# Split the data by Diagnosis
diagnosis_groups <- split(df, df$DX)

# Calculate Pearson correlation for each group
cor_results <- lapply(diagnosis_groups, function(group) {
  cor_test <- cor.test(group$adjusted_residuals, group$AGEyr,method="pearson")
  return(list(correlation = cor_test$estimate, p_value = cor_test$p.value))
})

# Convert the results into a data frame for easier viewing
cor_results_df <- do.call(rbind, cor_results)
cor_results_df <- as.data.frame(cor_results_df)
names(cor_results_df) <- c("Correlation", "P_Value")

# View the results
print(cor_results_df)

# Correlation    P_Value
# ASD   0.1111338  0.6136817
# CTL  -0.5293037 0.00782067
