rm(list=ls()); gc()
options(stringsAsFactors = F)

library(tidyverse)
library(readxl)
library(purrr)


setwd("C:/Users/Shelby_macdonaldlab/OneDrive - University of Pittsburgh/MacDonaldLab_Data/ASD_Development_March2023/Paper_files/Analysis_Code")

# load background proteins for Hom and Syn
# split by up and downregulation in this dataset

load("ASD_enrichment/WorkData/protein-input.RData")

###
# decide if processing for Hom up or Hom down or Syn up or Syn down
#FRAC="hom.up"
#FRAC="hom.down"

#FRAC="syn.up"
FRAC="syn.down"

# only select hom or syn
# no up/downregulation
#background=Hom.map$gene 
background=Syn.map$gene
###

# load cell-type specific transcript data from Wamsley 2024
# data is separated by multiple sheets

# load neuronal clusters or glia clusters
file_path1="ASD_enrichment/input/Wamsley_S6A_Neuronalclusters_DE_genes_2_21_23.xlsx"
file_path2="ASD_enrichment/input/Wamsley_S6B_Gliaclusters_DE_genes_2_21_23.xlsx"

# save each sheet as a dataframe
sheet_names <- excel_sheets(file_path1)
Neuronal <- map(sheet_names, ~ read_excel(file_path1, sheet = .x))
names(Neuronal) <- sub("Complete_geneList_", "", sheet_names)

sheet_names <- excel_sheets(file_path2)
Glia <- map(sheet_names, ~ read_excel(file_path2, sheet = .x))
names(Glia) <- sub("Complete_geneList_", "", sheet_names)

# paste all cell types together
tmp=paste0(c(names(Neuronal),names(Glia))) # 36


##################################################################
# STEPS
# load protein DE lists
# decide if filtering transcript lists for Hom or Syn
# load supplemental 6A and 6B data from Wamsley 2024 (downloaded from Biorxv because the data sheet was unavailable on science)
# determine the background for each cell cluster (n=36 cells, paper keeps referring 35 cell types?)
# calculate total overlap of measured proteins and transcripts per cell cluster
# save this list as numerical values so they can be plugged in during enrichment analysis as background
# save list of differentially expressed genes identified within each cell type

# but also need to adjust DE protein list for the number of measured proteins within each cell type
##################################################################
# example to check list
# 
#df=Neuronal[["INT_5_SST"]]
#df=df[df$gene %in% background, ] # 4383, 3994


##### loop
#
# select neuron and glia lists
list_of_lists=c(Neuronal, Glia)


# Initialize a list to save the results
#results_Hom.up <- list()
#background_Hom.up <- list()
#DE_Hom.up <- list()

#results_Hom.down <- list()
#background_Hom.down <- list()
#DE_Hom.down <- list()

#results_Syn.up <- list()
#background_Syn.up <- list()
#DE_Syn.up <- list()

results_Syn.down <- list()
background_Syn.down <- list()
DE_Syn.down <- list()

# Loop over each list in the list of lists
for (list_name in names(list_of_lists)) {
  
  # Extract the current list
  current_list <- list_of_lists[[list_name]]
  
  # Find the overlap with the background (Hom or Syn)
  overlapping_genes <- current_list$gene[current_list$gene %in% background]
  
  # Filter DE Hom or Syn proteins for overlapping_genes only
  #filtered_proteins <- overlapping_genes[overlapping_genes %in% Hom.up.genes$gene]
  #filtered_proteins <- overlapping_genes[overlapping_genes %in% Hom.down.genes$gene]
  
  #filtered_proteins <- overlapping_genes[overlapping_genes %in% Syn.up.genes$gene]
  filtered_proteins <- overlapping_genes[overlapping_genes %in% Syn.down.genes$gene]
  
  # Filter overlapping protein/transcript genes based on FDR < 0.1 and fold change up/downregulation
  # up
  #filtered_genes <- current_list$gene[current_list$gene %in% overlapping_genes & current_list$FDR_pVal_ASD_vs_CTL < 0.1 & current_list$logFC_ASD_vs_CTL > 0]
  
  # down
  filtered_genes <- current_list$gene[current_list$gene %in% overlapping_genes & current_list$FDR_pVal_ASD_vs_CTL < 0.1 & current_list$logFC_ASD_vs_CTL < 0]
  

  # Save the background number for each cell type protein/transcript overlap
  #background_Hom.up[[list_name]] <- as.numeric(length(overlapping_genes))
  #background_Hom.down[[list_name]] <- as.numeric(length(overlapping_genes))
  
  #background_Syn.up[[list_name]] <- as.numeric(length(overlapping_genes))
  background_Syn.down[[list_name]] <- as.numeric(length(overlapping_genes))
  
  # save the list of DE transcripts in each cell type with an FDR < 0.1
  #results_Hom.up[[list_name]] <- filtered_genes
  #results_Hom.down[[list_name]] <- filtered_genes
  
  #results_Syn.up[[list_name]] <- filtered_genes
  results_Syn.down[[list_name]] <- filtered_genes
  
  # edit the list of DE proteins to match what was identified in 
  #DE_Hom.up[[list_name]] <- filtered_proteins
  #DE_Hom.down[[list_name]] <- filtered_proteins
  
  #DE_Syn.up[[list_name]] <- filtered_proteins
  DE_Syn.down[[list_name]] <- filtered_proteins
  
}

# repeat for Hom.down and Syn.up and Syn.down

#####

save(results_Hom.up,
     results_Hom.down,
     
     results_Syn.up,
     results_Syn.down,
     
     background_Hom.up,
     background_Hom.down,
     
     background_Syn.up,
     background_Syn.down,
     
     DE_Hom.up,
     DE_Hom.down,
     
     DE_Syn.up,
     DE_Syn.down,
     
     
     file="ASD_enrichment/WorkData/transcript-celltype-Wamsley-list2.RData")
