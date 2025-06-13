rm(list=ls()); gc()
options(stringsAsFactors = F)

library(tidyverse)
library(readxl)
library(purrr)


setwd("C:/Users/Shelby_macdonaldlab/OneDrive - University of Pittsburgh/MacDonaldLab_Data/ASD_Development_March2023/Paper_files/Analysis_Code")

# load background proteins for Hom and Syn
# split by up and downregulation in this dataset

load("ASD_enrichment/WorkData/protein-input.RData")



######################################################################################
# load brain region transcript data from Gandal 2022
# supplementary Data 3
# data is separated by multiple columns
Gandal.list=read_excel("ASD_enrichment/input/Gandal_2022_supplementary3.xlsx",sheet="DEGene_Statistics")

# only retain relevant columns
# Define the patterns to match
patterns <- c("^ensembl_gene_id$", "^external_gene_name$", "^WholeCortex_ASD_", "^ASD_")

# Create a logical vector that matches the patterns
columns_to_retain <- grepl(paste(patterns, collapse = "|"), names(Gandal.list))

# Subset the dataframe to retain only the matching columns
Gandal.list <- Gandal.list[, columns_to_retain]
names(Gandal.list)[1]="ensg"
names(Gandal.list)[2]="gene"

# save each brain region as individual dataframes
# List of all brain regions
regions <- c("WholeCortex", "BA9", "BA44_45", "BA24", "BA4_6", "BA38", "BA20_37", "BA41_42_22", "BA3_1_2_5", "BA7", "BA39_40", "BA17")

# Initialize an empty list to store the dataframes
region_dfs <- list()

# Loop through each region to create individual dataframes
for (region in regions) {
  # Identify columns that match the current region
  if (region == "WholeCortex") {
    region_columns <- grep(paste0("^", region, "_"), names(Gandal.list), value = TRUE)
  } else {
    region_columns <- grep(paste0("^ASD_", region, "_"), names(Gandal.list), value = TRUE)
  }
  
  # Subset the dataframe to include the gene, ensg, and region-specific columns
  region_df <- Gandal.list[, c("ensg", "gene", region_columns)]
  
  # Store the dataframe in the list
  region_dfs[[region]] <- region_df
}


# edit the list so each column has the same name for FC and FDR value
# Function to rename columns
rename_columns <- function(df, new_names) {
  colnames(df) <- new_names
  return(df)
}

# Define new column names
new_column_names <- c("ensg", "gene","logFC","FDR")

# Iterate through the list and rename columns in each dataframe within the list
region_dfs <- lapply(region_dfs, rename_columns, new_names = new_column_names)


######################################################################################
##################################################################
# STEPS
# load protein DE lists
# decide if filtering transcript lists for Hom or Syn (DX) or (DX:AGE)
# also consider that syn dx:age up is more like down?

# load supplemental 3 data from Gandal 2022
# determine the background for each brain region (should be the same for everyone)

# calculate total overlap of measured proteins and transcripts per brain region
# save this list as numerical values so they can be plugged in during enrichment analysis as background

# save list of differentially expressed genes identified within each region

# adjust DE protein list to only include measured proteins within each brain region
##################################################################
# example to check list
# 
#df=regions_dfs[["BA17"]]
#df=df[df$gene %in% background, ] # 4383, 3994


# select list of all brain regions
list_of_lists=region_dfs

# make combined protein.map
protein.map=rbind(Hom.map,Syn.map)
protein.map=unique(protein.map)
########################



##### start loop
# here
###
# decide if processing for Hom up or Hom down or Syn up or Syn down
#FRAC="hom.up"
#FRAC="hom.down"

#FRAC="syn.up"
#FRAC="syn.down"

#FRAC="syn.int.up"
FRAC="syn.int.down"

# only select hom or syn
# no up/downregulation
#background=Hom.map$ensg
background=Syn.map$ensg
###

# Initialize a list to save the results
###
#results_Hom.up <- list()
#background_Hom.up <- list()
#DE_Hom.up <- list()

#results_Hom.down <- list()
#background_Hom.down <- list()
#DE_Hom.down <- list()

###
#results_Syn.up <- list()
#background_Syn.up <- list()
#DE_Syn.up <- list()

#results_Syn.down <- list()
#background_Syn.down <- list()
#DE_Syn.down <- list()

###
#results_Syn.int.up <- list()
#background_Syn.int.up <- list()
#DE_Syn.int.up <- list()

results_Syn.int.down <- list()
background_Syn.int.down <- list()
DE_Syn.int.down <- list()
##################################################################

# Loop over each list in the list of lists
for (list_name in names(list_of_lists)) {
  
  # Extract the current list
  current_list <- list_of_lists[[list_name]]
  
  # Find the overlap with the background (Hom or Syn)
  overlapping_genes <- current_list$ensg[current_list$ensg %in% background]
  
  #####################################
  # Filter DE Hom or Syn proteins for overlapping_genes only
  #filtered_proteins <- overlapping_genes[overlapping_genes %in% Hom.up.genes$ensg]
  #filtered_proteins <- overlapping_genes[overlapping_genes %in% Hom.down.genes$ensg]
  
  #filtered_proteins <- overlapping_genes[overlapping_genes %in% Syn.up.genes$ensg]
  #filtered_proteins <- overlapping_genes[overlapping_genes %in% Syn.down.genes$ensg]
  
  #filtered_proteins <- overlapping_genes[overlapping_genes %in% Syn.int.up.genes$ensg]
  filtered_proteins <- overlapping_genes[overlapping_genes %in% Syn.int.down.genes$ensg]
  
  # map filtered protein ensg names to gene names
  filtered_proteins=protein.map$gene[match(filtered_proteins, protein.map$ensg)]
  #####################################
  
  #####################################
  # Filter overlapping protein/transcript genes based on FDR < 0.1 and fold change up/downregulation
  # up
  #filtered_genes <- current_list$ensg[current_list$ensg %in% overlapping_genes & current_list$FDR < 0.1 & current_list$logFC > 0]
  
  # down
  filtered_genes <- current_list$ensg[current_list$ensg %in% overlapping_genes & current_list$FDR < 0.1 & current_list$logFC < 0]
  
  # map filtered transcript ensg names to gene names
  filtered_genes=protein.map$gene[match(filtered_genes, protein.map$ensg)]
  #####################################
  
  #####################################
  # Save the background number for each cell type protein/transcript overlap
  #background_Hom.up[[list_name]] <- as.numeric(length(overlapping_genes))
  #background_Hom.down[[list_name]] <- as.numeric(length(overlapping_genes))
  
  #background_Syn.up[[list_name]] <- as.numeric(length(overlapping_genes))
  #background_Syn.down[[list_name]] <- as.numeric(length(overlapping_genes))
  
  #background_Syn.int.up[[list_name]] <- as.numeric(length(overlapping_genes))
  background_Syn.int.down[[list_name]] <- as.numeric(length(overlapping_genes))
  #####################################
  
  #####################################
  # save the list of DE transcripts in each cell type with an FDR < 0.1
  #results_Hom.up[[list_name]] <- filtered_genes
  #results_Hom.down[[list_name]] <- filtered_genes
  
  #results_Syn.up[[list_name]] <- filtered_genes
  #results_Syn.down[[list_name]] <- filtered_genes
  
  #results_Syn.int.up[[list_name]] <- filtered_genes
  results_Syn.int.down[[list_name]] <- filtered_genes
  #####################################
  
  #####################################
  # edit the list of DE proteins to match what was identified in 
  #DE_Hom.up[[list_name]] <- filtered_proteins
  #DE_Hom.down[[list_name]] <- filtered_proteins
  
  #DE_Syn.up[[list_name]] <- filtered_proteins
  #DE_Syn.down[[list_name]] <- filtered_proteins
  
  #DE_Syn.int.up[[list_name]] <- filtered_proteins
  DE_Syn.int.down[[list_name]] <- filtered_proteins
  #####################################
}

###3 repeat loop until all data compartments are stored


save(results_Hom.up,
     results_Hom.down,
     
     results_Syn.up,
     results_Syn.down,
     
     results_Syn.int.up,
     results_Syn.int.down,
     
     background_Hom.up,
     background_Hom.down,
     
     background_Syn.up,
     background_Syn.down,
     
     background_Syn.int.up,
     background_Syn.int.down,
     
     DE_Hom.up,
     DE_Hom.down,
     
     DE_Syn.up,
     DE_Syn.down,
     
     DE_Syn.int.up,
     DE_Syn.int.down,
     
     
     file="ASD_enrichment/WorkData/transcript-regions-list.RData")



