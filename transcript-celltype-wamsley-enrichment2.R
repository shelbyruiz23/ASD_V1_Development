rm(list=ls()); gc()
options(stringsAsFactors = F)

library(ggplot2)
library(rjson)
library(knitr)
library(pheatmap)
library(openxlsx)
######################################
#Run enrichment tests with
# DE (q/FDR < 0.1) protein and transcripts (Wamsley 2024)
# split by up and downregulation

# setwd
setwd("C:/Users/Shelby_macdonaldlab/OneDrive - University of Pittsburgh/MacDonaldLab_Data/ASD_Development_March2023/Paper_files/Analysis_Code")

# Run Function
### run genelistoverlap.R function

######################################
# create names of all enrichment tests

######################################
########
########
########
# load each list for which there are total different backgrounds of total overlapping proteins and transcripts
# (1) Homogenate - 36 lists
# (2) Synaptosome - 36 lists

# load workdata of cell lists
load("ASD_enrichment/WorkData/transcript-celltype-Wamsley-list2.RData")
########
########
########
######################################

########################
#####
# SELECT 1
#####
############
# make lists of proteins and genes

# geneclasses are all cell types (n=36)
#geneclasses <- results_Hom.up
#geneclasses <- results_Hom.down
#geneclasses <- results_Syn.up
geneclasses <- results_Syn.down

# proteinclasses should be all protein types
# which is just hom and syn
#proteinclasses <- DE_Hom.up
#proteinclasses <- DE_Hom.down
#proteinclasses <- DE_Syn.up
proteinclasses <- DE_Syn.down

# assign names to gene lists
geneclassnames=names(geneclasses)
geneclassnames_abbrev=geneclassnames


###########################
###########################
# run 2x2 contingency analysis with each combination of geneclasses and Hom upregulated or downregulated proteins
#######


########################
#####
# SELECT 1
#####
############
## define result directory
#result_dir="ASD_enrichment/Results/transcript_Wamsley_hom_up" # Hom_up
#result_dir="ASD_enrichment/Results/transcript_Wamsley_hom_down" # Hom_down

#result_dir="ASD_enrichment/Results/transcript_Wamsley_syn_up" # Syn_up
result_dir="ASD_enrichment/Results/transcript_Wamsley_syn_down" # Syn_down

# create result directory folder
if (!dir.exists(result_dir)) {
  dir.create(result_dir, recursive = TRUE)
}
#####
############

########################
#####
# SELECT 1
#####

#####
############

##################################
# EDITING

########################
#####
# SELECT 1
#####
############
# Initialize matrices to store results
#ORmat <- matrix(nrow = length(geneclasses), ncol = 1, dimnames = list(names(geneclasses), as.character("Hom_up")))
#logPmat <- matrix(nrow = length(geneclasses), ncol = 1, dimnames = list(names(geneclasses), as.character("Hom_up")))
#Pmat <- matrix(nrow = length(geneclasses), ncol = 1, dimnames = list(names(geneclasses), as.character("Hom_up")))

#ORmat <- matrix(nrow = length(geneclasses), ncol = 1, dimnames = list(names(geneclasses), as.character("Hom_down")))
#logPmat <- matrix(nrow = length(geneclasses), ncol = 1, dimnames = list(names(geneclasses), as.character("Hom_down")))
#Pmat <- matrix(nrow = length(geneclasses), ncol = 1, dimnames = list(names(geneclasses), as.character("Hom_down")))

#ORmat <- matrix(nrow = length(geneclasses), ncol = 1, dimnames = list(names(geneclasses), as.character("Syn_up")))
#logPmat <- matrix(nrow = length(geneclasses), ncol = 1, dimnames = list(names(geneclasses), as.character("Syn_up")))
#Pmat <- matrix(nrow = length(geneclasses), ncol = 1, dimnames = list(names(geneclasses), as.character("Syn_up")))

ORmat <- matrix(nrow = length(geneclasses), ncol = 1, dimnames = list(names(geneclasses), as.character("Syn_down")))
logPmat <- matrix(nrow = length(geneclasses), ncol = 1, dimnames = list(names(geneclasses), as.character("Syn_down")))
Pmat <- matrix(nrow = length(geneclasses), ncol = 1, dimnames = list(names(geneclasses), as.character("Syn_down")))
#####
############

##################################
# EDITING

# Loop through all elements in geneclasses and proteinclasses
# i is equal to cell type

for (i in seq_along(geneclasses)) {
  #for (j in seq_along(proteinclasses)) {
  # Intersect with background list
  genes2 <- geneclasses[[i]]
  genes2 <- data.frame(genes2)
  
  # same, regardless of directionality
  # but remain consistent
  #backgroundTotal = background_Hom.up[[i]]
  #backgroundTotal = background_Hom.down[[i]]
  
  #backgroundTotal = background_Syn.up[[i]]
  backgroundTotal = background_Syn.down[[i]]
  
  overlap_res <- genelistOverlap(data.frame(proteinclasses[[i]]),
                                 genes2,
                                 backgroundTotal,
                                 print_result = FALSE,
                                 header = FALSE)
  ORmat[i] <- overlap_res[[1]]$OR
  logPmat[i] <- -log10(overlap_res[[1]]$hypergeo_p)
  Pmat[i] <- overlap_res[[1]]$hypergeo_p
  
  
  # Write overlapping genes to a file
  #write(overlap_res[[1]]$overlapping_genes,
   #    file = file.path(result_dir, sprintf("overlap_%s_%s.txt", names(geneclasses)[i],"Hom_up")))
  
  #write(overlap_res[[1]]$overlapping_genes,
   #    file = file.path(result_dir, sprintf("overlap_%s_%s.txt", names(geneclasses)[i],"Hom_down")))
  
  # Write overlapping genes to a file
  #write(overlap_res[[1]]$overlapping_genes,
   #     file = file.path(result_dir, sprintf("overlap_%s_%s.txt", names(geneclasses)[i],"Syn_up")))
  
  # Write overlapping genes to a file
  write(overlap_res[[1]]$overlapping_genes,
       file = file.path(result_dir, sprintf("overlap_%s_%s.txt", names(geneclasses)[i],"Syn_down")))
  #}
}

# Store workdata (end of each run)
save_path <- file.path(result_dir, "overlap_results.RData")
save(ORmat, logPmat, Pmat, file = save_path)




######################################################################################################
######################################################################################################
######################################################################################################

# after done looping through all conditions
# save overlapping genes to one dataframe 
### export overlapping genes of proteins and transcripts for supplemental table
# clear everything
# load overlapping genes from each comparison and export as excel sheet
rm(list=ls()); gc()

library(purrr)

# define folder to look at
main_directory="ASD_enrichment/Results"

# define the path for every folder in result_directories that begins with "transcript_Gandal"
# hom
subfolders=list.files(path = main_directory, pattern = "^transcript_Wamsley_hom_", full.names = TRUE)

# syn
#subfolders=list.files(path = main_directory, pattern = "^transcript_Wamsley_syn_", full.names = TRUE)

# Initialize an empty list to store file paths
all_text_files <- list()

# Loop through each relevant subfolder to get text files
for (subfolder in subfolders) {
  text_files <- list.files(path = subfolder, pattern = "\\.txt$", full.names = TRUE)
  all_text_files <- c(all_text_files, text_files)
}

# name lists by their ending file names
names(all_text_files)=all_text_files

# Read the contents of each text file
text_contents <- map(all_text_files, ~ readLines(.x, warn = FALSE))

# remove path from names
shortened_names <- sub("^.*(overlap_.*)", "\\1", all_text_files)
names(text_contents) <- shortened_names

#modify names
names(text_contents) <- sub("^.*overlap_", "", names(text_contents))
#names(text_contents) <- sub("\\.genes", "", names(text_contents))
names(text_contents) <- sub("\\.txt", "", names(text_contents))
#names(text_contents) <- sub("\\.genes", "", names(text_contents))


# Convert list of lists to list of dataframes
list_of_dataframes <- map(text_contents, function(content) {
  df <- as.data.frame(content)
  colnames(df) <- "gene"
  return(df)
})



# Add a value of NA if the list is empty
list_of_dataframes <- map(list_of_dataframes, function(df) {
  df$gene[df$gene == ""] <- NA
  return(df)
})


# Combine overlap lists into one dataframe

# Function 
# add column value to each dataframe for the corresponding comparison
add_source_column <- function(df, source) {
  df %>%
    mutate(comparison = source)
}

# Combine all dataframes into one and add a source column
combined_df <- map2_dfr(list_of_dataframes, names(list_of_dataframes), add_source_column)

# save the combined dataframe as supplemental table
#write.csv(combined_df, file="ASD_enrichment/Results/overlap/transcript-Wamsley-overlap-list2_hom.csv")
#write.csv(combined_df, file="ASD_enrichment/Results/overlap/transcript-Wamsley-overlap-list2_syn.csv")



######################################################################################################
######################################################################################################
######################################################################################################
