rm(list=ls()); gc()
options(stringsAsFactors = F)

# Describe identified co-regulated modules
####################################################################
# select one term representative of each module
# combine the selected term with the module descriptors


####################################################################

# load libraries
require(data.table)           # fast version of read.table
require(readxl)               # used for reading Excel format files
#require(viridis)              # colorblind friendly color palette
require(matrixStats)          # to get statistics on matrices
library(dplyr)
library(RColorBrewer)
library(tidyverse)

library(dplyr)
library(purrr)

####################################################################
# load Supplementary Table 22, Sheet "module_proteins" and separate into the 4 datasets
# load Supplementary Table 22, Sheet "module_ORA" 
#load("Resubmission/WorkData/modules.RData")
#load("Resubmission/WorkData/module_ORA.RData")


# within each dataset of the larger list
# replace matching 'old_module_name' data subset names with 'module' name from module_map dataframe
big_list=list(HOM.YNG_ORA, HOM.OLD_ORA, SYN.YNG_ORA, SYN.OLD_ORA)
names(big_list)=c("HOM.YNG_ORA", "HOM.OLD_ORA", "SYN.YNG_ORA", "SYN.OLD_ORA")

# Iterate through each dataset in the big_list
big_list <- lapply(names(big_list), function(ds_name) {
  
  # Get the specific dataset (list of modules)
  current_dataset <- big_list[[ds_name]]
  
  # Filter the map for only this dataset's modules
  # (Assuming module_map has a 'dataset' column to ensure correct matching)
  current_map <- module_map[module_map$dataset == ds_name, ]
  
  # Create a named vector for easy lookup: c("old_name" = "new_name")
  lookup_vec <- setNames(current_map$module, current_map$old_module_name)
  
  # Identify current names and find their replacements in the lookup vector
  current_names <- names(current_dataset)
  new_names <- lookup_vec[current_names]
  
  # Apply new names (keep old name if no match is found)
  names(current_dataset) <- ifelse(is.na(new_names), current_names, new_names)
  
  return(current_dataset)
})

# Re-apply the dataset names to the top-level list
names(big_list) <- c("HOM.YNG_ORA", "HOM.OLD_ORA", "SYN.YNG_ORA", "SYN.OLD_ORA")

SYN.YNG_ORA=big_list[["SYN.YNG_ORA"]]
####################################################################

select_representative_term <- function(df, subset_name) {
  
  # 1. Define preference order
  priority_levels <- c("curated; reactome", 
                       "GO; biological processes", 
                       "GO; molecular function", 
                       "GO; cellular components",
                       "curated; kegg")
  
  # 2. Check for NULL, empty, or missing column
  if (is.null(df) || nrow(df) == 0 || !("database" %in% colnames(df))) {
    # Create NA row while maintaining original column names if possible
    col_names <- if(!is.null(df)) colnames(df) else c("database", "adjP")
    blank_row <- as.data.frame(matrix(NA, nrow = 1, ncol = length(col_names)))
    colnames(blank_row) <- col_names
    return(cbind(subset_id = subset_name, blank_row))
  }
  
  # 3. Process and Sort
  result_row <- df %>%
    # Convert database to factor to enforce your specific priority order
    mutate(database = factor(database, levels = priority_levels)) %>%
    # Sort by database priority first, THEN by smallest adjP value
    arrange(database, adjP) %>%
    # Pick the top row
    slice(1) %>%
    # Add subset name as the first column
    mutate(subset_id = subset_name, .before = 1)
  
  return(result_row)
}

# Apply to your list (replace 'my_list' with your actual list name)
# imap_dfr automatically binds the results into one dataframe
HOM.OLD_Representative <- imap_dfr(HOM.OLD_ORA, ~select_representative_term(.x, .y))

# make rownames module names
rownames(SYN.OLD_Representative)=SYN.OLD_Representative$subset_id
SYN.OLD_Representative=SYN.OLD_Representative[,-1]

#save(HOM.YNG_Representative,
 #    HOM.OLD_Representative,
  #   SYN.YNG_Representative,
   #  SYN.OLD_Representative,
    # file="Resubmission/WorkData/modules_ORA_Representative.RData")

####################################################################
####################################################################