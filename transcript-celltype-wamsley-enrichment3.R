rm(list=ls()); gc()
options(stringsAsFactors = F)

library(ggplot2)
library(rjson)
library(knitr)
library(pheatmap)
library(openxlsx)


######################################
#Run enrichment tests with
# Upregulated and Downregulated protein and transcripts (Gandal 2022)

# setwd
setwd("C:/Users/Shelby_macdonaldlab/OneDrive - University of Pittsburgh/MacDonaldLab_Data/ASD_Development_March2023/Paper_files/Analysis_Code")

######################################
# load all subfolders to RData saved as OR, Pmat 

# Define the main folder and the pattern to match subfolders
main_folder <- "ASD_enrichment/Results"
pattern <- "transcript_Wamsley_"

# Get the list of subfolders that match the pattern
subfolders <- list.dirs(main_folder, recursive = FALSE, full.names = TRUE)
subfolders <- subfolders[grepl(pattern, basename(subfolders))]

# Iterate through each matching subfolder
for (subfolder in subfolders) {
  # Extract the suffix from the subfolder name
  suffix <- sub("\\Qtranscript_Wamsley_\\E", "", basename(subfolder))
  
  # Load the RData file
  rdata_file <- file.path(subfolder, "overlap_results.RData")
  load(rdata_file)
  
  # Rename the matrices
  assign(paste0(suffix, "_ORmat"), ORmat)
  assign(paste0(suffix, "_Pmat"), Pmat)
  assign(paste0(suffix, "_logPmat"), logPmat)
  
  # Clean up the original matrices from the environment
  rm(ORmat, Pmat, logPmat)
}

##############
# remove unwanted items from the global environment
########### # remove if checking up and downregulation
rm(hom_ORmat)
rm(hom_Pmat)
rm(hom_logPmat)

rm(syn_ORmat)
rm(syn_Pmat)
rm(syn_logPmat)

###########
########### # only remove if checking just Hom and Syn
#rm(hom_up_ORmat)
#rm(hom_up_Pmat)
#rm(hom_up_logPmat)

#rm(hom_down_ORmat)
#rm(hom_down_Pmat)
#rm(hom_down_logPmat)

#rm(syn_up_ORmat)
#rm(syn_up_Pmat)
#rm(syn_up_logPmat)

#rm(syn_down_ORmat)
#rm(syn_down_Pmat)
#rm(syn_down_logPmat)

##############
# make one df for logPmat

# Get all objects in the global environment
objects <- ls()

# Filter objects that end with "_logPmat"
logPmat_objects <- objects[grep("_logPmat$", objects)]

# Combine the columns of all _logPmat matrices into one dataframe
all_logPmat <- do.call(cbind, lapply(logPmat_objects, function(x) get(x)))


##############
# make one df for Pmat

# Get all objects in the global environment
objects <- ls()

# Filter objects that end with "_logPmat"
Pmat_objects <- objects[grep("_Pmat$", objects)]

# Combine the columns of all _logPmat matrices into one dataframe
all_Pmat <- do.call(cbind, lapply(Pmat_objects, function(x) get(x)))

##############
# make one df for ORmat

# Get all objects in the global environment
objects <- ls()

# Filter objects that end with "_logPmat"
ORmat_objects <- objects[grep("_ORmat$", objects)]

# Combine the columns of all _logPmat matrices into one dataframe
all_ORmat <- do.call(cbind, lapply(ORmat_objects, function(x) get(x)))

##############
# arrange the order of columns
print(colnames(all_ORmat))
#order=c("Hom.up.genes","Hom.down.genes","Syn.up.genes","Syn.down.genes","Syn.int.up.genes","Syn.int.down.genes")
order=c("Hom_up","Hom_down","Syn_up","Syn_down")
#order=c("Hom","Syn")


# fix column order for dfs
all_Pmat=all_Pmat[,order]
all_logPmat=all_logPmat[,order]
all_ORmat=all_ORmat[,order]

######################################################################
######
######
######
# create the total FDR table of all comparisons
# This will include 
# all brain regions from Gandal 2022 for transcripts split by upregulation and downregulation
# all DE proteins from Hom and Syn (DX) and Syn (DX:AGE)

# create empty FDR matrix
all_FDRmat <- matrix(nrow = length(rownames(all_Pmat)), ncol = length(colnames(all_Pmat)), dimnames = list(rownames(all_Pmat), colnames(all_Pmat)))

# compute FDR
for (i in 1:dim(all_Pmat)[2]){
  all_FDRmat[,i] = p.adjust(all_Pmat[,i], method = "fdr")
} # for (i in 1:dim(Pmat)[2]){

kable(all_FDRmat)


# define the order of the cell types
order=c("ASTRO_1",
        "ASTRO_2",
        "ASTRO_3",
        
        "BBB_Endo_1",
        "BBB_SM(Mural)",
        "BBB_MG3",
        "BBB_Pericytes",
        "MG1",
        "MG2",
        
        "EXT_3_L23",
        "EXT_6_L23",
        "EXT_10_L34",
        "EXT_11_L34",
        
        "EXT_1_L4",
        "EXT_13_L45",
        
        "EXT_7_L2356",
        "EXT_8_L6",
        
        "EXT_14_L5B",
        "EXT_12_L56",
        "EXT_2_L56",
        "EXT_4_L56",
        "EXT_5_L56",
        
        "EXT_9_L6",
        
        "INT_1_KIT",
        "INT_2_PV_SULF",
        "INT_3_RELN_ND",
        "INT_4_VIP_CR",
        "INT_5_SST",
        "INT_6_PV",
        
        "ODC_1",
        "ODC_2",
        "ODC_3",
        
        "OPCS_1",
        "OPCS_2",
        "OPCS_3",
        "OPCS_4"
        
)

# reorder all_logPmat
all_logPmat2=all_logPmat[order, ]
all_ORmat2=all_ORmat[order, ]

# make figure
pheatmap(all_logPmat2, display_numbers = round(all_ORmat2,digits=2), 
         number_color = "black", fontsize_number = 12,
         show_rownames=TRUE,
         labels_col = c("Hom_up","Hom_down","Syn_up","Syn_down"),
         color = colorRampPalette(c('white','red'))(100),
         cluster_rows = FALSE, cluster_cols = FALSE,
         breaks= seq(1.3,6, length=100))

######################################################################
######
######
######
library(purrr)

# combine stats outputs
# convert list of lists to list of dataframes
stats=list(all_ORmat = all_ORmat,
           all_Pmat = all_Pmat,
           all_logPmat = all_logPmat,
           all_FDRmat = all_FDRmat)




wb <- createWorkbook()

for (name in names(stats)) {
  addWorksheet(wb, sheetName = name)
  writeData(wb, sheet = name, stats[[name]], rowNames = TRUE)
}


#saveWorkbook(wb, "ASD_enrichment/Results/overlap/transcript-wamsley-enrichment-nodirection-stats.xlsx", overwrite = TRUE)
# tmp
#saveWorkbook(wb, "ASD_enrichment/Results/overlap/transcript-wamsley-enrichment-nodirection-stats-tmp.xlsx", overwrite = TRUE)
######
######
######
######################################################################

######################################################################
######
######
######
# plot data
library(tibble)
library(tidyverse)

######################################
# plot enrichment p values of enrichment
logPmat=as.data.frame(all_logPmat)
#####
# plot -log10(p value)

# convert rownames to columns

logPmat=rownames_to_column(logPmat)
names(logPmat)[1]="transcript"

logPmat=logPmat %>%
  pivot_longer(!transcript, names_to = "protein", values_to = "log10p")

#############
#############
# hom and syn only
# just use up name
#up=logPmat
#up$protein=factor(up$protein,levels=c("Syn","Hom"))
#############
#############

# make  up specific list
up=logPmat[grep("_up", logPmat$protein), ]


# make  down specific list
down=logPmat[grep("_down", logPmat$protein), ]


##########################################################
# add in fdr threshold
#########################
# create line on plot for designating significance cutoffs
# calculate p value cutoff for an FDR of < 0.05

# Hom
p_val=all_Pmat[c(1:36),1]
fdr_values <- p.adjust(p_val, method = "BH")
threshold_p_value <- max(p_val[fdr_values < 0.05], na.rm = TRUE) # 0.000598
threshold_logp_value=-log10(threshold_p_value) # 3.223

# Syn
p_val=all_Pmat[c(1:36),3]
fdr_values <- p.adjust(p_val, method = "BH")
threshold_p_value <- max(p_val[fdr_values < 0.05], na.rm = TRUE) # 0.0017
threshold_logp_value=-log10(threshold_p_value) # 2.76

# use homogenate FDR threshold
# as long as plot does not mislead syn findings

##########################################################
# Create a horizontal bar plot
up$protein=factor(up$protein,levels=c("Syn_up","Hom_up"))

up$transcript=factor(up$transcript,levels=c(unique(up$transcript)))

order=c("ASTRO_1",
        "ASTRO_2",
        "ASTRO_3",
        
        "BBB_Endo_1",
        "BBB_SM(Mural)",
        "BBB_MG3",
        "BBB_Pericytes",
        "MG1",
        "MG2",
        
        "EXT_3_L23",
        "EXT_6_L23",
        "EXT_10_L34",
        "EXT_11_L34",
        
        "EXT_1_L4",
        "EXT_13_L45",
        
        "EXT_7_L2356",
        "EXT_8_L6",
        
        "EXT_14_L5B",
        "EXT_12_L56",
        "EXT_2_L56",
        "EXT_4_L56",
        "EXT_5_L56",
        
        "EXT_9_L6",
        
        "INT_1_KIT",
        "INT_2_PV_SULF",
        "INT_3_RELN_ND",
        "INT_4_VIP_CR",
        "INT_5_SST",
        "INT_6_PV",
        
        "ODC_1",
        "ODC_2",
        "ODC_3",
        
        "OPCS_1",
        "OPCS_2",
        "OPCS_3",
        "OPCS_4"
        
        )


rev_order=rev(order)


ggplot(up, aes(x=log10p, y=factor(transcript, levels=rev_order), fill=protein)) + 
  geom_bar(stat="identity", position=position_dodge(),width=0.75)+
  scale_fill_manual(values=c("#50A8EB","#C909EB"))+
  theme_classic()+
  geom_vline(xintercept = 1.30103, linetype = "dashed", color = "red")+
  geom_vline(xintercept = 2.55, linetype = "dashed", color = "red")+
  xlim(0,4)+
  theme(legend.position = "none")

##########################################################
# Create a horizontal bar plot
down$protein=factor(down$protein,levels=c("Syn_down","Hom_down"))

down$transcript=factor(down$transcript,levels=c(unique(down$transcript)))

order=c("ASTRO_1",
        "ASTRO_2",
        "ASTRO_3",
        
        "BBB_Endo_1",
        "BBB_SM(Mural)",
        "BBB_MG3",
        "BBB_Pericytes",
        "MG1",
        "MG2",
        
        "EXT_3_L23",
        "EXT_6_L23",
        "EXT_10_L34",
        "EXT_11_L34",
        
        "EXT_1_L4",
        "EXT_13_L45",
        
        "EXT_7_L2356",
        "EXT_8_L6",
        
        "EXT_14_L5B",
        "EXT_12_L56",
        "EXT_2_L56",
        "EXT_4_L56",
        "EXT_5_L56",
        
        "EXT_9_L6",
        
        "INT_1_KIT",
        "INT_2_PV_SULF",
        "INT_3_RELN_ND",
        "INT_4_VIP_CR",
        "INT_5_SST",
        "INT_6_PV",
        
        "ODC_1",
        "ODC_2",
        "ODC_3",
        
        "OPCS_1",
        "OPCS_2",
        "OPCS_3",
        "OPCS_4"
        
)


rev_order=rev(order)

ggplot(down, aes(x=log10p, y=factor(transcript, levels=rev_order), fill=protein)) + 
  geom_bar(stat="identity", position=position_dodge(),width=0.75)+
  scale_fill_manual(values=c("#50A8EB","#C909EB"))+
  theme_classic()+
  geom_vline(xintercept = 1.30103, linetype = "dashed", color = "red")+
  geom_vline(xintercept = 2.55, linetype = "dashed", color = "red")+
  xlim(0,4)+
  theme(legend.position = "none")


##########################################################



######
######
######
######################################################################



