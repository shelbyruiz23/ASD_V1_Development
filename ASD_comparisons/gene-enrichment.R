rm(list=ls()); gc()
options(stringsAsFactors = F)

library(ggplot2)
library(rjson)
library(knitr)
library(pheatmap)
library(openxlsx)


######################################
#Run enrichment tests with
# ASD risk genes
# 185 ASD genes (Fu 2022)
# 185 ASD genes + 373 DD genes (unique only) (Fu 2022)
# SFARI (S,1,2)

# setwd
setwd("C:/Users/Shelby_macdonaldlab/OneDrive - University of Pittsburgh/MacDonaldLab_Data/ASD_Development_March2023/Paper_files/Analysis_Code")

# Run Function
### run genelistoverlap.R function

######################################

######################################
# create names of all enrichment tests

######################################
########
########
########
# load each list for which there are total different backgrounds of total overlapping proteins and transcripts
# (1) Homogenate - 3 lists
# (2) Synaptosome - 3 lists

# load workdata of cell lists
load("ASD_enrichment/WorkData/gene-list.RData")
########
########
########


###########################
###########################
# run 2x2 contingency analysis with each set of genes and proteins
#######


############
## define result directory
result_dir="ASD_enrichment/Results/gene-enrichment"


# create result directory folder
if (!dir.exists(result_dir)) {
  dir.create(result_dir, recursive = TRUE)
}
#####
############



############

geneclasses=list(Fu_185_Hom.genes,
                 Fu_558_Hom.genes,
                 SFARI_Hom.genes,
                 Fu_185_Syn.genes,
                 Fu_558_Syn.genes,
                 SFARI_Syn.genes)

names(geneclasses)=c("Fu_185_Hom",
                     "Fu_558_Hom",
                     "SFARI_Hom",
                     "Fu_185_Syn",
                     "Fu_558_Syn",
                     "SFARI_Syn")

geneclassnames=names(geneclasses)

geneclassnames_abbrev=geneclassnames

############
background=list(Fu_185_Hom.background,
                     Fu_558_Hom.background,
                     SFARI_Hom.background,
                     Fu_185_Syn.background,
                     Fu_558_Syn.background,
                     SFARI_Syn.background)

names(background)=c("Fu_185_Hom",
                     "Fu_558_Hom",
                     "SFARI_Hom",
                     "Fu_185_Syn",
                     "Fu_558_Syn",
                     "SFARI_Syn")

#######################

############
proteinclasses=list(Fu_185_Hom.proteins,
                     Fu_558_Hom.proteins,
                     SFARI_Hom.proteins,
                     Fu_185_Syn.proteins,
                     Fu_558_Syn.proteins,
                     SFARI_Syn.proteins)

names(proteinclasses)=c("Fu_185_Hom",
                         "Fu_558_Hom",
                         "SFARI_Hom",
                         "Fu_185_Syn",
                         "Fu_558_Syn",
                         "SFARI_Syn")

#######################
# Initialize matrices to store results

ORmat <- matrix(nrow = length(geneclasses), ncol = 1, dimnames = list(names(geneclasses), as.character("comparison")))
logPmat <- matrix(nrow = length(geneclasses), ncol = 1, dimnames = list(names(geneclasses), as.character("comparison")))
Pmat <- matrix(nrow = length(geneclasses), ncol = 1, dimnames = list(names(geneclasses), as.character("comparison")))
DEprotpercent <- matrix(nrow = length(geneclasses), ncol = 1, dimnames = list(names(geneclasses), as.character("comparison")))
DEgenepercent <- matrix(nrow = length(geneclasses), ncol = 1, dimnames = list(names(geneclasses), as.character("comparison")))
#######################################################################
for (i in seq_along(geneclasses)) {
  
  # Intersect with background list
  genes2 <- geneclasses[[i]]
  genes2 <- data.frame(genes2)
  
  backgroundTotal=background[[i]]
  
  overlap_res <- genelistOverlap(data.frame(proteinclasses[[i]]),
                                 genes2,
                                 backgroundTotal,
                                 print_result = FALSE,
                                 header = FALSE)
  ORmat[i] <- overlap_res[[1]]$OR
  logPmat[i] <- -log10(overlap_res[[1]]$hypergeo_p)
  Pmat[i] <- overlap_res[[1]]$hypergeo_p
  
  #DEprotpercent[i] <- overlap_res[[1]]$percent_list1_background
  #DEgenepercent[i] <- overlap_res[[1]]$percent_list2_background
  
  DEprotpercent[i] <- overlap_res[[1]]$percent_overlap_list1
  DEgenepercent[i] <- overlap_res[[1]]$percent_overlap_list2
  
  # Write overlapping genes to a file
  write(overlap_res[[1]]$overlapping_genes,
        file = file.path(result_dir, sprintf("overlap_%s.txt", names(geneclasses)[i])))

}

# Store workdata (end of each run)
save_path <- file.path(result_dir, "overlap_results.RData")
save(ORmat, logPmat, Pmat, file = save_path)

#######################################################################
# load data
load("ASD_enrichment/Results/gene-enrichment/overlap_results.RData")

# create empty FDR matrix
FDRmat <- matrix(nrow = length(rownames(Pmat)), ncol = length(colnames(Pmat)), dimnames = list(rownames(Pmat), colnames(Pmat)))

# compute FDR
for (i in 1:dim(Pmat)[2]){
  FDRmat[,i] = p.adjust(Pmat[,i], method = "fdr")
} # for (i in 1:dim(Pmat)[2]){

kable(FDRmat)

# define order of comparisons
order=c("Fu_185_Hom",
        "Fu_185_Syn",
        "Fu_558_Hom",
        "Fu_558_Syn",
        "SFARI_Hom",
        "SFARI_Syn")

# reorder
logPmat=as.matrix(logPmat[order, ])
colnames(logPmat)="comparison"
ORmat=as.matrix(ORmat[order, ])
colnames(ORmat)="comparison"


# make figure
pheatmap(logPmat, display_numbers = round(ORmat,digits=2), 
         number_color = "black", fontsize_number = 12,
         show_rownames=TRUE,
         labels_col = c("Comparison"),
         color = colorRampPalette(c('white','red'))(100),
         cluster_rows = FALSE, cluster_cols = FALSE,
         breaks= seq(1.3,6, length=100))


DEprotpercent=as.matrix(DEprotpercent[order, ])
colnames(DEprotpercent)="DE protein %"

# make figure
pheatmap(DEprotpercent, display_numbers = round(DEprotpercent,digits=3), 
         number_color = "black", fontsize_number = 12,
         show_rownames=TRUE,
         labels_col = c("DE Protein % overlap"),
         color = colorRampPalette(c('yellow','purple'))(100),
         cluster_rows = FALSE, cluster_cols = FALSE,
         breaks= seq(0,0.10, length=100))



DEgenepercent=as.matrix(DEgenepercent[order, ])
colnames(DEgenepercent)="DE gene %"

# make figure
pheatmap(DEgenepercent, display_numbers = round(DEgenepercent,digits=3), 
         number_color = "black", fontsize_number = 12,
         show_rownames=TRUE,
         labels_col = c("DE Gene % overlap"),
         color = colorRampPalette(c('yellow','purple'))(100),
         cluster_rows = FALSE, cluster_cols = FALSE,
         breaks= seq(0,0.10, length=100))


#######################################################################
######################################################################
######
######
######
library(purrr)

# combine stats outputs
# convert list of lists to list of dataframes
stats=list(all_ORmat = ORmat,
           all_Pmat = Pmat,
           all_logPmat = logPmat,
           all_FDRmat = FDRmat)




wb <- createWorkbook()

for (name in names(stats)) {
  addWorksheet(wb, sheetName = name)
  writeData(wb, sheet = name, stats[[name]], rowNames = TRUE)
}


#saveWorkbook(wb, "ASD_enrichment/Results/overlap/gene-enrichment-stats.xlsx", overwrite = TRUE)


######
######
######
######################################################################
####
#########################################################
# make venn diagram

library(venneuler)

# Fu 185 ASD genes overlap with DE proteins
Fu.genes=c(Fu_185_Hom.genes,Fu_185_Syn.genes) # 160
Fu.genes=unique(Fu.genes) # 85

A=length(Fu.genes) # 85

B=length(Fu_185_Hom.proteins) # 216
C=length(Fu_185_Syn.proteins) # 284

A_B=length(intersect(Fu.genes,Fu_185_Hom.proteins)) # 1
A_C=length(intersect(Fu.genes, Fu_185_Syn.proteins)) # 3
B_C=length(intersect(Fu_185_Hom.proteins, Fu_185_Syn.proteins)) # 29

MyVenn=venneuler(c(
  A=A,
  B=B,
  C=C,
  "A&B"=A_B,
  "A&C"=A_C,
  "B&C"=B_C
))

MyVenn$labels <- c("Fu","Hom","Syn")
plot(MyVenn)

###################

# Fu 185 ASD & 377 NDD genes overlap with DE proteins
Fu.genes=c(Fu_558_Hom.genes,Fu_558_Syn.genes) # 383
Fu.genes=unique(Fu.genes) # 206

A=length(Fu.genes) # 206

B=length(Fu_558_Hom.proteins) # 216
C=length(Fu_558_Syn.proteins) # 284

A_B=length(intersect(Fu.genes,Fu_558_Hom.proteins)) # 5
A_C=length(intersect(Fu.genes, Fu_558_Syn.proteins)) # 8
B_C=length(intersect(Fu_185_Hom.proteins, Fu_558_Syn.proteins)) # 29

MyVenn=venneuler(c(
  A=A,
  B=B,
  C=C,
  "A&B"=A_B,
  "A&C"=A_C,
  "B&C"=B_C
))

MyVenn$labels <- c("Fu","Hom","Syn")
plot(MyVenn)

###################

###################

# SFARI (S,1,2) genes overlap with DE proteins
Fu.genes=c(SFARI_Hom.genes,SFARI_Syn.genes) # 853
Fu.genes=unique(Fu.genes) # 468

A=length(Fu.genes) # 468

B=length(SFARI_Hom.proteins) # 229
C=length(SFARI_Syn.proteins) # 303

A_B=length(intersect(Fu.genes,SFARI_Hom.proteins)) # 19
A_C=length(intersect(Fu.genes, SFARI_Syn.proteins)) # 25
B_C=length(intersect(SFARI_Hom.proteins, SFARI_Syn.proteins)) # 31

MyVenn=venneuler(c(
  A=A,
  B=B,
  C=C,
  "A&B"=A_B,
  "A&C"=A_C,
  "B&C"=B_C
))

MyVenn$labels <- c("SFARI","Hom","Syn")
plot(MyVenn)

###################
  
