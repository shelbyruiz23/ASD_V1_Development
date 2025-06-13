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
# load ASD risk gene data from Fu 185 (Fu 2022)
# supplementary Data 11

############
# (1) Fu et al, 185 genes (Hom)
Fu.list=read_excel("ASD_enrichment/input/Fu_2022_supplementary.xlsx",sheet="Supplementary Table 11")
names(Fu.list)[3]="ensg"

# Find the overlap with the background (Hom)
Fu.list=Fu.list[Fu.list$ensg %in% Hom.map$ensg, ] # 4475
Fu_185_Hom.background=length(Fu.list$ensg)
  
# From identified proteins, filter for genes that were DE for ASD in Fu
Fu_185=filter(Fu.list, ASD185 == "TRUE")
Fu_185_Hom.genes=as.character(Fu_185$gene) # 83

Fu_185_Hom.proteins=Hom.genes[Hom.genes$ensg %in% Fu.list$ensg, ] # 216
Fu_185_Hom.proteins=as.character(Fu_185_Hom.proteins$gene)

############

############
# (2) Fu et al, 185 genes (Syn)
Fu.list=read_excel("ASD_enrichment/input/Fu_2022_supplementary.xlsx",sheet="Supplementary Table 11")
names(Fu.list)[3]="ensg"

# Find the overlap with the background (Syn)
Fu.list=Fu.list[Fu.list$ensg %in% Syn.map$ensg, ] # 4055
Fu_185_Syn.background=length(Fu.list$ensg)

# From identified proteins, filter for genes that were DE for ASD in Fu
Fu_185=filter(Fu.list, ASD185 == "TRUE")
Fu_185_Syn.genes=as.character(Fu_185$gene) # 77

Fu_185_Syn.proteins=Syn.genes[Syn.genes$ensg %in% Fu.list$ensg, ] # 284
Fu_185_Syn.proteins=as.character(Fu_185_Syn.proteins$gene)

############

############
# (3) Fu et al, 185 genes (ASD) + 373 genes (NDD) (unique only) (Hom)
Fu.list=read_excel("ASD_enrichment/input/Fu_2022_supplementary.xlsx",sheet="Supplementary Table 11")
names(Fu.list)[3]="ensg"

# Find the overlap with the background (Hom)
Fu.list=Fu.list[Fu.list$ensg %in% Hom.map$ensg, ] # 4475
Fu_558_Hom.background=length(Fu.list$ensg)

# From identified proteins, filter for genes that were DE for ASD in Fu
Fu_558=filter(Fu.list, ASD185 == "TRUE" | NDD373 == TRUE)
Fu_558_Hom.genes=as.character(Fu_558$gene) # 200

Fu_558_Hom.proteins=Hom.genes[Hom.genes$ensg %in% Fu.list$ensg, ] # 216
Fu_558_Hom.proteins=as.character(Fu_558_Hom.proteins$gene)

############

############
# (4) Fu et al, 185 genes (ASD) + 373 genes (NDD) (unique only) (Syn)
Fu.list=read_excel("ASD_enrichment/input/Fu_2022_supplementary.xlsx",sheet="Supplementary Table 11")
names(Fu.list)[3]="ensg"

# Find the overlap with the background (Syn)
Fu.list=Fu.list[Fu.list$ensg %in% Syn.map$ensg, ] # 4055
Fu_558_Syn.background=length(Fu.list$ensg)

# From identified proteins, filter for genes that were DE for ASD in Fu
Fu_558=filter(Fu.list, ASD185 == "TRUE" | NDD373 == TRUE)
Fu_558_Syn.genes=as.character(Fu_558$gene) # 183

Fu_558_Syn.proteins=Syn.genes[Syn.genes$ensg %in% Fu.list$ensg, ] # 284
Fu_558_Syn.proteins=as.character(Fu_558_Syn.proteins$gene)

############


############
# (5) SFARI genes (S,1,2) (Hom)
# use all identified proteins as background

SFARI.list=read.csv("ASD_enrichment/input/SFARI-Gene_genes_11-03-2023release_12-05-2023export.csv")
names(SFARI.list)[4]="ensg"

# Find the overlap with the background (Hom)
SFARI_Hom.background=4699

# From identified proteins, filter for genes that have a gene score of S, 1, or 2
SFARI.list=SFARI.list[SFARI.list$ensg %in% Hom.map$ensg, ] # 476
SFARI=filter(SFARI.list, gene.score == 1 | gene.score == 2 | syndromic == 1)
SFARI_Hom.genes=as.character(SFARI$gene.symbol) # 438

SFARI_Hom.proteins=Hom.genes # 229
SFARI_Hom.proteins=as.character(SFARI_Hom.proteins$gene)


###############################

############
# (6) SFARI genes (S,1,2) (Syn)
# use all identified proteins as background

SFARI.list=read.csv("ASD_enrichment/input/SFARI-Gene_genes_11-03-2023release_12-05-2023export.csv")
names(SFARI.list)[4]="ensg"

# Find the overlap with the background (Syn)
SFARI_Syn.background=4259

# From identified proteins, filter for genes that have a gene score of S, 1, or 2
SFARI.list=SFARI.list[SFARI.list$ensg %in% Syn.map$ensg, ] # 453
SFARI=filter(SFARI.list, gene.score == 1 | gene.score == 2 | syndromic == 1) # 415
SFARI_Syn.genes=as.character(SFARI$gene.symbol) # 415

SFARI_Syn.proteins=Syn.genes # 303
SFARI_Syn.proteins=as.character(SFARI_Syn.proteins$gene)


###############################

# save workdata
save(Fu_185_Hom.background,
     Fu_185_Hom.genes,
     Fu_185_Hom.proteins,
     
     Fu_185_Syn.background,
     Fu_185_Syn.genes,
     Fu_185_Syn.proteins,
     
     Fu_558_Hom.background,
     Fu_558_Hom.genes,
     Fu_558_Hom.proteins,
     
     Fu_558_Syn.background,
     Fu_558_Syn.genes,
     Fu_558_Syn.proteins,
     
     SFARI_Hom.background,
     SFARI_Hom.genes,
     SFARI_Hom.proteins,
     
     SFARI_Syn.background,
     SFARI_Syn.genes,
     SFARI_Syn.proteins,
     
     
     
     file="ASD_enrichment/WorkData/gene-list.RData")
