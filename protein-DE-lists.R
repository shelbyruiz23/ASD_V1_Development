rm(list=ls()); gc()
options(stringsAsFactors = F)

library(tidyverse)
library(readxl)

setwd("C:/Users/Shelby_macdonaldlab/OneDrive - University of Pittsburgh/MacDonaldLab_Data/ASD_Development_March2023/Paper_files")

##################################################################
####################
### Overview of Steps
# load DE results for Hom (DX), Syn (DX), and Syn (DX:AGE)
# remove proteins that map to none or multiple ensembl ids
# count the total number of ensembl genes to consider in Hom and Syn
# extract protein ensembl ids for DE protein lists
####################
##################################################################

##################################################################
####################
# Hom (DX)
# load homogenate results by diagnosis
load("Analysis_Code/Code_Homogenate/WorkData/RESULTS-HOMOGENATE.RData")
Hom=RE.PROTEIN$MLM.DX
# remove proteins that map to multiple ensemblids
Hom=filter(Hom, ensg != "" & !grepl("/",ensg)) # 4700
#check for proteins that map to the same gene
length(unique(Hom$ensg)) #4699
length(unique(Hom$protein)) #4700
# only retain duplicated gene - protein pair with the lowest q value
Hom <- Hom %>%
  group_by(ensg) %>%
  filter(q == min(q)) %>%
  ungroup()

Hom.map=Hom[,c(1,8,9)] # 4699

####################
# Syn (DX)
# load synaptosome results by diagnosis
load("Analysis_Code/Code_Synaptosome/WorkData/RESULTS-SYNAPTOSOME-updateINT.RData")
Syn=RE.PROTEIN$MLM.DX
# remove proteins that map to multiple ensemblids
Syn=filter(Syn, ensg != "" & !grepl("/",ensg)) # 4259
#check for proteins that map to the same gene
length(unique(Syn$ensg)) #4259
length(unique(Syn$protein)) #4259
Syn.map=Syn[,c(1,8,9)]

####################
# Syn (interaction of DX and age)
# load synaptosome results by diagnosis
load("Analysis_Code/Code_Synaptosome/WorkData/RESULTS-SYNAPTOSOME-updateINT.RData")
Syn2=RE.PROTEIN$MLM.DX.AGE
# remove proteins that map to multiple ensemblids
Syn2=filter(Syn2, ensg != "" & !grepl("/",ensg)) # 4259
##################################################################

######################################
# create list of all Hom and Syn proteins to overlap with enrichment lists
# union of Hom or Syn and enrichment list will be the background for all analyses
#Hom.list=as.character(Hom.map$ensg) # 4699
#Syn.list=as.character(Syn.map$ensg) # 4259

######################################
##################################################################

##################################################################
####################
#make list of DE proteins across Hom and Syn
# (1) Hom DE DX (q<0.1)
# (2) Hom DE DX (up) (q<0.1)
# (3) Hom DE DX (down) (q<0.1)
# (4) Syn DE DX (q<0.1)
# (5) Syn DE DX (up) (q<0.1)
# (6) Syn DE DX (down) (q<0.1)
# (7) Syn DE DX:AGE (q<0.1)
# (8) Syn DE DX:AGE (up) (q<0.1)
# (9) Syn DE DX:AGE (down) (q<0.1)

##################################################################


###############################
# (1) Hom DE DX (q<0.1)
Hom.genes=filter(Hom, q < 0.1) # 229
#Hom.genes=as.character(Hom.genes$ensg)

###############################

###############################
# (2) Hom DE DX (up) (q<0.1)
Hom.up.genes=filter(Hom, q < 0.1 & Estimate > 0) # 139
#Hom.up.genes=as.character(Hom.up.genes$ensg) 

###############################

###############################
# (3) Hom DE DX (down) (q<0.1)
Hom.down.genes=filter(Hom, q < 0.1 & Estimate < 0) # 90
#Hom.down.genes=as.character(Hom.down.genes$ensg) 

###############################

###############################
# (4) Syn DE DX (q<0.1)
Syn.genes=filter(Syn, q < 0.1) # 303
#Syn.genes=as.character(Syn.genes$ensg) 

###############################

###############################
# (5) Syn DE DX (up) (q<0.1)
Syn.up.genes=filter(Syn, q < 0.1 & Estimate > 0)# 172
#Syn.up.genes=as.character(Syn.up.genes$ensg) 

###############################

###############################
# (6) Syn DE DX (down) (q<0.1)
Syn.down.genes=filter(Syn, q < 0.1 & Estimate < 0)# 131
#Syn.down.genes=as.character(Syn.down.genes$ensg)

###############################

###############################
# (7) Syn DE DX:AGE (q<0.1)
Syn.int.genes=filter(Syn2, q < 0.1) # 304
#Syn.int.genes=as.character(Syn.int.genes$ensg)

###############################

###############################
# (8) Syn DE DX:AGE (up) (q<0.1)
Syn.int.up.genes=filter(Syn2, q < 0.1 & Estimate > 0) # 142
#Syn.int.up.genes=as.character(Syn.int.up.genes$ensg) 

###############################

###############################
# (9) Syn DE DX:AGE (down) (q<0.1)
Syn.int.down.genes=filter(Syn2, q < 0.1 & Estimate < 0)# 162
#Syn.int.down.genes=as.character(Syn.int.down.genes$ensg)

###############################


###############################
save(Hom.genes,
     Hom.up.genes,
     Hom.down.genes,
     Syn.genes,
     Syn.up.genes,
     Syn.down.genes,
     Syn.int.genes,
     Syn.int.up.genes,
     Syn.int.down.genes,
     #Hom.list,
     #Syn.list,
     Hom.map,
     Syn.map,
     file="Analysis_Code/ASD_enrichment/WorkData/protein-input.RData")


