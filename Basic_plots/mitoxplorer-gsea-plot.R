# plot mitoXplorer GSEA enrichment results

rm(list=ls()); gc()
options(stringsAsFactors = F)



require(data.table)           # fast version of read.table
require(readxl)               # used for reading Excel format files
require(viridis)              # colorblind friendly color palette
library(xlsx)
library(ggplot2)
library(ggpubr)
library(tidyverse)

##########################################################################

###################
# (1) load DE by DX results for Homogenate data
#load(paste0("../WorkData/RESULTS-HOMOGENATE.RData"))
Hom.DX=RE.PROTEIN$MLM.DX

# (2) load mitoxplorer gsea enrichment results for Homogenate data
Hom.mit.DX=as.data.frame(read_excel("Supplemental_tables/MitoXplorer_GSEA/MitoXplorer_GSEA.xlsx",sheet="Homogenate_DX"))

# (3/4) load DE by DX results for Synaptosome data
#load(paste0("../WorkData/RESULTS-SYNAPTOSOME-updateINT.RData"))
Syn.DX=RE.PROTEIN$MLM.DX
#Syn.DX.AGE=RE.PROTEIN$MLM.DX.AGE
Syn.DX.AGE=as.data.frame(read_excel("Supplemental_tables/Results/Protein/Results_protein.xlsx",sheet="Synaptosome_DEwithclustering"))

# (5/6) load mitoxplorer gsea enrichment results for Homogenate data
Syn.mit.DX=as.data.frame(read_excel("Supplemental_tables/MitoXplorer_GSEA/MitoXplorer_GSEA.xlsx",sheet="Synaptosome_DX"))
Syn.mit.DX.AGE=as.data.frame(read_excel("Supplemental_tables/MitoXplorer_GSEA/MitoXplorer_GSEA.xlsx",sheet="Synaptosome_DXandAGE"))

##########################################################################

##set up example dataframe
df.gsea=Syn.mit.DX.AGE
names(df.gsea)[11]="LEDGE_GENES"
df.DX=Syn.DX.AGE

names(df.DX)[24]="q_use"

# Create new columns in df.gsea
# (1) Number of upregulated proteins DE (q<0.1)
# (2) Number of downregulated proteins DE (q<0.1)
# (3) Number of proteins not DE (q>0.1)
# (4) Number of proteins in term

#########################
df.gsea$No_upregulated <- sapply(strsplit(df.gsea$LEDGE_GENES, ";"), function(genes) {
  sum(sapply(genes, function(gene) {
    match_result <- match(gene, df.DX$gene, nomatch = 0)
    if (match_result > 0 && df.DX$Estimate[match_result] > 0 && df.DX$q[match_result] < 0.1) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }))
})


df.gsea$No_downregulated <- sapply(strsplit(df.gsea$LEDGE_GENES, ";"), function(genes) {
  sum(sapply(genes, function(gene) {
    match_result <- match(gene, df.DX$gene, nomatch = 0)
    if (match_result > 0 && df.DX$Estimate[match_result] < 0 && df.DX$q[match_result] < 0.1) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }))
})

#########################
# just for Syn DX:AGE interaction

df.gsea$No_upregulated <- sapply(strsplit(df.gsea$LEDGE_GENES, ";"), function(genes) {
  sum(sapply(genes, function(gene) {
    match_result <- match(gene, df.DX$gene, nomatch = 0)
    if (match_result > 0 && df.DX$cluster.ASD[match_result] == "increasing" && df.DX$q_use[match_result] < 0.1) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }))
})


df.gsea$No_downregulated <- sapply(strsplit(df.gsea$LEDGE_GENES, ";"), function(genes) {
  sum(sapply(genes, function(gene) {
    match_result <- match(gene, df.DX$gene, nomatch = 0)
    if (match_result > 0 && df.DX$cluster.ASD[match_result] == "decreasing" && df.DX$q_use[match_result] < 0.1) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }))
})

df.gsea$No_stableregulated <- sapply(strsplit(df.gsea$LEDGE_GENES, ";"), function(genes) {
  sum(sapply(genes, function(gene) {
    match_result <- match(gene, df.DX$gene, nomatch = 0)
    if (match_result > 0 && df.DX$cluster.ASD[match_result] == "stable" && df.DX$q_use[match_result] < 0.1) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }))
})






#########################

df.gsea$No_null <- sapply(strsplit(df.gsea$LEDGE_GENES, ";"), function(genes) {
  sum(sapply(genes, function(gene) {
    match_result <- match(gene, df.DX$gene, nomatch = 0)
    if (match_result > 0 && df.DX$q_use[match_result] > 0.1) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }))
})

df.gsea$No_genes <- sapply(strsplit(df.gsea$LEDGE_GENES, ";"), function(genes) {
  length(genes)
})

#########################
#calculate the percentage of the total proteins are up, down, or null
#this value will be multiplied to the plotted enrichment score for color visualization

df.gsea$No_upregulated.percent=df.gsea$No_upregulated/df.gsea$No_genes
df.gsea$No_downregulated.percent=df.gsea$No_downregulated/df.gsea$No_genes
df.gsea$No_stableregulated.percent=df.gsea$No_stableregulated/df.gsea$No_genes
df.gsea$No_null.percent=df.gsea$No_null/df.gsea$No_genes


#########################
#calculate values for barplot
# use combined enrichment score and multiply the value by percentage of total genes for making stacked barplot
df.gsea$CS_up=df.gsea$`COMBINED SCORE`*df.gsea$No_upregulated.percent
df.gsea$CS_null=df.gsea$`COMBINED SCORE`*df.gsea$No_null.percent
df.gsea$CS_down=df.gsea$`COMBINED SCORE`*df.gsea$No_downregulated.percent
df.gsea$CS_stable=df.gsea$`COMBINED SCORE`*df.gsea$No_stableregulated.percent

# use combined enrichment score and multiply the value by percentage of total genes for making stacked barplot
#df.gsea$logPVAL=-log10(df.gsea$PVAL)

#df.gsea$CS_up=df.gsea$logPVAL*df.gsea$No_upregulated.percent
#df.gsea$CS_null=df.gsea$logPVAL*df.gsea$No_null.percent
#df.gsea$CS_down=df.gsea$logPVAL*df.gsea$No_downregulated.percent

#########################
#save info as file
#Hom.mit.DX=df.gsea
#Syn.mit.DX=df.gsea
Syn.mit.DX.AGE=df.gsea




##########################################################################
#make barplots

#collapse data into long format where term, condition (null, up, down) and value is individual
Hom.tmp=Hom.mit.DX[,c(1,6,19,20,21)]
#Hom.tmp=Hom.mit.DX[,c(1,22,19,20,21)]
names(Hom.tmp)[2]="CS"
Hom.tmp$PROCESS=factor(Hom.tmp$PROCESS)
Hom.tmp.long=gather(Hom.tmp, group,value, CS_up:CS_down,factor_key = T)

Syn.tmp=Syn.mit.DX[,c(1,6,19,20,21)]
#Syn.tmp=Syn.mit.DX[,c(1,22,19,20,21)]
names(Syn.tmp)[2]="CS"
Syn.tmp$PROCESS=factor(Syn.tmp$PROCESS)
Syn.tmp.long=gather(Syn.tmp, group,value, CS_up:CS_down,factor_key = T)

#SynINT.tmp=Syn.mit.DX.AGE[,c(1,6,19,20,21)]
#SynINT.tmp=Syn.mit.DX.AGE[,c(1,22,19,20,21)]

SynINT.tmp=Syn.mit.DX.AGE[,c(1,6,21:24)]
names(SynINT.tmp)[2]="CS"
SynINT.tmp$PROCESS=factor(SynINT.tmp$PROCESS)
SynINT.tmp.long=gather(SynINT.tmp, group,value, CS_up:CS_stable,factor_key = T)






######################
#select dataset
#df=Hom.tmp.long 
#df=Syn.tmp.long 
df=SynINT.tmp.long 

df=arrange(df, -CS)
#df$group=factor(df$group, levels = c("CS_up","CS_down","CS_null"))
df$group=factor(df$group, levels = c("CS_up","CS_down","CS_stable","CS_null"))

######################
# make barplot
df %>%
  mutate(PROCESS=fct_reorder(PROCESS, CS)) %>%
  ggplot( aes(x=PROCESS,y=value, fill=group))+
  geom_bar(position="stack",stat="identity")+
  scale_fill_manual(values=c("red","blue","grey"))+
  coord_flip()+
  theme_classic()+
  theme(legend.position = "none")+
  ylim(0,15.5)

#save as pdf
#portrat, 12 x 8

###save workData
#save(Hom.mit.DX,Syn.mit.DX,Syn.mit.DX.AGE,file="Analysis_Code/MitoXplorer/WorkData/MitoXplorer-gsea-results-q.RData")

######################
# make barplot
df %>%
  mutate(PROCESS=fct_reorder(PROCESS, CS)) %>%
  ggplot( aes(x=PROCESS,y=value, fill=group))+
  geom_bar(position="stack",stat="identity")+
  scale_fill_manual(values=c("#27AD81FF","#2C728EFF","#472D7BFF","grey"))+
  coord_flip()+
  theme_classic()+
  theme(legend.position = "none")+
  ylim(0,15.5)

