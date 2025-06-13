#Plot GSEA plots
#test out subcellular locations
#2023_0403


###################################################


rm(list=ls()); gc()
options(stringsAsFactors = F)

library(tidyverse)
library(readxl)
library(data.table)
library(tidyr)

###################################################
#setwd("C:/Users/Shelby_macdonaldlab/OneDrive - University of Pittsburgh/MacDonaldLab_Data/ASD_Development_March2023")
###################################################
#### Part 0 - Load these
##load .RData
#load('GSEA/2023_0327/GSEA-enrichment.RData')

###################################################



###################################################
#keep only MF,BP, 

HOM.GSEA.DX.downreg_q0.1.short=filter(HOM.GSEA.DX.downreg_q0.1, database == "GO; cellular components" | 
database == "GO; molecular functions" | database == "GO; biological processes")

HOM.GSEA.DX.upreg_q0.1.short=filter(HOM.GSEA.DX.upreg_q0.1, database == "GO; cellular components" | 
                                        database == "GO; molecular functions" | database == "GO; biological processes")


SYN.GSEA.DX.downreg_q0.1.short=filter(SYN.GSEA.DX.downreg_q0.1, database == "GO; cellular components" | 
                                        database == "GO; molecular functions" | database == "GO; biological processes")

SYN.GSEA.DX.upreg_q0.1.short=filter(SYN.GSEA.DX.upreg_q0.1, database == "GO; cellular components" | 
                                      database == "GO; molecular functions" | database == "GO; biological processes")

SYN.GSEA.DXandAGE.downreg_q0.1.short=filter(SYN.GSEA.DXandAGE.downreg_q0.1, database == "GO; cellular components" | 
                                        database == "GO; molecular functions" | database == "GO; biological processes")

SYN.GSEA.DXandAGE.upreg_q0.1.short=filter(SYN.GSEA.DXandAGE.upreg_q0.1, database == "GO; cellular components" | 
                                      database == "GO; molecular functions" | database == "GO; biological processes")

###################################################
#select top 3 enrichment terms for BP, CC, MF for all groups based on adjP value
#if terms are completely redundant, choose the highest and move to the next most enriched term
#if there are not 3 terms for each enrichment term (BP,CC,MF), evenly disperse enrichments until you reach a total of 9 terms

####
#SYN Downreg
df=SYN.GSEA.DX.downreg_q0.1.short

df2=df %>%
  arrange(adjP) %>%
  group_by(database) %>%
  slice(1:3)

SYN.downreg=df2
###

####
#SYN Upreg
#redundant terms
df=SYN.GSEA.DX.upreg_q0.1.short

df2=df %>%
  arrange(adjP) %>%
  group_by(database) %>%
  slice(1:4)

SYN.upreg=df2
SYN.upreg=SYN.upreg[c(-2,-8,-11),]
###

####
#HOM Downreg
#redundant terms
df=HOM.GSEA.DX.downreg_q0.1.short

df2=df %>%
  arrange(adjP) %>%
  group_by(database) %>%
  slice(1:9)

HOM.downreg=df2
HOM.downreg=HOM.downreg[-9,]
###

####
#HOM Upreg

df=HOM.GSEA.DX.upreg_q0.1.short

df2=df %>%
  arrange(adjP) %>%
  group_by(database) %>%
  slice(1:6)

HOM.upreg=df2

###

####log2 convert the overrepresentation values
SYN.downreg$ORlog2=log2(SYN.downreg$OR)
SYN.upreg$ORlog2=log2(SYN.upreg$OR)

HOM.downreg$ORlog2=log2(HOM.downreg$OR)
HOM.upreg$ORlog2=log2(HOM.upreg$OR)


###insert value of 10 for infinity oddsratio in Syn upregulated
SYN.upreg$ORlog2[SYN.upreg$GSterm == "GOMF_AMINOPHOSPHOLIPID_FLIPPASE_ACTIVITY"] <- 10





#plotting

library(ggplot2)
library(forcats)



##theme function
theme_dose <- function(font.size=14) {
  theme_bw() +
    theme(axis.text.x = element_text(colour = "black",
                                     size = font.size, vjust =1 ),
          axis.text.y = element_text(colour = "black",
                                     size = font.size, hjust =1 ),
          axis.title = element_text(margin=margin(10, 5, 0, 0),
                                    color = "black",
                                    size = font.size),
          axis.title.y = element_text(angle=90)
    )
}




##############################################################
#bind all dataframes into one plot
df=rbind(HOM.downreg,HOM.upreg,SYN.downreg,SYN.upreg)


ggplot(df, showCategory = 10,
       
       aes(ORlog2, fct_reorder(GSterm, ORlog2))) +
  
  geom_segment(aes(xend=0, yend = GSterm)) +
  
  geom_point(aes(color=adjP, size = DEgenes)) +
  
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
                        
                        trans = "log10",
                        
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  
  scale_size_continuous(range=c(2, 10)) +
  
  theme_dose(6) +
  
  xlab("log2(Odds Ratio)") +
  xlim(0,10)+
  ylab(NULL) +
  
  
  
  ggtitle("Top 9 Enrichment terms for HOM and SYN")

##############################################################






