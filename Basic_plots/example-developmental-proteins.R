# plot example developmental proteins

rm(list=ls()); gc()
options(stringsAsFactors = F)

# this portion performs a tracjectory analysis for the age effects

require(data.table)           # fast version of read.table
require(readxl)               # used for reading Excel format files
require(viridis)              # colorblind friendly color palette
library(xlsx)
library(ggplot2)
library(ggpubr)
library(tidyverse)

#tissue="HOMOGENATE-PROTEIN"
#tissue="SYNAPTOSOME-PROTEIN"

##########################################################################
setwd("C:/Users/Shelby_macdonaldlab/OneDrive - University of Pittsburgh/MacDonaldLab_Data/ASD_Development_March2023/Paper_files")
###################
# (1) load regression results for Homogenate data
# (2) load residual data for the effects of age only for Homogenate data
load(paste0("Analysis_Code/Code_Homogenate/WorkData/RESULTS-HOMOGENATE.RData"))
resid_Hom <- as.data.frame(read_excel("Supplemental_tables/Results/Residuals/Residual_results.xlsx", sheet="R.H.2"))
rownames(resid_Hom)=resid_Hom$...1
resid_Hom=resid_Hom[,-1]

# (3) load residual data for the effects of age and DX:age int added back only for Synaptosome data
resid_Syn <- as.data.frame(read_excel("Supplemental_tables/Results/Residuals/Residual_results.xlsx", sheet="R.S.5"))
rownames(resid_Syn)=resid_Syn$...1
resid_Syn=resid_Syn[,-1]

##########################################################################

#plot example proteins by age
#chose example proteins that were most regulated by age (increasing or decreasing) or no relationship change across age
#CFL2 (Q9Y281) - non-transitional in brain var

#PPT1 (P50897)
#NDUFA8 (P51970)
#EphA5 (P54756)

#select dataset to plot
#RESage=resid_Hom
RESage=resid_Syn

RESage.select=RESage[c("P50897","P51970","P54756"),]
RESage.select=rownames_to_column(RESage.select)
names(RESage.select)[1]="protein"
RESage.select$protein=as.factor(RESage.select$protein)

select=gather(RESage.select,ID,Value, S1866:S1672,factor_key = T)
select=left_join(select,meta.samples[,c(4,7,6)],by=c("ID"="ID"))
##########################################################################
#remove ASD subjects
select=filter(select, DX == "CTL")

# Plotting multiple Regression Lines
ggplot(select,aes(x=AGEyr,y=Value,color=protein, fill=protein))+
  geom_point(size=3,stroke=0.1,shape=21)+
  theme_classic()+
  geom_smooth(data=select, method=lm,se=T,fullrange=TRUE,
              aes(color=protein))+
  
  scale_fill_manual(values=c("#27AD81FF","#472D7BFF","#2C728EFF"))+
  scale_color_manual(values=c("#27AD81FF","#472D7BFF","#2C728EFF"))+
  
  ylim(-2,2)+
  theme(legend.position = c(.15, .75))+
  scale_x_continuous(breaks=seq(5,30,by=5))
