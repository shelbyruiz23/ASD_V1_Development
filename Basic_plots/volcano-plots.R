#Volcano plots- ASD Development 

#################


rm(list=ls()); gc()
options(stringsAsFactors = F)

#setwd
#setwd("C:/Users/Shelby_macdonaldlab/OneDrive - University of Pittsburgh/MacDonaldLab_Data/ASD_Development_March2023")

#load libraries
library(data.table)           # fast version of read.table
library(readxl)               # used for reading Excel format files
library(viridis)              # colorblind friendly color palette
library(matrixStats)          # to get statistics on matrices
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(ggrepel)
###################

# Volcano plot list
# (1) Homogenate - DX main fig
# (2) Synaptosome - DX main fig
# (3) Synaptosome - DX:AGE supp fig
# (4) Homogenate Unimputed data - DX supp fig
# (5) Synaptosome Unimputed data - DX supp fig


###################
# (1) load regression results for Homogenate data
#load(paste0("../WorkData/RESULTS-HOMOGENATE.RData"))

#name Homogenate files
HOM.MLM.DX=RE.PROTEIN$MLM.DX

# (2/3) load regression results for Synaptosome data
#load(paste0("../WorkData/RESULTS-SYNAPTOSOME-updateINT.RData"))

#name Synaptosome files
SYN.MLM.DX=RE.PROTEIN$MLM.DX
SYN.MLM.DX.AGE=RE.PROTEIN$MLM.DX.AGE

# (4) load regression results for Unimputed Homogenate data
#load(paste0("../WorkData/LMER-REGRESSION-HOMOGENATE-PROTEIN.RData"))

#name Homogenate files
U.HOM.MLM.DX=RE.DX

# (5) load regression results for unimputed Synaptosome data
#load(paste0("../WorkData/LMER-REGRESSION-SYNAPTOSOME-PROTEIN.RData"))

#name Synaptosome files
U.SYN.MLM.DX=RE.DX




#####################################################
# (1) select dataset to plot 
###############
#df=HOM.MLM.DX 
# 4727 proteins
# q<0.1 up (139), down (91)

#df=SYN.MLM.DX
# 4287 proteins
# q<0.1 up (173), down (131)

#df=SYN.MLM.DX.AGE
#changed xlim (-0.1,0.1)
# 4287 proteins
# q<0.1 up (142), down (163)

#df=U.HOM.MLM.DX
# 2578 proteins
# q<0.1 up (22), down (11)

df=U.SYN.MLM.DX
# 2455 proteins
# q<0.1 up (105), down (94)



t=filter(df, q < 0.1 & Estimate > 0) #up
t=filter(df, q < 0.1 & Estimate < 0) #down
###############

#(2) add labeling variables
df$Group=ifelse(df$q < 0.05, "q<0.05",
                    ifelse(df$q < 0.1 & df$q > 0.05, "q<0.1",
                           "nonsignificant"))


#(3) volcano plot

#define colors of points by group
cols=c("q<0.05"="#1037D1","q<0.1"="#7ad151","nonsignificant"="grey")

ggplot(df, 
       aes(x=Estimate, y=-log10(p)),
       fill=Group)+
  geom_point(shape=21,
             size=3,
             color="black",aes(fill=Group))+
  scale_fill_manual(values=cols)+
  
  xlim(-2,2)+
  ylim(0,6)+
  theme_light()+
  theme(legend.position = "none")+

geom_text_repel(aes(label=ifelse(q < 0.05 | q < 0.1 & Estimate > 0.3 | q < 0.1 & Estimate < -0.3, as.character(gene),'')),
                max.overlaps = 300)

#####################################################

#(3) volcano plot
#for Synaptosome interaction of DX:AGE
#effect size is much different than main effects

#define colors of points by group
cols=c("q<0.05"="#1037D1","q<0.1"="#7ad151","nonsignificant"="grey")

ggplot(df, 
       aes(x=Estimate, y=-log10(p)),
       fill=Group)+
  geom_point(shape=21,
             size=3,
             color="black",aes(fill=Group))+
  scale_fill_manual(values=cols)+
  
  xlim(-0.1,0.1)+
  ylim(0,6)+
  theme_light()+
  theme(legend.position = "none")+
  
  geom_text_repel(aes(label=ifelse(q < 0.1 & Estimate > 0.04 | q < 0.1 & Estimate < -0.02, as.character(gene),'')),
                  max.overlaps = 300)
