# (1) Bar plot of DE proteins for effect of DX and AGE and DX:AGE interaction
# (2) p value distributions for selected effects

#################
rm(list=ls()); gc()
options(stringsAsFactors = F)

#setwd
#setwd("C:/Users/Shelby_macdonaldlab/OneDrive - University of Pittsburgh/MacDonaldLab_Data/ASD_Development_March2023/Paper_files/Figure_code")

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

###################
#load regression results for Homogenate data
#load(paste0("../WorkData/LMER-REGRESSION-HOMOGENATE-PROTEIN.RData"))

#name Homogenate files
HOM.MLM.DX=reDX
HOM.MLM.AGE=reAGE

#load regression results for Synaptosome data
#load(paste0("../WorkData/LMER-REGRESSION-SYNAPTOSOME-PROTEIN.RData"))

#name Synaptosome files
SYN.MLM.DX=reDX
SYN.MLM.AGE=reAGE.noint
SYN.MLM.AGEandDX=reDX.AGE

###################
#Hom.DX
#q < 0.05 (32)
#q < 0.1  (230)
df=filter(HOM.MLM.DX, q < 0.05)
df=filter(HOM.MLM.DX, q < 0.1)

#Syn.DX
#q < 0.05 (93)
#q < 0.1  (304)
df=filter(SYN.MLM.DX, q < 0.05)
df=filter(SYN.MLM.DX, q < 0.1)

#Hom.AGE
#q < 0.05 (1261)
#q < 0.1  (1672)
df=filter(HOM.MLM.AGE, q < 0.05)
df=filter(HOM.MLM.AGE, q < 0.1)

#Syn.AGE
#q < 0.05 (912)
#q < 0.1  (1291)
df=filter(SYN.MLM.AGE, q < 0.05)
df=filter(SYN.MLM.AGE, q < 0.1)

#Syn.DX.AGE
#q < 0.05 (0)
#q < 0.1  (305)
df=filter(SYN.MLM.AGEandDX, q < 0.05)
df=filter(SYN.MLM.AGEandDX, q < 0.1)
###################

###################
#Bar plots
HOM.X=c("unique proteins",
        "proteins q < 0.1 | AGE",
        "proteins q < 0.05 | AGE",
        "proteins q < 0.1 | DX",
        "proteins q < 0.05 | DX")

HOM.Y=c(4727,1672,1261,230,32)
SYN.Y=c(4287,1291,912,304,93)

tmp=data.frame(HOM.X,HOM.Y,SYN.Y)

tmp.long=gather(tmp,label,number,HOM.Y:SYN.Y)

tmp.long$label <- factor(tmp.long$label, levels=c("HOM.Y", "SYN.Y"))

bp <- ggbarplot(
  tmp.long, x = "HOM.X", y = "number", 
  fill= "label",color="black", palette = c("#fee6ce","#fdae6b"),
  stat="identity",
  width=0.9,
  position = position_dodge(0.9),
  orientation="vertical",
  order = c("unique proteins",
            "proteins q < 0.1 | AGE",
            "proteins q < 0.05 | AGE",
            "proteins q < 0.1 | DX",
            "proteins q < 0.05 | DX"
  )
) +
  ylab("No. Proteins") +
  xlab("") +
  geom_text(aes(label=number), vjust=-0.2)+
  theme_light()+
  theme(axis.text.x = element_text(angle =45, hjust=1),legend.position=c(.6,0.8))



plot(bp)

#####frequency plots of p value distribution


#setup dataframe
HOM.DX.P=as.data.frame(HOM.MLM.DX$p)
HOM.AGE.P=as.data.frame(HOM.MLM.AGE$p)

SYN.DX.P=as.data.frame(SYN.MLM.DX$p)
SYN.AGE.P=as.data.frame(SYN.MLM.AGE$p)

#SYN.AGEandDX.P=as.data.frame(SYN.MLM.AGEandDX$p)



#plot
hist(HOM.DX.P$`HOM.MLM.DX$p`,breaks=seq(0,1,.01),las=1,xlab="empirical P",main="",col="#fee6ce",ylim=c(0,900))
hist(HOM.AGE.P$`HOM.MLM.AGE$p`,breaks=seq(0,1,.01),las=1,xlab="empirical P",main="",col="#fee6ce",ylim=c(0,900))

hist(SYN.DX.P$`SYN.MLM.DX$p`,breaks=seq(0,1,.01),las=1,xlab="empirical P",main="",col="#fdae6b",ylim=c(0,900))
hist(SYN.AGE.P$`SYN.MLM.AGE$p`,breaks=seq(0,1,.01),las=1,xlab="empirical P",main="",col="#fdae6b",ylim=c(0,900))


hist(SYN.AGEandDX.P$`SYN.MLM.AGEandDX$p`,breaks=seq(0,1,.01),las=1,xlab="empirical P",main="",col="#fdae6b",ylim=c(0,900))




######################################
