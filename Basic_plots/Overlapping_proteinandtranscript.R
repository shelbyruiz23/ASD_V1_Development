#Gandal, et al. (2022), overlap plots

rm(list=ls()); gc()
options(stringsAsFactors = F)

##########################################################################
#setwd
setwd("C:/Users/Shelby_macdonaldlab/OneDrive - University of Pittsburgh/MacDonaldLab_Data/ASD_Development_March2023")
##########################################################################

# Libraries
library(tidyverse)
library(viridis)
library(patchwork)
library(hrbrthemes)
library(ggplot2)
library(readxl)
library(ggrepel)


#####################################################################################
#####################################################################################
# comparison for V1 DE proteins and transcripts
#load homogenate MLM.DX data
RE.DX=read_excel("Bert_ASDpipeline_dontalter_2022_09/ASD-Development-Homogenate-Package-Sept2022/Results/HOMOGENATE-PROTEIN.xlsx",sheet="MLM.DX")
DX.q=filter(RE.DX,q < 0.1)

###########################
######################
# filter for DE proteins
#Protein=DX.q[,c(9,1,2,7,8)] # 230
# or retain all proteins
Protein=RE.DX[,c(9,1,2,7,8)] #4,727
###########################
######################

#load Gandal V1 supplemental table
T.DE=read_excel("Gandal_overlap/input/Haney-DE-V1-BioRX2020.xlsx")

#overlap of IDs in general
names(T.DE)[1]="ensg"
test=merge(RE.DX,T.DE,by="ensg") #4,619

T.DE.q=filter(T.DE, ASD_BA17_FDR < 0.1) # 4,782

###########################
######################
# filter for DE transcripts
#Transcript=T.DE.q[,c(1,9,10,11)]
# or retain all transcripts
Transcript=T.DE[,c(1,9,10,11)] # 24,836
###########################
######################

#merge overlapping protein and transcript data
Overlap=merge(Protein,Transcript,by="ensg") # all (4,619)

#add group column based on q value
Overlap$group=ifelse(Overlap$q < 0.05, "q0.05",
                     ifelse(Overlap$q < 0.1 & Overlap$q > 0.05,"q0.1",
                            "nonsig"))

Overlap$group=factor(Overlap$group,levels=c("nonsig","q0.1","q0.05"))
Overlap=arrange(Overlap, group)

#define colors of points by group
cols=c("q0.05"="#1037D1","q0.1"="#7ad151","nonsig"="grey")

ggplot(Overlap,
       aes(x=Estimate,y=ASD_BA17_logFC,fill=group))+
  geom_point(size=6,stroke=0.1,shape=21)+
  geom_abline(slope=1,color="darkgrey")+
  theme_light()+
  
  xlim(-1.65,1.65)+
  ylim(-1.65,1.65)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  
  scale_fill_manual(values=cols)+
  ylab("Transcript log2FC")+
  xlab("Protein log2FC") +
  #geom_text_repel(aes(label=gene),max.overlaps = 100)+
  theme(legend.position = "none")+
  labs(title = "protein, transcript all overlap", 
       subtitle = "n=4619, r=0.30, p=3.2E-94")

#correlation
res2=cor.test(Overlap$Estimate, Overlap$ASD_BA17_logFC,
              method="pearson")

# just DE proteins/transcripts
print(res2$estimate) #0.8149507 
print(res2$p.value) # 4.409405e-22

# all overlapping proteins/transcripts
print(res2$estimate) #0.2962546  
print(res2$p.value) # 3.173195e-94

#####################################################################################
#load synaptosome clustering by age data
RE.DX=read_excel("Bert_ASDpipeline_dontalter_2022_09/ASD-Development-Synaptosome-Package-Sept2022/Results/SYNAPTOSOME-PROTEIN.xlsx",sheet="MLM.DX")
DX.q=filter(RE.DX,q < 0.1)
Protein=DX.q[,c(9,1,2,7,8)]


#merge overlapping protein and transcript data
Overlap=merge(Protein,Transcript,by="ensg")

#add group column based on q value
Overlap$group=ifelse(Overlap$q < 0.05, "q0.05",
                     "q0.1")

Overlap$group=factor(Overlap$group,levels=c("q0.1","q0.05"))
Overlap=arrange(Overlap, group)

#define colors of points by group
cols=c("q0.05"="#1037D1","q0.1"="#7ad151")

ggplot(Overlap,
       aes(x=Estimate,y=ASD_BA17_logFC,fill=group))+
  geom_point(size=6,stroke=0.1,shape=21)+
  geom_abline(slope=1,color="darkgrey")+
  theme_light()+
  
  xlim(-1.55,1.55)+
  ylim(-1.55,1.55)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  
  scale_fill_manual(values=cols)+
  ylab("Transcript log2FC")+
  xlab("Protein log2FC") +
  geom_text_repel(aes(label=gene),max.overlaps = 100)+
  theme(legend.position = "none")
  



#####################################################################################
#############################
###############
#look at correlation of effect size by DX for all overlapping proteins

#load homogenate clustering by age data
RE.DX=read_excel("Bert_ASDpipeline_dontalter_2022_09/ASD-Development-Homogenate-Package-Sept2022/Results/HOMOGENATE-PROTEIN.xlsx",sheet="MLM.DX")
#DX.q=filter(RE.DX,q < 0.1)
Protein=RE.DX[,c(9,1,2,7,8)]




#load Gandal V1 supplemental table
T.DE=read_excel("Gandal_overlap/input/Haney-DE-V1-BioRX2020.xlsx")

#overlap of IDs in general
names(T.DE)[1]="ensg"
test=merge(RE.DX,T.DE,by="ensg") #4,619

#T.DE.q=filter(T.DE, ASD_BA17_FDR < 0.1)
Transcript=T.DE[,c(1,9,10,11)]


#merge overlapping protein and transcript data
Overlap=merge(Protein,Transcript,by="ensg")



ggplot(Overlap,
       aes(x=Estimate,y=ASD_BA17_logFC))+
  geom_point(size=3,stroke=0.1,shape=21,fill="grey")+
  geom_abline(slope=1,color="darkgrey")+
  theme_light()+
  
  xlim(-1.55,1.55)+
  ylim(-1.55,1.55)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  

  ylab("Transcript log2FC")+
  xlab("Protein log2FC") +
  #geom_text_repel(aes(label=gene),max.overlaps = 100)+
  theme(legend.position = "none")


#####################################################################################
#load synaptosome clustering by age data
RE.DX=read_excel("Bert_ASDpipeline_dontalter_2022_09/ASD-Development-Synaptosome-Package-Sept2022/Results/SYNAPTOSOME-PROTEIN.xlsx",sheet="MLM.DX")
#DX.q=filter(RE.DX,q < 0.1)
Protein=RE.DX[,c(9,1,2,7,8)]


#merge overlapping protein and transcript data
Overlap=merge(Protein,Transcript,by="ensg")



ggplot(Overlap,
       aes(x=Estimate,y=ASD_BA17_logFC))+
  geom_point(size=3,stroke=0.1,shape=21,fill="grey")+
  geom_abline(slope=1,color="darkgrey")+
  theme_light()+
  
  xlim(-1.55,1.55)+
  ylim(-1.55,1.55)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  
  scale_fill_manual(values=cols)+
  ylab("Transcript log2FC")+
  xlab("Protein log2FC") +

  theme(legend.position = "none")




#####################################################################################