# (1) Plot scatterplot comparison of overlapping Homogenate and Synaptosome effects for DX
# (2) Plot scatterplot comparison of overlapping Homogenate and Synaptosome effects for AGE
# (3) Plot scatterplot comparison of overlapping Homogenate and Synaptosome effects for DX:AGE

#################
rm(list=ls()); gc()
options(stringsAsFactors = F)

#setwd
setwd("C:/Users/Shelby_macdonaldlab/OneDrive - University of Pittsburgh/MacDonaldLab_Data/ASD_Development_March2023/Paper_files/Analysis_Code")

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
load(paste0("Code_Homogenate/WorkData/RESULTS-HOMOGENATE.RData"))

#name Homogenate files
HOM.MLM.DX=RE.PROTEIN$MLM.DX
HOM.MLM.AGE=RE.PROTEIN$MLM.AGE


#load regression results for Homogenate data with interaction forced into the model
load(paste0("Code_Homogenate/WorkData/LMER-REGRESSION-DXAGEint-HOMOGENATE-PROTEIN.RData"))

#name Homogenate files
HOM.MLM.AGEandDX=reDX_AGE

#load regression results for Synaptosome data
load(paste0("Code_Synaptosome/WorkData/RESULTS-SYNAPTOSOME-updateINT.RData"))

#name Synaptosome files
SYN.MLM.DX=RE.PROTEIN$MLM.DX
SYN.MLM.AGE=RE.PROTEIN$MLM.AGE.noint
SYN.MLM.AGEandDX=RE.PROTEIN$MLM.DX.AGE

############################################################################
#Estimate differences by overlapping proteins by diagnosis
Overlapping=merge(HOM.MLM.DX, SYN.MLM.DX,by="protein") #3852
Overlapping.Estimates=Overlapping[,c(1,2,10,7,15,8)]
names(Overlapping.Estimates)[2]="HOM.Estimate"
names(Overlapping.Estimates)[3]="SYN.Estimate"
names(Overlapping.Estimates)[4]="HOM.q"
names(Overlapping.Estimates)[5]="SYN.q"
names(Overlapping.Estimates)[6]="gene"

df=Overlapping.Estimates

#add label for changes by DX occurring in the same directio
df$q.overlap=ifelse(df$HOM.q < 0.1 & df$SYN.q < 0.1, "Hom & Syn",
                    ifelse(df$HOM.q < 0.1 & df$SYN.q > 0.1, "Hom only",
                           ifelse(df$SYN.q < 0.1 & df$HOM.q > 0.1, "Syn only",
                                  "neither")))

#df$q.overlap=factor(df$q.overlap, levels=c("Hom & Syn","Hom only","Syn only","neither"))
df$q.overlap=factor(df$q.overlap, levels=c("neither","Syn only","Hom only","Hom & Syn"))
df=arrange(df, q.overlap)


ggplot(df, aes(x=HOM.Estimate,y=SYN.Estimate,fill=q.overlap))+
  geom_point(size=3,stroke=0.1,shape=21)+
  scale_fill_manual(values=c("grey","#50A8EB","#C909EB","#180EE6"))+
  xlim(-2,2)+
  ylim(-2,2)+
  geom_abline(slope=1)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept=0)+
  theme_light()+
  theme(legend.position = "none")



res=cor.test(df$HOM.Estimate, df$SYN.Estimate,
             method="spearman")

print(res$estimate) #0.2730203 
print(res$p.value) #8.128653e-67
############################################################################

######################################
#Estimate differences by overlapping proteins by AGE
Overlapping=merge(HOM.MLM.AGE, SYN.MLM.AGE,by="protein")
Overlapping.Estimates=Overlapping[,c(1,2,10,7,15,8)]
names(Overlapping.Estimates)[2]="HOM.Estimate"
names(Overlapping.Estimates)[3]="SYN.Estimate"
names(Overlapping.Estimates)[4]="HOM.q"
names(Overlapping.Estimates)[5]="SYN.q"
names(Overlapping.Estimates)[6]="gene"

tmp=Overlapping.Estimates

#tmp$HOM4=tmp$HOM.Estimate*4
#tmp$HOM32=tmp$HOM.Estimate*32
#tmp$HOM.Estimate.adj=tmp$HOM32 - tmp$HOM4

#tmp$SYN4=tmp$SYN.Estimate*4
#tmp$SYN32=tmp$SYN.Estimate*32
#tmp$SYN.Estimate.adj=tmp$SYN32 - tmp$SYN4


df=tmp

#add label for changes by DX occurring in the same directio
df$q.overlap=ifelse(df$HOM.q < 0.1 & df$SYN.q < 0.1, "Hom & Syn",
                    ifelse(df$HOM.q < 0.1 & df$SYN.q > 0.1, "Hom only",
                           ifelse(df$SYN.q < 0.1 & df$HOM.q > 0.1, "Syn only",
                                  "neither")))

df$q.overlap=factor(df$q.overlap, levels=c("neither","Syn only","Hom only","Hom & Syn"))
df=arrange(df, q.overlap)


ggplot(df, aes(x=HOM.Estimate,y=SYN.Estimate,fill=q.overlap))+
  geom_point(size=3,stroke=0.1,shape=21)+
  scale_fill_manual(values=c("grey","#50A8EB","#C909EB","#180EE6"))+
  xlim(-0.15,0.15)+
  ylim(-0.15,0.15)+
  geom_abline(slope=1)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept=0)+
  theme_light()+
  theme(legend.position = "none")

#correlation
res2=cor.test(df$HOM.Estimate, df$SYN.Estimate,
              method="spearman")

print(res2$estimate) #0.6597851 
print(res2$p.value) #0

############################################################################

############################################################################

######################################
# if you force DX:AGE into Homogenate model

#Estimate differences by overlapping proteins by DX:AGE

Overlapping=merge(HOM.MLM.AGEandDX, SYN.MLM.AGEandDX,by="protein")
Overlapping.Estimates=Overlapping[,c(1,2,8,7,13,14)]
names(Overlapping.Estimates)[2]="HOM.Estimate"
names(Overlapping.Estimates)[3]="SYN.Estimate"
names(Overlapping.Estimates)[4]="HOM.q"
names(Overlapping.Estimates)[5]="SYN.q"
names(Overlapping.Estimates)[6]="gene"

tmp=Overlapping.Estimates



df=tmp

#add label for changes by DX:AGE occurring in the same directio
df$q.overlap=ifelse(df$SYN.q < 0.1 & df$HOM.q > 0.1, "Syn only",
                                  "neither")

df$q.overlap=factor(df$q.overlap, levels=c("neither","Syn only"))
df=arrange(df, q.overlap)


ggplot(df, aes(x=HOM.Estimate,y=SYN.Estimate,fill=q.overlap))+
  geom_point(size=3,stroke=0.1,shape=21)+
  scale_fill_manual(values=c("grey","#50A8EB"))+
  xlim(-0.15,0.15)+
  ylim(-0.15,0.15)+
  geom_abline(slope=1)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept=0)+
  theme_light()+
  theme(legend.position = "none")

#correlation
res2=cor.test(df$HOM.Estimate, df$SYN.Estimate,
              method="spearman")

print(res2$estimate) #0.324
print(res2$p.value) #6.77x10-95

############################################################################


df.Hom=read_excel("C:/Users/Shelby_macdonaldlab/OneDrive - University of Pittsburgh/MacDonaldLab_Data/ASD_Development_March2023/Paper_files/Supplemental_tables/Age_flexmix/Age_flexmix_results.xlsx",sheet="Homogenate_flexmix")
df.Syn=read_excel("C:/Users/Shelby_macdonaldlab/OneDrive - University of Pittsburgh/MacDonaldLab_Data/ASD_Development_March2023/Paper_files/Supplemental_tables/Age_flexmix/Age_flexmix_results.xlsx",sheet="Synaptosome_flexmix")


df.Hom=df.Hom[,c(1,5,8)]
df.Syn=df.Syn[,c(1,5,8)]

names(df.Hom)[2]="b.asd_Hom"
names(df.Hom)[3]="b.ctl_Hom"

names(df.Syn)[2]="b.asd_Syn"
names(df.Syn)[3]="b.ctl_Syn"

Overlapping=merge(df.Hom, df.Syn,by="protein")

################
# add age effects here
HOM.MLM.AGE=HOM.MLM.AGE[,c(1,7)]
names(HOM.MLM.AGE)[2]="q.age_hom"
SYN.MLM.AGE=SYN.MLM.AGE[,c(1,7)]
names(SYN.MLM.AGE)[2]="q.age_syn"

Overlapping=left_join(Overlapping,HOM.MLM.AGE,by="protein")
Overlapping=left_join(Overlapping,SYN.MLM.AGE,by="protein")
##################

df=Overlapping

#add label for changes by DX occurring in the same directio
df$q.overlap=ifelse(df$q.age_hom < 0.1 & df$q.age_syn < 0.1, "Hom & Syn",
                    ifelse(df$q.age_hom < 0.1 & df$q.age_syn > 0.1, "Hom only",
                           ifelse(df$q.age_syn < 0.1 & df$q.age_hom > 0.1, "Syn only",
                                  "neither")))

df$q.overlap=factor(df$q.overlap, levels=c("neither","Syn only","Hom only","Hom & Syn"))
df=arrange(df, q.overlap)


ggplot(df, aes(x=b.ctl_Hom,y=b.ctl_Syn,fill=q.overlap))+
  geom_point(size=3,stroke=0.1,shape=21)+
  scale_fill_manual(values=c("grey","#50A8EB","#C909EB","#180EE6"))+
  xlim(-0.15,0.15)+
  ylim(-0.15,0.15)+
  geom_abline(slope=1)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept=0)+
  theme_light()+
  theme(legend.position = "none")

#correlation
res2=cor.test(df$b.ctl_Hom, df$b.ctl_Syn,
              method="spearman")

print(res2$estimate) #0.631
print(res2$p.value) #0


#####################################################
ggplot(df, aes(x=b.asd_Hom,y=b.asd_Syn,fill=q.overlap))+
  geom_point(size=3,stroke=0.1,shape=21)+
  scale_fill_manual(values=c("grey","#50A8EB","#C909EB","#180EE6"))+
  xlim(-0.15,0.15)+
  ylim(-0.15,0.15)+
  geom_abline(slope=1)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept=0)+
  theme_light()+
  theme(legend.position = "none")


res2=cor.test(df$b.asd_Hom, df$b.asd_Syn,
              method="spearman")

print(res2$estimate) #0.570
print(res2$p.value) #0
