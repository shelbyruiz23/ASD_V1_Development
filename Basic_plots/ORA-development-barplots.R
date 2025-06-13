#ORA from webgestalt for developmental protein clusters
#start: 05/15/2023

rm(list=ls()); gc()
options(stringsAsFactors = F)

library(tidyverse)
library(ggplot2)
library(readxl)

#criteria
#everything has an FDR < 0.05 for enrichment against the total proteins in the synaptosome dataset
#first check GO:Biological Process noRedundant.
#If no terms have an FDR < 0.05, use reactome pathway.
#If no reactome terms, list all proteins in group

#setwd
setwd("C:/Users/Shelby_macdonaldlab/OneDrive - University of Pittsburgh/MacDonaldLab_Data/ASD_Development_March2023")

I=read_excel("Paper_files/Supplemental_tables/ORA/Syn_age_ORA.xlsx",sheet="increasing")
I.stable=read_excel("Paper_files/Supplemental_tables/ORA/Syn_age_ORA.xlsx",sheet="I.stable")
I.dec=read_excel("Paper_files/Supplemental_tables/ORA/Syn_age_ORA.xlsx",sheet="I.decreasing")
S.inc=read_excel("Paper_files/Supplemental_tables/ORA/Syn_age_ORA.xlsx",sheet="S.increasing")
S=read_excel("Paper_files/Supplemental_tables/ORA/Syn_age_ORA.xlsx",sheet="stable")
S.dec=read_excel("Paper_files/Supplemental_tables/ORA/Syn_age_ORA.xlsx",sheet="S.decreasing")
D.inc=read_excel("Paper_files/Supplemental_tables/ORA/Syn_age_ORA.xlsx",sheet="D.increasing")
D.stable=read_excel("Paper_files/Supplemental_tables/ORA/Syn_age_ORA.xlsx",sheet="D.stable")
D=read_excel("Paper_files/Supplemental_tables/ORA/Syn_age_ORA.xlsx",sheet="decreasing")


I=filter(I, FDR < 0.05)
I.stable=filter(I.stable, FDR < 0.05)
I.dec=filter(I.dec, FDR < 0.05)
S.inc=filter(S.inc, FDR < 0.05)
S=filter(S, FDR < 0.05)
S.dec=filter(S.dec, FDR < 0.05)
D.inc=filter(D.inc, FDR < 0.05)
D.stable=filter(D.stable, FDR < 0.05)
D=filter(D, FDR < 0.05)

#if more than 10 terms, reduce to 5 terms based on pValue
S.dec=filter(S.dec, pValue < 0.0002)
S.inc=filter(S.inc, pValue < 0.0000235)

#combine terms all into one df
df=rbind(I,I.stable,I.dec,S.inc,S.dec,D.stable)
df$p=-log10(df$pValue)

#plot
ggplot(df, aes(x=description, y=p))+
  geom_bar(stat="identity",fill='#D3D3D3',width=0.6)+
  coord_flip()+
  theme_classic()
