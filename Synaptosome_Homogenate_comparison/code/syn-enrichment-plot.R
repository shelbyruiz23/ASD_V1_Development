rm(list=ls()); gc()
options(stringsAsFactors = F)

tissue="HOM-SYN-COMP-POOLS"
########################################################
# performed syngo enrichment for proteins with a FC > 1
# syngo plot looks great

# overlap results with synaptic proteins

# annotate regression results with synaptic or nonsynaptic terms
# annotate regression results with pre or post or both synaptic terms

# make volcano plot annotated by synaptic or nonsynaptic terms
# make volcano plot annotated by pre, post, both terms
# make density plot annotated by synaptic or nonsynaptic terms
# make density plot annotated by pre, post, both terms
library(readxl)
library(ggplot2)
library(viridisLite)
library(scales)
library(tidyverse)
library(ggrepel)
library(limma)
########################################################

# load Regression results
load(paste0("Analysis_Code/Synaptosome_enrichment/WorkData/LMER-REGRESSION-",tissue,".RData"))

# export bulk syngo database release (version=20231201, export date=20231216)
syngo.onto=as.data.frame(read_excel("Analysis_Code/Synaptosome_enrichment/input/syngo_ontologies.xlsx"))
########################################################

########################################################
#################################
##########
# make list of syngo ontologies split by pre and postsynaptic

# filter for cellular location only
syngo.onto=filter(syngo.onto, domain == "CC"| domain == "BP" & shortname == "synaptic"| name == "process in the presynapse" & shortname == "presynaptic" | name == "process in the postsynapse" & shortname == "postsynaptic")

# filter for parent terms synapse, presynapse, postsynapse, synaptic
# other additional terms for later are parent terms like synaptic membrane
syngo.onto=filter(syngo.onto, shortname == "synapse" | shortname == "presynapse" | shortname == "postsynapse" | shortname == "synaptic" | shortname == "presynaptic" | shortname == "postsynaptic")

# filter for relevant columns
syngo.onto=syngo.onto[,c(5,9)]
names(syngo.onto)[2]="ensg"

# Separate the ensg column into separate rows
syngo.onto <- syngo.onto %>%
  separate_rows(ensg, sep = ",\\s*") %>%
  # Add a row number for each shortname group
  group_by(shortname) %>%
  mutate(row_id = row_number()) %>%
  # Pivot the dataframe to wide format
  pivot_wider(names_from = shortname, values_from = ensg) %>%
  # Remove the row_id column
  select(-row_id)

#concatenate BP and CC categories
pre=c(unique(syngo.onto$presynapse),unique(syngo.onto$presynaptic)) #1011
post=c(unique(syngo.onto$postsynapse),unique(syngo.onto$postsynaptic)) #1193
syn=c(unique(syngo.onto$synapse),unique(syngo.onto$synaptic)) #2614

syn.list=data.frame(ensg=unique(syn))

# Assign corresponding values based on conditions
syn.categories <- syn.list %>%
  mutate(
    category = case_when(
      ensg %in% pre & !(ensg %in% post) ~ "presynapse",
      ensg %in% post & !(ensg %in% pre) ~ "postsynapse",
      !(ensg %in% pre) & !(ensg %in% post) | ensg %in% pre & ensg %in% post ~ "synapse",
      TRUE ~ NA_character_
    )
  )

#summary_table <- as.data.frame(table(syn.categories$category))
# post 599
# pre 364
# syn 640


# add synaptic annotations to regression dataframe
RETissue=merge(RETissue,syn.categories,by="ensg",all.x = T)

#check uniprot for synaptic annotation of any of the 14 proteins that did not map to ensmbl genes
# Q96NW7 = postsynapse, LRRC7
# Q9BYB0 = postsynapse, SHANK3
# add these annotations to the df

#loop through new values
new_value="LRRC7"
row_index=which(RETissue$protein == "Q96NW7")

RETissue[row_index, "gene"]=new_value

#save version
RETissue.set=RETissue

#add labels
RETissue.set[["category"]]=ifelse(is.na(RETissue.set[["category"]]),"nonsynapse",RETissue.set[["category"]])
RETissue.set$group= ifelse(RETissue.set$category == "presynapse" | RETissue.set$category == "postsynapse" | RETissue.set$category == "synapse", "synaptic",
                          "nonsynaptic")
########################################################

###################################################
#plotting by basic synaptic or not
######
df=RETissue.set

df$group=factor(df$group,levels=c("nonsynaptic","synaptic"))
df=arrange(df,group)




ggplot(df,
       aes(x=Estimate,y=-log10(p),fill=group))+
  geom_point(size=3,stroke=0.1,shape=21,alpha=0.5)+
  xlim(-6,6)+
  geom_vline(xintercept = 0,color="darkgrey")+
  
  scale_fill_manual(values=c("grey","#440154FF"))+
  ylab("-log10(p)")+
  xlab("Synaptosome/Homogenate")+
  
  theme_light()

#save as 8x8 landscape

###################################################
#plotting with labeling genes
###

#for individual plots by subject distribution
# syp
# syn1
# syngap1
# kif5c
# icam5
#shank3 

# downregulated expected
#pvalb

df$gene.label=ifelse(df$gene == "GAP43"|
                       df$gene == "EIF4B"|
                       df$gene == "MARCKS"|
                       df$gene == "SCN2A"|
                       df$gene == "MTOR"|
                       df$gene == "SLC17A7"|
                       df$gene == "KALRN"|
                       #df$gene == "SCAMP5"|
                       df$gene == "EPHB2"|
                       df$gene == "SEPTIN5"|
                       df$gene == "GRM2"|
                       #df$gene == "FLOT2"|
                       #df$gene == "SNAP25"|
                                 df$gene == "NAPG"|
                                 df$gene == "NRCAM"|
                                 df$gene == "ATP2B2"|
                                 df$gene == "RAP2A"|
                                 df$gene == "SYN1"|
                                 df$gene == "MAP2"|
                                 df$gene == "SYP"|
                                 df$gene == "GRM5"|
                                 df$gene == "SV2A"|
                                 df$gene == "SYT3"|
                                 df$gene == "VAMP2"|
                                 #df$gene == "ACTN2"|
                                 #df$gene == "NCAM1"|
                                 df$gene == "SYNGAP1"|
                                 #df$gene == "GRIN1"|
                       #df$gene == "GRIN2A"|
                       df$gene == "SYP"|
                       df$gene == "SYN2"|
                       df$gene == "NRXN3"|
                       df$gene == "NLGN3"|
                       df$gene == "KIF5C"|
                       df$gene == "ICAM5"|
                       df$gene == "PVALB"|
                                 df$gene == "CAMK2A"|
                                 df$gene == "GRIA1"|
                                 df$gene == "DLG4"|
                                 #df$gene == "SHANK3"|
                                 #df$gene == "CNR1"|
                                 #df$gene == "HCN1"|
                                 #df$gene == "GRIN2B"|
                                 df$gene == "MTOR"|
                                 df$gene == "COL4A2"|
                                 df$gene == "DDX3Y"|
                                 #df$gene == "H1-2"|
                                 df$gene == "PVALB"|
                                 #df$gene == "SPTAN1"|
                                 #df$gene == "SLC1A2"|
                                 df$gene == "KCNC3"|
                                 #df$gene == "ARHGEF9"|
                                 #df$gene == "STXBP1"|
                                 df$gene == "GFAP"|
                                 #df$gene == "GABRA3"|
                       df$gene == "SNRPD1"|
                       df$gene == "TMEM63B"|
                       df$gene == "FAM234A"|
                       df$gene == "SLC22A23"|
                       df$gene == "CHGA"|
                       df$gene == "PTMA"|
                                 df$gene == "GPHN","YES",
                               "NO")



ggplot(df,
       aes(x=Estimate,y=-log10(p),fill=group))+
  geom_point(size=3,stroke=0.1,shape=21,alpha=0.6)+
  xlim(-6,6)+
  geom_vline(xintercept = 0,color="darkgrey")+
  
  scale_fill_manual(values=c("grey","#440154FF"))+
  ylab("-log10(p)")+
  xlab("Synaptosome/Homogenate")+
  
  theme_light() +
  
  geom_text_repel(aes(label=ifelse(gene.label == "YES", as.character(gene),"")),
                  max.overlaps=2000,
                  min.segment.length = 0
  )


###################################################
#plotting by more specific synaptic location
#######
df$category=factor(df$category,levels=c("nonsynapse","synapse","presynapse","postsynapse"))
df=arrange(df,category)

ggplot(df,
       aes(x=Estimate,y=-log10(p),fill=category))+
  geom_point(size=3,stroke=0.1,shape=21,alpha=0.6)+
  xlim(-6,6)+
  geom_vline(xintercept = 0,color="darkgrey")+
  
  scale_fill_manual(values=c("grey","#440154FF","#29AF7FFF","#FDE725FF"))+
  ylab("-log10(p)")+
  xlab("Synaptosome/Homogenate")+
  
  theme_light()

#######
df2=filter(df, category != "nonsynapse")
df2$category=factor(df2$category,levels=c("synapse","presynapse","postsynapse"))
df2=arrange(df2,category)

ggplot(df2,
       aes(x=Estimate,y=-log10(p),fill=category))+
  geom_point(size=3,stroke=0.1,shape=21,alpha=0.6)+
  xlim(-6,6)+
  geom_vline(xintercept = 0,color="darkgrey")+
  
  scale_fill_manual(values=c("#440154FF","#29AF7FFF","#FDE725FF"))+
  ylab("-log10(p)")+
  xlab("Synaptosome/Homogenate")+
  
  theme_light()



#######
df3=filter(df2, category != "synapse")
df3$category=factor(df3$category,levels=c("presynapse","postsynapse"))
df3=arrange(df3,category)

ggplot(df3,
       aes(x=Estimate,y=-log10(p),fill=category))+
  geom_point(size=3,stroke=0.1,shape=21,alpha=0.6)+
  xlim(-6,6)+
  geom_vline(xintercept = 0,color="darkgrey")+
  
  scale_fill_manual(values=c("#29AF7FFF","#FDE725FF"))+
  ylab("-log10(p)")+
  xlab("Synaptosome/Homogenate")+
  
  theme_light()

##################

####add density plot
ggplot(df,
       aes(x=Estimate,fill=category))+
  geom_density(alpha=0.6)+
  xlim(-6,6)+
  geom_vline(xintercept = 0,color="darkgrey")+
  theme_light()+
  scale_fill_manual(values=c("grey","#440154FF","#29AF7FFF","#FDE725FF"))+
  xlab("log2 (Synaptosome/Homogenate Ratio average / subject)")


ggplot(df,
       aes(x=Estimate,fill=group))+
  geom_density(alpha=0.6)+
  xlim(-6,6)+
  geom_vline(xintercept = 0,color="darkgrey")+
  geom_vline(xintercept = 1.0,color="black")+
  geom_vline(xintercept = 0.7,color="black")+
  theme_light()+
  scale_fill_manual(values=c("grey","#440154FF"))+
  xlab("log2 (Synaptosome/Homogenate Ratio average / subject)")



######################################################################################################
###################################################
###################################################
# now see if you can plot individual subject comparisons using the residuals
res.Tissue.full=res.Tissue
#res.Tissue=res.Tissue.full

#########
# major plex effect
# consider only the plexes that are consistent with each other (based on MDS plot)

# make new meta file where only these plexes are kept
#Hom outlier plexes are 
# F53, F65 = F1
# F54, F66, F67 = F2
# F30, F27 = F3
# F26, F29, F40 = F4

#tmp=meta.samples
#tmp$PLEX=ifelse(tmp$PLEX == "F53" | tmp$PLEX == "F65","F1",
 #               ifelse(tmp$PLEX == "F54" | tmp$PLEX == "F66" | tmp$PLEX == "F67", "F2",
  #                     ifelse(tmp$PLEX == "F30" | tmp$PLEX == "F27", "F3",
   #                           "F4")))

#meta.samples=filter(tmp, PLEX == "F2" | PLEX == "F4")

#res.Tissue=res.Tissue[,colnames(res.Tissue) %in% rownames(meta.samples)]

########


########
# calculate the within subject enrichment values
x_columns=grep("\\.h$",colnames(res.Tissue),value=TRUE)

# Loop through the ".h" columns, subtract the corresponding ".s" column and create a new column with the result
for (col in x_columns) {
  y_col <- gsub("\\.h$", ".s", col)
  new_col <- gsub("\\.h$", ".diff", col)
  res.Tissue[[new_col]] <- res.Tissue[[y_col]] - res.Tissue[[col]]
}


#subset Syn and HOM ratio into new dataframe
df.PRT.Ratio=res.Tissue[,c(57:84)]
#df.PRT.Ratio=res.Tissue[,c(33:48)]
#df.PRT.Ratio.tmp=res.Tissue[,c(57:84)]

# add protein names to df
df.PRT.Ratio=rownames_to_column(df.PRT.Ratio)
names(df.PRT.Ratio)[1]="protein"

#### add gene names to df
df.PRT.Ratio=left_join(df.PRT.Ratio,protein.map,by="protein")
df.PRT.Ratio=df.PRT.Ratio[,c(1,30,31,2:29)]
#df.PRT.Ratio=df.PRT.Ratio[,c(1,18,19,2:17)]




select=filter(df.PRT.Ratio, 
              gene == "SYP" |
                gene == "DLG4" |
                gene == "CAMK2A" |
                gene == "MAP2" |
                gene == "SYT1" |
                gene == "DDX3Y" |
                gene == "SNRPD1" |
                gene == "CASP14"|
                gene == "SRSF6"|
                gene == "PVALB"|
                protein == "Q9BYB0"|
                gene == "GPHN"|
                gene == "SLC17A7"|
                gene == "GRM5"|
                gene == "GFAP")

##################################################
# add new value for 
#loop through new values
new_value="SHANK3"
row_index=which(select$protein == "Q9BYB0")

select[row_index, "gene"]=new_value
##################################################
select$gene=factor(select$gene,levels=c("SRSF6","SNRPD1","DDX3Y","CASP14","GFAP","PVALB","GPHN","SHANK3","DLG4","CAMK2A","GRM5","MAP2","SLC17A7","SYT1","SYP"))

select=left_join(select,RETissue.set[,c(2,11)],by="protein")

select.long = gather(select, subject, value, S1866.diff:S1672.diff,factor_key = T)

#####
####add box plot

ggplot(select.long,
       aes(x=value,y=gene,fill=group))+
  geom_boxplot(alpha=0.6)+
  xlim(-8,8)+
  geom_vline(xintercept = 0,color="darkgrey")+
  theme_light()+
  scale_fill_manual(values=c("grey","#440154FF"))

#plotMDS(res.Tissue)
