rm(list=ls()); gc()
options(stringsAsFactors = F)

tissue="HOM-SYN-COMP-POOLS"

############################################
### This versions is designed for protein level work
### Start with VAE.full where pooled controls were imputed along with samples
### Retain only pools and NT subjects present in both Hom and Syn datasets
### Remove outlier pooled controls
### Retain only unique peptides that correspond to overlapping proteins in both Hom and Syn datasets
### Roll up to protein level and combine dataframes
### Perform internal reference scaling (pools) normalization on the observed scale
### save workdata and perform model selection in next script


#####################################

#### libraries used
require(data.table)           # fast version of read.table
require(readxl)               # used for reading Excel format files
require(viridis)              # colorblind friendly color palette
require(matrixStats)          # to get statistics on matrices
require(limma)


###################################################
#load imputed values that were imputed with pooled controls

####################
### Synaptosome data
load("Analysis_code/Synaptosome_enrichment/WorkData/VAE-SYNAPTOSOME-PEPTIDE.RData")

#prepare synaptosome files for merging with homogenate data
meta.pep.syn=meta.pep

vae.syn=vae.full

meta.syn=filter(meta.full, DX == "CTL" | DX == "POOL") # 46 x 12
meta.syn$ID.Tissue=paste(meta.syn$ID,".s",sep="")
meta.syn$Tissue=rep("Synaptosome",each=46)

####################
### Homogenate data
load("Analysis_code/Synaptosome_enrichment/WorkData/VAE-HOMOGENATE-PEPTIDE.RData")

meta.pep.hom=meta.pep

vae.hom=vae.full

meta.hom=filter(meta.full, DX == "CTL" | DX == "POOL") # 49 x 12
meta.hom=meta.hom[rownames(meta.hom) %in% rownames(meta.syn), ] # 46 x 12
#remove outlier pools
meta.hom=filter(meta.hom, ID != "P3008")
meta.hom=filter(meta.hom, ID != "P4011") # 44 x 12
meta.hom$ID.Tissue=paste(meta.hom$ID,".h",sep="")
meta.hom$Tissue=rep("Homogenate",each=44)

####################
#####
# clean vae.pep dfs with only retained subjects
vae.hom=vae.hom[ ,colnames(vae.hom) %in% rownames(meta.hom)] # 35103 x 44
vae.syn=vae.syn[ ,colnames(vae.syn) %in% rownames(meta.syn)] # 29311 x 46

#####
# swap out rownames on meta subject files and vae.pep dfs

rownames(meta.hom)=meta.hom$ID.Tissue
rownames(meta.syn)=meta.syn$ID.Tissue

#####
# swap out colnames of vae.pep dfs with new IDs
colnames(vae.hom)=rownames(meta.hom)
colnames(vae.syn)=rownames(meta.syn)

####################
#####
# remove nonunique peptides

#Hom
i.multi=grep(";",meta.pep.hom$ROLLUP) # 2221
if(length(i.multi) > 0){meta.pep.hom=meta.pep.hom[-i.multi,]} # 32,882
length(unique(meta.pep.hom$ROLLUP)) # 4727

#Syn
i.multi=grep(";",meta.pep.syn$ROLLUP) # 1820
if(length(i.multi) > 0){meta.pep.syn=meta.pep.syn[-i.multi,]} # 27,491
length(unique(meta.pep.syn$ROLLUP)) # 4287

#####
# find overlapping proteins
tmp=intersect(meta.pep.syn$ROLLUP,meta.pep.hom$ROLLUP) # 3852
meta.pep.syn=meta.pep.syn[meta.pep.syn$ROLLUP %in% tmp, ] # 26,921
meta.pep.hom=meta.pep.hom[meta.pep.hom$ROLLUP %in% tmp, ] # 31,369

#####
# retain only peptides that correspond to overlapping proteins
vae.hom=vae.hom[rownames(meta.pep.hom), ] # 31369 x 28
vae.syn=vae.syn[rownames(meta.pep.syn), ] # 26921 x 28


##############
#####
# roll up peptides to proteins

#Hom
vae.hom.prt=aggregate(vae.hom~meta.pep.hom$ROLLUP,FUN="sum")
rownames(vae.hom.prt)=vae.hom.prt[,1]
vae.hom.prt=as.matrix(vae.hom.prt[,-1])
dim(vae.hom.prt) # 3852 x 44

#Syn
vae.syn.prt=aggregate(vae.syn~meta.pep.syn$ROLLUP,FUN="sum")
rownames(vae.syn.prt)=vae.syn.prt[,1]
vae.syn.prt=as.matrix(vae.syn.prt[,-1])
dim(vae.syn.prt) # 3852 x 46


##############
#####
# combine Hom and Syn dfs
meta.samples=rbind(meta.hom,meta.syn) # 90 x 14

vae.prt=cbind(vae.hom.prt,vae.syn.prt) # 3852 x 90

######################################################################
### turn categorical variables into factors
meta.samples$SEX=factor(c("Male","Female")[meta.samples$SEX], levels=c("Male","Female"))

meta.samples$RACE=factor(c("White","Black or African-American","Unknown")[meta.samples$RACE],levels=c("White","Black or African-American","Unknown"),labels=c("CAU","AFR","UNK"))

meta.samples$PLEX=factor(meta.samples$PLEX, levels=unique(meta.samples$PLEX))

meta.samples$Tissue=factor(meta.samples$Tissue, levels=unique(meta.samples$Tissue))



### some additional information that will be used throughout the pipeline
samples=rownames(meta.samples)[meta.samples$ROLE == "SAMPLE"]
plexes=unique(meta.samples$PLEX)

# now some things for for plotting
plex.border=cumsum(table(meta.samples$PLEX))
plex.border=plex.border[-length(plex.border)]+.5
plex.mid=c(0,cumsum(table(meta.samples$PLEX)))
plex.mid=.5+(plex.mid[-length(plex.mid)]+plex.mid[-1])/2
names(plex.mid)=plexes

pchTissue=c(15,17);names(pchTissue)=unique(meta.samples$Tissue)
colPLEX=viridis(length(unique(meta.samples$PLEX))); names(colPLEX)=unique(meta.samples$PLEX)



##### PLOTS
pdf(paste0("Analysis_Code/Synaptosome_enrichment/Plots/total-intensity-ROLLEDUP-",tissue,".pdf"),height=8,width=12)
plot(colSums(vae.prt,na.rm=T)/1e6,col=colPLEX[meta.samples$PLEX],pch=pchTissue[meta.samples$Tissue],
     las=1,xaxt="n",xlab="",ylab="")
abline(v=plex.border,lty=3,col="grey50")
title(ylab="total-intensity (millions)",line = 2.5)
title(main="Total intensity for samples by plex")
axis(side = 1, at=plex.mid, labels = paste("Plex", 1:10))
legend("bottomleft",legend = levels(meta.samples$Tissue),pch=pchTissue,col="grey75",bty="n")
dev.off()


pdf(paste0("Analysis_Code/Synaptosome_enrichment/Plots/MDS-ROLLEDUP-",tissue,".pdf"),height = 8, width =12)
par(mfrow=c(1,2))
par(oma=c(0,0,2,0))
par(mar=c(5,5,1,2))
plotMDS(log2(vae.prt),col=colPLEX[meta.samples$PLEX],las=1,main="",dim.plot = c(1,2))
legend("topleft",legend=paste("Plex",1:10),pch=15,pt.cex = 2,col=colPLEX,bty="n")
plotMDS(log2(vae.prt),col=colPLEX[meta.samples$PLEX],las=1,main="",dim.plot = c(3,4))
mtext(text="MDS of normalized log2(intensity) data",side=3,outer=T,cex=1.5,line=0)

dev.off()


######################################################
#########
#now do IRS TMT scaling at the protein level
#####

PC_flag <- grepl("POOL", meta.samples$ROLE, fixed=TRUE)
Sample_flag <- !grepl("POOL", meta.samples$ROLE, fixed=TRUE)

#recode plex values to be 1 through 10
t=meta.samples
t$PLEX2=ifelse(t$PLEX == "F54",1,
               ifelse(t$PLEX == "F53",2,
                      ifelse(t$PLEX == "F65",3,
                             ifelse(t$PLEX == "F66",4,
                                    ifelse(t$PLEX == "F67",5,
                                           ifelse(t$PLEX == "F29", 6,
                                                  ifelse(t$PLEX == "F30",7,
                                                         ifelse(t$PLEX == "F27",8,
                                                                ifelse(t$PLEX == "F26",9,
                                                                       10)))))))))


t$PLEX=t$PLEX2

tmp_meta.test=t

#################
# IRS
irs <- c()
for(i in seq(10)){
  flag_tmp <- PC_flag & (tmp_meta.test$PLEX==i)
  df_tmp <- vae.prt[,flag_tmp]
  df_tmp[df_tmp==0] <- NA
  irs <- cbind(irs,rowMeans(df_tmp,na.rm = T))
  colnames(irs)[i] <- paste("pc",i,sep="")
}


irs_scaling <- apply(irs, 1, function(x) exp(mean(log(x),na.rm=T)))

irs1 <- irs_scaling/irs
dim(irs1) # 3852 x 10


for(i in seq(10)){
  flag_tmp <- !PC_flag & (tmp_meta.test$PLEX==i)
  vae.prt[,flag_tmp] <- vae.prt[,flag_tmp]*irs1[,i]
}

norm.prt=vae.prt

#####################################################################

##### PLOTS
pdf(paste0("Analysis_Code/Synaptosome_enrichment/Plots/total-intensity-IRSnorm-",tissue,".pdf"),height=8,width=12)
plot(colSums(norm.prt,na.rm=T)/1e6,col=colPLEX[meta.samples$PLEX],pch=pchTissue[meta.samples$Tissue],
     las=1,xaxt="n",xlab="",ylab="")
abline(v=plex.border,lty=3,col="grey50")
title(ylab="total-intensity (millions)",line = 2.5)
title(main="Total intensity for samples by plex after IRS norm")
axis(side = 1, at=plex.mid, labels = paste("Plex", 1:10))
legend("bottomleft",legend = levels(meta.samples$Tissue),pch=pchTissue,col="grey75",bty="n")
dev.off()


pdf(paste0("Analysis_Code/Synaptosome_enrichment/Plots/MDS-IRSnorm-",tissue,".pdf"),height = 8, width =12)
par(mfrow=c(1,2))
par(oma=c(0,0,2,0))
par(mar=c(5,5,1,2))
plotMDS(log2(norm.prt),col=colPLEX[meta.samples$PLEX],las=1,main="",dim.plot = c(1,2))
legend("topright",legend=paste("Plex",1:10),pch=15,pt.cex = 2,col=colPLEX,bty="n")
plotMDS(log2(norm.prt),col=colPLEX[meta.samples$PLEX],las=1,main="",dim.plot = c(3,4))
mtext(text="MDS of IRS normalized log2(intensity) data",side=3,outer=T,cex=1.5,line=0)

dev.off()


######################################################
### remove pooled controls and replot
meta.samples=filter(meta.samples,DX == "CTL") # 56 x 14
norm.prt=norm.prt[ ,colnames(norm.prt) %in% rownames(meta.samples)] # 3852 x 56

#redefine sample names before saving
# ^^^^^^^^^^^^^^^^^^

##### PLOTS
pdf(paste0("Analysis_Code/Synaptosome_enrichment/Plots/total-intensity-IRSnormnoPOOLs-",tissue,".pdf"),height=8,width=12)
plot(colSums(norm.prt,na.rm=T)/1e6,col=colPLEX[meta.samples$PLEX],pch=pchTissue[meta.samples$Tissue],
     las=1,xaxt="n",xlab="",ylab="")
abline(v=plex.border,lty=3,col="grey50")
title(ylab="total-intensity (millions)",line = 2.5)
title(main="Total intensity for samples by plex after IRS norm")
axis(side = 1, at=plex.mid, labels = paste("Plex", 1:10))
legend("bottomleft",legend = levels(meta.samples$Tissue),pch=pchTissue,col="grey75",bty="n")
dev.off()


pdf(paste0("Analysis_Code/Synaptosome_enrichment/Plots/MDS-IRSnormnoPOOLs-",tissue,".pdf"),height = 8, width =12)
par(mfrow=c(1,2))
par(oma=c(0,0,2,0))
par(mar=c(5,5,1,2))
plotMDS(log2(norm.prt),col=colPLEX[meta.samples$PLEX],las=1,main="",dim.plot = c(1,2))
legend("bottomright",legend=paste("Plex",1:10),pch=15,pt.cex = 2,col=colPLEX,bty="n")
plotMDS(log2(norm.prt),col=colPLEX[meta.samples$PLEX],las=1,main="",dim.plot = c(3,4))
mtext(text="MDS of IRS normalized log2(intensity) data",side=3,outer=T,cex=1.5,line=0)

dev.off()


###################################################
# adjust any additional files
meta.pep=meta.pep.syn



#########################
#### save the resulting data
save(meta.samples,meta.pep,norm.prt,samples,plexes,plex.border,plex.mid,pchDX,colPLEX,file=paste0("Analysis_Code/Synaptosome_enrichment/WorkData/IRS-NORMALIZED-",tissue,".RData"))


########################################

