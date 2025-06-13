#updated combining results to include synaptosome data with 
#age and dx interaction in the model

rm(list=ls()); gc()
options(stringsAsFactors = F)

require(data.table)
require(WriteXLS)
require(doParallel)
require(viridis)

### load the peptide-protein-gene mapping information
load("WorkData/peptide-protein-gene-maps.RData")

### read the results information for the peptides
tissue="SYNAPTOSOME-PEPTIDE"

load(paste0("LMER-REGRESSION-INT-",tissue,".RData"))

peptide.MLM.DX=reDX
peptide.MLM.DX$protein=peptide.map[peptide.MLM.DX$peptide,"protein"]
peptide.MLM.DX$gene=peptide.map[peptide.MLM.DX$peptide,"gene"]
peptide.MLM.DX$ensg=peptide.map[peptide.MLM.DX$peptide,"ensg"]

peptide.MLM.DX.AGE=reDX.AGE
names(peptide.MLM.DX.AGE)[1]="peptide"
peptide.MLM.DX.AGE$protein=peptide.map[peptide.MLM.DX.AGE$peptide,"protein"]
peptide.MLM.DX.AGE$gene=peptide.map[peptide.MLM.DX.AGE$peptide,"gene"]
peptide.MLM.DX.AGE$ensg=peptide.map[peptide.MLM.DX.AGE$peptide,"ensg"]

load(paste0("LMER-REGRESSION-",tissue,".RData"))

peptide.MLM.AGE=reAGE
peptide.MLM.AGE$protein=peptide.map[peptide.MLM.AGE$peptide,"protein"]
peptide.MLM.AGE$gene=peptide.map[peptide.MLM.AGE$peptide,"gene"]
peptide.MLM.AGE$ensg=peptide.map[peptide.MLM.AGE$peptide,"ensg"]

RE.PEPTIDE=list()
RE.PEPTIDE[["MLM.DX"]]=peptide.MLM.DX
RE.PEPTIDE[["MLM.AGE"]]=peptide.MLM.AGE
RE.PEPTIDE[["MLM.DX.AGE"]]=peptide.MLM.DX.AGE


#WriteXLS(RE.PEPTIDE,paste0("SYNAPTOSOME-PEPTIDE-INT.xlsx"))
#write.csv(peptide.MLM.DX,file="SYNAPTOSOME-PEPTIDE-INT.csv")




#####################################################
#################################################
###########################################
#proteins
### read the results information for the proteins
tissue="SYNAPTOSOME-PROTEIN"

load(paste0("WorkData/LMER-REGRESSION-",tissue,".RData"))

protein.MLM.DX=reDX
protein.MLM.DX$gene=protein.map[protein.MLM.DX$protein,"gene"]
protein.MLM.DX$ensg=protein.map[protein.MLM.DX$protein,"ensg"]

protein.MLM.AGE=reAGE
protein.MLM.AGE$gene=protein.map[protein.MLM.AGE$protein,"gene"]
protein.MLM.AGE$ensg=protein.map[protein.MLM.AGE$protein,"ensg"]

protein.MLM.AGE.noint=reAGE.noint
protein.MLM.AGE.noint$gene=protein.map[protein.MLM.AGE.noint$protein,"gene"]
protein.MLM.AGE.noint$ensg=protein.map[protein.MLM.AGE.noint$protein,"ensg"]

protein.MLM.DX.AGE=reDX.AGE
protein.MLM.DX.AGE$gene=protein.map[protein.MLM.DX.AGE$protein,"gene"]
protein.MLM.DX.AGE$ensg=protein.map[protein.MLM.DX.AGE$protein,"ensg"]


RE.PROTEIN=list()
RE.PROTEIN[["MLM.DX"]]=protein.MLM.DX
RE.PROTEIN[["MLM.AGE"]]=protein.MLM.AGE
RE.PROTEIN[["MLM.AGE.noint"]]=protein.MLM.AGE.noint
RE.PROTEIN[["MLM.DX.AGE"]]=protein.MLM.DX.AGE

#WriteXLS(RE.PROTEIN,paste0("Results/",tissue,".xlsx"))

save(norm.pep,norm.prt,meta.samples,RE.PEPTIDE,RE.PROTEIN,
     colPLEX,pchDX,plex.border,plex.mid,plexes,pools,samples,selected.model,file="WorkData/RESULTS-SYNAPTOSOME-updateINT.RData")




### make a plot compare results from PTT and MLM.DX
pdf("Plots/comparing-PTT-MLM-DX-SYNAPTOSOME.pdf",height=12,width=12)
par(mfrow=c(2,2))
par(mar=c(3,3,3,1))
par(oma=c(3,3,1,1))
plot(peptide.PTT$effect,peptide.MLM.DX$Estimate,las=1,pch=19,col=viridis(5)[2])
mtext("Peptide",side=2,line=3)
mtext("Fold Change",side=3,line=1)
abline(a=0,b=1,col="grey25")

plot(-log10(peptide.PTT$p.emperical),-log10(peptide.MLM.DX$p),las=1,pch=19,col=viridis(5)[2],xlim=range(-log10(peptide.PTT$p.emperical)),ylim=range(-log10(peptide.PTT$p.emperical)))
mtext("-log10(P)",side=3,line=1)
abline(a=0,b=1,col="grey25")

plot(protein.PTT$effect,protein.MLM.DX$Estimate,las=1,pch=19,col=viridis(5)[3])
mtext("Protein",side=2,line=3)
abline(a=0,b=1,col="grey25")

plot(-log10(protein.PTT$p.emperical),-log10(protein.MLM.DX$p),las=1,pch=19,col=viridis(5)[3],xlim=range(-log10(protein.PTT$p.emperical)),ylim=range(-log10(protein.PTT$p.emperical)),
     xlab="paired t-test",ylab="")
abline(a=0,b=1,col="grey25")

mtext("Paired T-Test",side=1,line=1,outer=T)
mtext("Mixed Linear Model",side=2,line=1,outer=T)
dev.off()