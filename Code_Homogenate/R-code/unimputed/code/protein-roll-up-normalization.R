rm(list=ls()); gc()
options(stringsAsFactors = F)

require(data.table)
tissue="HOMOGENATE-PROTEIN"

load("../WorkData/QC-HOMOGENATE-PEPTIDE.RData")

### select the samples and complete peptides
qc.pep=qc.pep[,samples]
meta.samples=meta.samples[samples,]  # 62

qc.pep=qc.pep[rowSums(is.na(qc.pep)) == 0,]
meta.pep=meta.pep[rownames(qc.pep),]  # 11,777

### selet the peptides linked to a single protein
remove.pep=rownames(meta.pep)[grep("; ",meta.pep$ROLLUP)]   # 1005
meta.pep=meta.pep[!(rownames(meta.pep) %in% remove.pep),]   # 10,772
qc.pep=qc.pep[rownames(meta.pep),]

# roll-up into proteins
qc.prt=aggregate(qc.pep~meta.pep$ROLLUP,FUN="sum")
rownames(qc.prt)=qc.prt[,1]
qc.prt=as.matrix(qc.prt[,-1])

# SL-normalize
sl.mean=colMeans(qc.prt)
all.mean=mean(sl.mean)
norm.prt=t(t(qc.prt)*(all.mean/sl.mean))

### save the necesary information
save(norm.prt,meta.samples,meta.pep,samples,file=paste0("WorkData/NORMALIZED-",tissue,".RData"))
