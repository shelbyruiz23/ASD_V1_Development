rm(list=ls()); gc()
options(stringsAsFactors = F)

tissue="SYNAPTOSOME-PEPTIDE"
########################################################
#### THIS VERSION SKIPS THE SURROGATE VARIABLES AND 
#### TAKES ITS INPUT FROM THE MODEL SELECTION STEP
########################################################

# this part performs the linear regression on the data given the selected covariates and sva
# at the end the use can select the results for the covariate of interest. In most cases this will be DX.
# as one woudll expect, regression will only be performed for the samples
# this part uses parallel programming, which will speed up calcualtions. This only works on linux and mac machines



#### libraries used
require(data.table)           # fast version of read.table
require(readxl)               # used for reading Excel format files
require(viridis)              # colorblind friendly color palette
require(matrixStats)          # to get statistics on matrices
require(doParallel)           # parallel computing
require(limma)
require(sva)
require(qvalue)
require(lmerTest)


### load the VAE imptued data
load(paste0("WorkData/MODEL-SELECTION-",tissue,".RData"))

selected.model=c(selected.model,"DX:AGEyr")

# functions needed
### linear regression using lm
lmerRegression <- function(feature,Y,XX,covariates,random.effect){
  # function performs a linear regression given a set of covariates and return the BIC for the model
  model.formula<-formula(paste("Y['",feature,"',] ~ ",paste(c(covariates),collapse=" + "),paste0(" + (1|",random.effect,")"),sep=""))
  re.lmer=lmer(model.formula,data=XX,REML = T,control = lmerControl(optimizer="bobyqa"))
  return(summary(re.lmer)$coefficient)
}
  

### select the results for a specific effect 
selectResults<-function(feature,results,effect){
  irow=match(paste0(effect),rownames(results[[feature]]))
  re.out=data.frame(feature,t(results[[feature]][irow,]))
  return(re.out)
}

# run the mixedlinear regression, this step saves the covariate coefficients for each feature
system.time(reLMER <- mclapply(rownames(norm.pep),lmerRegression,Y=log2(norm.pep[,samples]),XX=meta.samples[samples,],
                                       covariates=selected.model,random.effect="PAIR",mc.cores=1,mc.preschedule=T))
names(reLMER)=rownames(norm.pep)

# select the results for DX
system.time(reDX<-as.data.frame(rbindlist(mclapply(names(reLMER),selectResults,results=reLMER,effect="DXASD",mc.cores = 1,mc.preschedule = T))))
colnames(reDX)=c("peptide","Estimate","Std.Error","df","t.value","p") # give friendly names to the columns

1-qvalue(reDX$p)$pi0  # 0.259834
reDX$q=qvalue(reDX$p)$qvalues # use qvalue to account for multiple testing

###
#DX & AGE - model
sum(reDX$q < 0.05) # 2
sum(reDX$q < 0.10) # 11
sum(reDX$q < 0.20) # 759
###

###
#DX,AGE,DX:AGE - model
sum(reDX$q < 0.05) # 1 
sum(reDX$q < 0.10) # 736
sum(reDX$q < 0.20) # 2451
###

pdf(paste0("Plots/DX-mixed-model-p-values-",tissue,".pdf"),height=8,width=12)
  hist(reDX$p,breaks=seq(0,1,.01),las=1,xlab="P",ylab="Number of peptides",main="",col=viridis(3)[2])
dev.off()

# select the results for AGE
system.time(reAGE<-as.data.frame(rbindlist(mclapply(names(reLMER),selectResults,results=reLMER,effect="AGEyr",mc.cores = 1,mc.preschedule = T))))
colnames(reAGE)=c("peptide","Estimate","Std.Error","df","t.value","p") # give friendly names to the columns

1-qvalue(reAGE$p)$pi0  # .1492732
reAGE$q=qvalue(reAGE$p)$qvalues # use qvalue to account for multiple testing

###
#DX & AGE - model
sum(reAGE$q < 0.05) # 3162
sum(reAGE$q < 0.10) # 5107
sum(reAGE$q < 0.20) # 8426
###

###
#DX,AGE,DX:AGE - model
sum(reAGE$q < 0.05) #347 
sum(reAGE$q< 0.10) # 667
sum(reAGE$q < 0.20) # 1322
###

pdf(paste0("Plots/AGE-mixed-model-p-values-",tissue,".pdf"),height=8,width=12)
hist(reAGE$p,breaks=seq(0,1,.01),las=1,xlab="P",ylab="Number of peptides",main="",col=viridis(5)[4])
dev.off()

# select the results for the interaction
system.time(reDX.AGE<-as.data.frame(rbindlist(mclapply(names(reLMER),selectResults,results=reLMER,effect="DXASD:AGEyr",mc.cores = 1,mc.preschedule = T))))
colnames(reDX.AGE)=c("protein","Estimate","Std.Error","df","t.value","p") # give friendly names to the columns

1-qvalue(reDX.AGE$p)$pi0 # .3225723
reDX.AGE$q=qvalue(reDX.AGE$p)$qvalues # use qvalue to account for multiple testing


###
#DX,AGE,DX:AGE - model
sum(reDX.AGE$q < 0.05) # 0
sum(reDX.AGE$q < 0.10) # 0
sum(reDX.AGE$q < 0.20) # 2525
###

pdf(paste0("Plots/DX.AGE-mixed-model-p-values-",tissue,".pdf"),height=8,width=12)
hist(reDX.AGE$p,breaks=seq(0,1,.01),las=1,xlab="P",ylab="Number of proteins",main="",col=viridis(5)[2])
dev.off()




save(reLMER,reDX,reAGE,reDX.AGE,meta.samples,meta.pep,norm.pep,selected.model,samples,pools,plexes,plex.border,plex.mid,pchDX,colPLEX,file=paste0("WorkData/LMER-REGRESSION-INT-",tissue,".RData"))


write.csv(reDX,"Syn_peptide_MLM.DX_INT.csv")
