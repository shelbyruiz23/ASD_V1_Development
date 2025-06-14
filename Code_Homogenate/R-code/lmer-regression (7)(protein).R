rm(list=ls()); gc()
options(stringsAsFactors = F)

tissue="HOMOGENATE-PROTEIN"
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
#require(lme4)
require(lmerTest)
require(sva)
require(qvalue)


### load thenormalized data and the informaation on the model
load(paste0("WorkData/MODEL-SELECTION-",tissue,".RData"))
selected.model=c(selected.model,"DX:AGEyr") # only for testing adding interaction of age and dx

########################################################
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
system.time(reLMER <- mclapply(rownames(norm.prt),lmerRegression,Y=log2(norm.prt[,samples]),XX=meta.samples[samples,],
                                       covariates=selected.model,random.effect="PAIR",mc.cores=1,mc.preschedule=T))


names(reLMER)=rownames(norm.prt)

# select the results for DX
system.time(reDX<-as.data.frame(rbindlist(mclapply(names(reLMER),selectResults,results=reLMER,effect="DXASD",mc.cores = 1,mc.preschedule = T))))
colnames(reDX)=c("protein","Estimate","Std.Error","df","t.value","p") # give friendly names to the columns


1-qvalue(reDX$p)$pi0
reDX$q=qvalue(reDX$p)$qvalues # use qvalue to account for multiple testing


sum(reDX$q < 0.05) # 32
sum(reDX$q < 0.10) # 230
sum(reDX$q < 0.20) # 658

##############################
# for forcing interaction of DX:AGE into model
sum(reDX$q < 0.05) # 0
sum(reDX$q < 0.10) # 0
sum(reDX$q < 0.20) # 7

# select the results for DX:AGE interaction
system.time(reDX_AGE<-as.data.frame(rbindlist(mclapply(names(reLMER),selectResults,results=reLMER,effect="DXASD:AGEyr",mc.cores = 1,mc.preschedule = T))))
colnames(reDX_AGE)=c("protein","Estimate","Std.Error","df","t.value","p") # give friendly names to the columns


1-qvalue(reDX_AGE$p)$pi0
reDX_AGE$q=qvalue(reDX_AGE$p)$qvalues # use qvalue to account for multiple testing


# for forcing interaction of DX:AGE into model

# DX effect
sum(reDX$q < 0.05) # 0
sum(reDX$q < 0.10) # 0
sum(reDX$q < 0.20) # 7

# DX:AGE effect
sum(reDX_AGE$q < 0.05) # 0
sum(reDX_AGE$q < 0.10) # 0
sum(reDX_AGE$q < 0.20) # 0

##############################

pdf(paste0("Plots/DX-mixed-model-p-values-",tissue,".pdf"),height=8,width=12)
hist(reDX$p,breaks=seq(0,1,.01),las=1,xlab="P",ylab="Number of proteins",main="",col=viridis(5)[3])
dev.off()

# select the results for AGE
system.time(reAGE<-as.data.frame(rbindlist(mclapply(names(reLMER),selectResults,results=reLMER,effect="AGEyr",mc.cores = 4,mc.preschedule = T))))
colnames(reAGE)=c("protein","Estimate","Std.Error","df","t.value","p") # give friendly names to the columns

1-qvalue(reAGE$p)$pi0
reAGE$q=qvalue(reAGE$p)$qvalues # use qvalue to account for multiple testing


sum(reAGE$q < 0.05) # 1261
sum(reAGE$q < 0.10) # 1672
sum(reAGE$q < 0.20) # 2300

pdf(paste0("Plots/AGE-mixed-model-p-values-",tissue,".pdf"),height=8,width=12)
hist(reAGE$p,breaks=seq(0,1,.01),las=1,xlab="P",ylab="Number of proteins",main="",col=viridis(5)[4])
dev.off()


save(reLMER,reDX,reAGE,meta.samples,meta.pep,norm.prt,selected.model,samples,pools,plexes,plex.border,plex.mid,pchDX,colPLEX,file=paste0("WorkData/LMER-REGRESSION-",tissue,".RData"))

#save(reLMER,reDX_AGE,meta.samples,meta.pep,norm.prt,selected.model,samples,pools,plexes,plex.border,plex.mid,pchDX,colPLEX,file=paste0("LMER-REGRESSION-DXAGEint-",tissue,".RData"))
