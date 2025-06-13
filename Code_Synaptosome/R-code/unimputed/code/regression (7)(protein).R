rm(list=ls()); gc()
options(stringsAsFactors = F)

tissue="SYNAPTOSOME-PROTEIN"
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
library(tidyverse)

# this determines whether LM (random.effect = NULL) or LMER (random.effect = some variable)results are used 
random.effect="PAIR"
### load the unimputed data
if(!is.null (random.effect)){
  load(paste0("WorkData/LMER-MODEL-SELECTION-",tissue,".RData"))
}else{  
  load(paste0("WorkData/LM-MODEL-SELECTION-",tissue,".RData"))
}
# functions needed
### linear regression using lm
regression <- function(feature,Y,XX,covariates,random.effect){
  # function performs a linear regression given a set of covariates and return the BIC for the model
  if(!is.null(random.effect)){
    model.formula<-formula(paste("Y['",feature,"',] ~ ",paste(c(covariates),collapse=" + "),paste0(" + (1|",random.effect,")"),sep=""))
    re.lmer=lmer(model.formula,data=XX,REML = T,control = lmerControl(optimizer="bobyqa"))
  }else{
    model.formula<-formula(paste("Y['",feature,"',] ~ ",paste(c(covariates),collapse=" + "),sep=""))
    re.lmer=lm(model.formula,data=XX)
  }
  return(summary(re.lmer)$coefficient)
}
  

### select the results for a specific effect 
selectResults<-function(feature,results,effect){
  irow=match(paste0(effect),rownames(results[[feature]]))
  re.out=data.frame(feature,t(results[[feature]][irow,]))
  return(re.out)
}

# run the mixedlinear regression, this step saves the covariate coefficients for each feature step
system.time(reREGR <- mclapply(rownames(norm.prt),regression,Y=log2(norm.prt[,samples]),XX=meta.samples[samples,],
                                       covariates=selected.model,random.effect=random.effect,mc.cores=1,mc.preschedule=T))
names(reREGR)=rownames(norm.prt)
reREGR[[1]]

effects=c("DXASD","AGEyr","DXASD:AGEyr") # effects of interest as identified in the regression coefficients
RE=list()
# select the results for DX
for(effect in effects){
  print(effect)
  system.time(re<-as.data.frame(rbindlist(mclapply(names(reREGR),selectResults,results=reREGR,effect=effect,mc.cores = 1,mc.preschedule = T))))
  if(!is.null(random.effect)){
      colnames(re)=c("protein","Estimate","Std.Error","df","t.value","p") # give friendly names to the columns
  }else{
    colnames(re)=c("protein","Estimate","Std.Error","t.value","p") # give friendly names to the columns
  }
  RE[[effect]]=re
}

# summarize and plotting
re.summary=NULL
colors=viridis(length(effects)+2);colors=colors[-c(1,length(colors))];names(colors)=effects
for(effect in effects){
  RE[[effect]]$q=qvalue(RE[[effect]]$p)$qvalues # use qvalue to account for multiple testing
  re.summary=rbind.data.frame(re.summary,data.frame(effect,pi1=1-qvalue(RE[[effect]]$p)$pi0,
                                                    q.lt.05=sum(RE[[effect]]$q < 0.05),
                                                    q.lt.10=sum(RE[[effect]]$q < 0.10),
                                                    q.lt.20=sum(RE[[effect]]$q < 0.20)))
  
  if(!is.null(random.effect)){
    pdf(paste0("Plots/LMER-",effect,"-p-values-",tissue,".pdf"),height=8,width=12)
  }else{
    pdf(paste0("Plots/LM-",effect,"-p-values-",tissue,".pdf"),height=8,width=12)
  }
  hist(RE[[effect]]$p,breaks=seq(0,1,.01),las=1,xlab="P",ylab="Number of proteins",main="",col=colors[effect])
  dev.off()
}
re.summary

#pdf saving did not want to save with the colon in the file name
hist(RE[["DXASD:AGEyr"]]$p,breaks=seq(0,1,.01),las=1,xlab="P",ylab="Number of proteins",main="",col=colors["DXASD:AGEyr"])


#save individual effect regression results
RE.DX=RE$DXASD
RE.AGE=RE$AGEyr
RE.DX.AGE=RE$`DXASD:AGEyr`

#map proteins to genes
#load WorkData from peptide-protein-gene-maps (1.1)(peptide)
#ignore meta.pep file now (old)
RE.DX=left_join(RE.DX,protein.map,by="protein")
RE.AGE=left_join(RE.AGE,protein.map,by="protein")
RE.DX.AGE=left_join(RE.DX.AGE,protein.map,by="protein")



# save the results
if(!is.null(random.effect)){
  save(reREGR,RE.DX,RE.AGE,RE.DX.AGE,re.summary,meta.samples,samples,norm.prt,selected.model,file=paste0("WorkData/LMER-REGRESSION-",tissue,".RData"))
}else{
  save(reREGR,RE.DX,RE.AGE,RE.DX.AGE,re.summary,meta.samples,samples,norm.prt,selected.model,file=paste0("WorkData/LM-REGRESSION-",tissue,".RData"))
}


#output the results
write.csv(RE.DX,"Results/LMER-REGRESSION-DX-UNIMPUTED-HOMOGENATE-PROTEIN-2023-11-18.csv")
write.csv(RE.AGE,"Results/LMER-REGRESSION-AGE-UNIMPUTED-HOMOGENATE-PROTEIN-2023-11-18.csv")
write.csv(RE.DX.AGE,"Results/LMER-REGRESSION-DX-AGE-UNIMPUTED-HOMOGENATE-PROTEIN-2023-11-18.csv")
