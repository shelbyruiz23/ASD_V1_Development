rm(list=ls()); gc()
options(stringsAsFactors = F)

tissue="HOM-SYN-COMP-POOLS"
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
require(lme4)


### load the normalized data and the information on the model

### load the IRS normalized data
load(paste0("Analysis_Code/Synaptosome_enrichment/WorkData/IRS-NORMALIZED-",tissue,".RData"))
### no model selection step


#################################

#if no model selection, choose the most likely differences
selected.model=c("Tissue","AGEyr")

# this determines whether LM (random.effect = NULL) or LMER (random.effect = some variable)results are used 
# ID is necessary to use as a random effect because Hom and Syn should be within subject comparisons
#random.effect=NULL
random.effect="ID"

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

#################
######
#####

effects=c("TissueSynaptosome","AGEyr") # effects of interest as identified in the regression coefficients
RE=list()
# select the results for TissueSynaptosome
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

########################################
#### save RETissue and add gene names
RETissue=RE$TissueSynaptosome

### add q value
RETissue$q=qvalue(RETissue$p)$qvalues

load("Analysis_Code/Code_Synaptosome/WorkData/peptide-protein-gene-maps.RData")
RETissue=left_join(RETissue,protein.map,by="protein")


########################################
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
    pdf(paste0("Analysis_Code/Synaptosome_enrichment/Plots/LMER-",effect,"-p-values-",tissue,".pdf"),height=8,width=12)
  }else{
    pdf(paste0("Analysis_Code/Synaptosome_enrichment/Plots/LM-",effect,"-p-values-",tissue,".pdf"),height=8,width=12)
  }
  hist(RE[[effect]]$p,breaks=seq(0,1,.01),las=1,xlab="P",ylab="Number of proteins",main="",col=colors[effect])
  dev.off()
}
re.summary

############################################################################################################################
# return residuals with only the effect of TissueSynaptosome added back
# residuals will be used to calculate within subject enrichment distributions of select proteins

###########
# functions needed
lmerResiduals <- function(feature,Y,XX,covariates,random.effect,add.back){
  # function performs a linear regression given a set of covariates and returns the residuals
  # add.back are the effect as specified in the lmer coefficient that you want to backin (use NULL if you don't want to anything back)
  model.formula<-formula(paste("Y['",feature,"',] ~ ",paste(c(covariates),collapse=" + "),paste0(" + (1|",random.effect,")"),sep=""))
  re.lmer=lmer(model.formula,data=XX,REML = T,control = lmerControl(optimizer="bobyqa"))
  re.residuals=residuals(re.lmer)
  re.coef=summary(re.lmer)$coefficients
  
  # add back the effect of interest
  if(!is.null(add.back)){
    # set up the matrix for the covariates
    covMatrix=model.matrix(formula(paste("~",paste(selected.model,collapse=" + "))),data=XX)
    # find the regression coefficients of interest
    beta=as.matrix(re.coef[add.back,"Estimate"])
    # determine the adjustment
    cov.adjust=covMatrix[,add.back]%*%beta
    # add this effect back in
    re.residuals=re.residuals+cov.adjust
  }
  
  # add back the add.back effect
  return(data.frame(t(re.residuals)))
}


### calcualte the residuals adding back the effect of age and age:dx interaction
system.time(res.Tissue <- mclapply(rownames(norm.prt),lmerResiduals,Y=log2(norm.prt[,samples]),XX=meta.samples[samples,],
                               covariates=selected.model,random.effect="ID",add.back=c("TissueSynaptosome"),mc.cores=1,mc.preschedule=T))

res.Tissue=as.data.frame(rbindlist(res.Tissue))
rownames(res.Tissue)=rownames(norm.prt)
colnames(res.Tissue)=colnames(norm.prt)
#########################################################

############################################################################################################################
# save workdata
# norm.prt
# meta.pep
# meta.samples
# selected.model
# RE
# RETissue
# res.Tissue
# protein.map

save(norm.prt,
     meta.pep,
     meta.samples,
     selected.model,
     RE,
     RETissue,
     res.Tissue,
     protein.map,
     samples,
     plexes,
     plex.border,
     plex.mid,
     pchDX,
     colPLEX,
     file=paste0("Analysis_Code/Synaptosome_enrichment/WorkData/LMER-REGRESSION-",tissue,".RData"))

#write.csv(RETissue,"RETissue_lmer_results.csv")
