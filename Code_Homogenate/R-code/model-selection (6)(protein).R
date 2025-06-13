rm(list=ls()); gc()
options(stringsAsFactors = F)

tissue="HOMOGENATE-PROTEIN"

#### this piece of code performs the model selection procedure
# the user supplies the possible covariates to consider in the model 
# the program then determines whicb of these covariates should be added to the model based on BIC
# this is sequential process where at each step the covariate that improves the most features (proteins)
# gets added to the model until less than 5% of the features are improved.
# all decisions are based on log2 transformed data
# model selection is based on samples only (pools are ignored)
# THIS CODE IS FOR A MIXED LINEAR MODEL.

#### libraries used
require(data.table)           # fast version of read.table
require(readxl)               # used for reading Excel format files
require(viridis)              # colorblind friendly color palette
require(matrixStats)          # to get statistics on matrices
require(doParallel)           # parallel computing
require(lmerTest)


### load the VAE imptued data
load(paste0("WorkData/NORMALIZED-",tissue,".RData"))

#### fucntion
determineBIC<-function(feature,Y,XX,covariates,random.effect){
  # function performs a linear regression given a set of covariates and return the BIC for the model
  if(!is.null(random.effect)){
    model.formula<-formula(paste("Y['",feature,"',] ~ ",paste(c(covariates),collapse=" + "),paste0(" + (1|",random.effect,")"),sep=""))
    re.regr=lmer(model.formula,data=XX,REML = F,control = lmerControl(optimizer="bobyqa"))
  }else{
    model.formula<-formula(paste("Y['",feature,"',] ~ ",paste(c(covariates),collapse=" + "),sep=""))
    re.regr=lm(model.formula,data=XX)
  }
  return(BIC(re.regr))
}



# possible model covariates
#meta.samples$AGEyr2=meta.samples$AGEyr^2 # create a quadratic age effect
possible.covariates=c("PLEX","SEX","DX","AGEyr","AGEyr2","PMIhr","RACE")
# random effect (usualyy something like PAIR for a mixed model, when NULL a regular linear model is used)
random.effect="PAIR" 

# initial run uses and intercept only
null.model.covariates="DX"    # DX will be included in the model because it is the covariate of main interest
model.proteins=rownames(norm.prt) # list of peptides to model



# here the cut-offs in improvement in BIC are determined for the different covariates. These values are equivalent to a likelihood ratio test with p < 0.005
# continuous and 2 factor covariates
re.model.selection=NULL   # results of the model selection step
delta.cut=rep(-(qchisq(0.005,df=1,lower.tail=F)-log(nrow(meta.samples))),length(possible.covariates)) 
names(delta.cut)=possible.covariates
# for the factor covariates
delta.cut["PLEX"]=-(qchisq(0.005,df=nlevels(meta.samples$PLEX)-1,lower.tail=F)-log(nrow(meta.samples)))
delta.cut["RACE"]=-(qchisq(0.005,df=nlevels(meta.samples$RACE)-1,lower.tail=F)-log(nrow(meta.samples)))


repeat{
  # model to compare to
  system.time(re.null.BIC<-unlist(mclapply(model.proteins,determineBIC,Y=log2(norm.prt[model.proteins,samples]),XX=meta.samples[samples,],
                                           covariates=null.model.covariates,random.effect=random.effect,
                                           mc.cores=1,mc.preschedule = T,mc.silent = T)))
  
  if(length(grep("1",null.model.covariates)) != 0 ){null.model.covariates = NULL}   # remove the intercept from the null model covariates
  # this is no longer needed
  
  # alternative models
  delta.BIC=NULL
  for(covariate in setdiff(possible.covariates,null.model.covariates)){
    alt.model.covariates=c(null.model.covariates,covariate)
    print(alt.model.covariates); flush.console()
    re.alt.BIC<-unlist(mclapply(model.proteins,determineBIC,Y=log2(norm.prt[model.proteins,samples]),XX=meta.samples[samples,],
                                covariates=alt.model.covariates,random.effect=random.effect,
                                mc.cores=1,mc.preschedule = T,mc.silent = T))
    delta.BIC=cbind(delta.BIC,re.alt.BIC - re.null.BIC < delta.cut[covariate]) #  improved TRUE/FALSE
  }
  colnames(delta.BIC)=setdiff(possible.covariates,null.model.covariates)
  rownames(delta.BIC)=model.proteins
  (improved.BIC=as.data.frame(t(colSums(delta.BIC))))
  (improved.BIC=cbind.data.frame(improved.BIC,data.frame(selected=colnames(improved.BIC)[which.max(improved.BIC)]))) # covariates that most improved the model
  (improved.BIC=cbind.data.frame(improved.BIC,data.frame(improved=sum(rowSums(delta.BIC) > 0))))      # number of features where an added covariate improved the model
  (improved.BIC=cbind.data.frame(improved.BIC,data.frame(not.improved=sum(rowSums(delta.BIC) == 0)))) # number of features not showing improvement
  
  # update the null model
  null.model.covariates=c(null.model.covariates,improved.BIC$selected) 
  
  # select the features that showed improvement
  model.proteins=rownames(delta.BIC)[rowSums(delta.BIC) > 0]
  
  # store the results from this iteration
  if(is.null(re.model.selection)){
    re.model.selection=improved.BIC
  }else{
    # add the results 
    re.model.selection=rbind.data.frame(re.model.selection,re.model.selection[nrow(re.model.selection),])
    re.model.selection[nrow(re.model.selection),]=NA
    re.model.selection[nrow(re.model.selection),colnames(improved.BIC)]=improved.BIC
  }
  print(re.model.selection); flush.console()
  if(max(improved.BIC[1:(ncol(improved.BIC)-3)]) < nrow(norm.prt)*0.05 | nrow(re.model.selection) == length(possible.covariates)){break}
}

# some clean-up
selected.model=null.model.covariates
if(nrow(re.model.selection) < length(possible.covariates) | re.model.selection[nrow(re.model.selection),re.model.selection$selected[nrow(re.model.selection)]] < nrow(norm.prt)*0.05){
  selected.model=null.model.covariates[-length(null.model.covariates)]
}
selected.model
re.model.selection

### save the results with the input
if(!is.null(random.effect)){
  save(meta.samples,meta.pep,norm.prt,selected.model,re.model.selection,file=paste0("WorkData/LMER-MODEL-SELECTION-",tissue,".RData"))
}else{
  save(meta.samples,meta.pep,norm.prt,selected.model,re.model.selection,file=paste0("WorkData/LM-MODEL-SELECTION-",tissue,".RData"))
}

### look for additional models including interactions with DX that might be of interest

re.null.BIC=unlist(mclapply(rownames(norm.prt),determineBIC,Y=log2(norm.prt[,samples]),XX=meta.samples[samples,],
                            covariates=selected.model,random.effect=random.effect,
                            mc.cores=64,mc.preschedule = T,mc.silent = T))
re.alt.BIC=unlist(mclapply(rownames(norm.prt),determineBIC,Y=log2(norm.prt[,samples]),XX=meta.samples[samples,],
                           covariates=c(selected.model,"AGEyr2:DX"),random.effect=random.effect,
                           mc.cores=64,mc.preschedule = T,mc.silent = T))
sum(re.alt.BIC - re.null.BIC < delta.cut["AGEyr2"])  # 4 < .05*nrow(norm.prt) = 128.9



# no need to add interactions at this point
(selected.model=c(selected.model))

# neither meets the criteria of 0.05*nrow(norm.prt) (1644.1)
# stay with the main effects model

### save the results with the input
if(!is.null(random.effect)){
  save(meta.samples,meta.pep,norm.prt,selected.model,re.model.selection,samples,file=paste0("WorkData/LMER-MODEL-SELECTION-",tissue,".RData"))
}else{
  save(meta.samples,meta.pep,norm.prt,selected.model,re.model.selection,samples,file=paste0("WorkData/LM-MODEL-SELECTION-",tissue,".RData"))
}