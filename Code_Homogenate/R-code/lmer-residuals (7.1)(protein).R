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
require(sva)
require(qvalue)
require(lmerTest)


### load thenormalized data and the informaation on the model
load(paste0("../WorkData/MODEL-SELECTION-",tissue,".RData"))
selected.model=selected.model

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


system.time(res.prt <- mclapply(rownames(norm.prt),lmerResiduals,Y=log2(norm.prt[,samples]),XX=meta.samples[samples,],
                                covariates=selected.model,random.effect="PAIR",add.back=c("DXASD","AGEyr"),mc.cores=1,mc.preschedule=T))
res.prt=as.data.frame(rbindlist(res.prt))
rownames(res.prt)=rownames(norm.prt)
colnames(res.prt)=colnames(norm.prt)

#################
#residual=res.prt
#residual_AGE=res.prt
#residual_DX=res.prt
#residual_DXandAGE=res.prt

#save(residual,
 #    residual_DX,
  #   residual_DXandAGE,
   #  meta.samples,
    # norm.prt,
    # file="Analysis_Code/Code_Homogenate/WorkData/residuals.RData")
############


WriteXLS::WriteXLS(res.prt,paste0(tissue,"-RESIDUALS+DX+AGE+INTERACTION_11-22-2022.xlsx"),SheetNames = tissue,row.names = T)

res.prt=as.matrix(res.prt)
### plotting
pch=c(15,19); names(pch)=c("CTL","ASD")
col=viridis(5)[c(2,4)]; names(col)=c("CTL","ASD")

#### example
prts=c("P29218","Q53HC0","Q9P0S9","P25713","Q9BUT1","Q9P2U7","Q08209","Q9NQ66","Q9NTI2")
pdf("examples.pdf",height=12,width=12)
par(mfrow=c(3,3))
par(mar=c(2,2,1,1))
par(oma=c(3,3,1,1))
for(prt in prts){

  plot(meta.samples$AGEyr,res.prt[prt,],pch=pch[meta.samples$DX],col=col[meta.samples$DX],ylim=c(-1,1),las=1,xlab="",ylab="",xaxt="n",yaxt="n")
  if(grep(prt,prts)%%3 == 1){
    axis(2,at=seq(-1,1,.5),las=1)  
  }else{
    axis(2,at=seq(-1,1,.5),labels = rep("",5),las=1)
  }
  if(grep(prt,prts)> 6){
    axis(1,at=seq(5,30,5),las=1)  
  }else{
    axis(1,at=seq(5,30,5),labels = rep("",6),las=1)
  }
  

  legend("topleft",legend=prt,bty="n")
  
  # use lm to get the coefficients for DX, Age, and interaction
  model.formula<-formula(paste("res.prt['",prt,"',] ~ ",paste(c("DX*AGEyr"),collapse=" + "),sep=""))
  betas=lm(model.formula,data=meta.samples)$coefficients
  # regression line for NT
  a.ctl=betas["(Intercept)"]
  b.ctl=betas["AGEyr"]
  # regression line for ASD
  a.asd=betas["(Intercept)"]+betas["DXASD"]
  b.asd=betas["AGEyr"]+betas["DXASD:AGEyr"]
  
  # plot both lines
  abline(a=a.ctl,b=b.ctl,col=col["CTL"],lwd=2)
  abline(a=a.asd,b=b.asd,col=col["ASD"],lwd=2)
  }
mtext("Age (yr)",side=1,line=1,outer=T)
mtext(expression(paste("log"[2],"(PLEX adjusted abundance)")),side=2,line=1,outer=T)

dev.off()