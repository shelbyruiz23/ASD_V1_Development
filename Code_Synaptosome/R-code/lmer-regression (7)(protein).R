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


### load the normalized data and the information on the model
load(paste0("WorkData/MODEL-SELECTION-",tissue,".RData"))
selected.model=c(selected.model,"DX:AGEyr")

# for calculating the effect of age in synaptosome data do not include interaction (DX:AGEyr) in the model
# calculate age with main effects only

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
                                       covariates=selected.model,random.effect="PAIR",mc.cores=64,mc.preschedule=T))
names(reLMER)=rownames(norm.prt)

############################
# select the results for DX
system.time(reDX<-as.data.frame(rbindlist(mclapply(names(reLMER),selectResults,results=reLMER,effect="DXASD",mc.cores = 4,mc.preschedule = T))))
colnames(reDX)=c("protein","Estimate","Std.Error","df","t.value","p") # give friendly names to the columns

1-qvalue(reDX$p)$pi0 # 0.265
reDX$q=qvalue(reDX$p)$qvalues # use qvalue to account for multiple testing


sum(reDX$q < 0.05) # 93
sum(reDX$q < 0.10) # 304
sum(reDX$q < 0.20) # 711

pdf(paste0("Plots/DX-mixed-model-p-values-",tissue,".pdf"),height=8,width=12)
hist(reDX$p,breaks=seq(0,1,.01),las=1,xlab="P",ylab="Number of proteins",main="",col=viridis(5)[3])
dev.off()

############################
# select the results for AGE
# if needed, reselect the model to include only the main effects
system.time(reAGE<-as.data.frame(rbindlist(mclapply(names(reLMER),selectResults,results=reLMER,effect="AGEyr",mc.cores = 4,mc.preschedule = T))))
colnames(reAGE)=c("protein","Estimate","Std.Error","df","t.value","p") # give friendly names to the columns

1-qvalue(reAGE$p)$pi0 # .146
reAGE$q=qvalue(reAGE$p)$qvalues # use qvalue to account for multiple testing


sum(reAGE$q < 0.05) # 144
sum(reAGE$q < 0.10) # 246
sum(reAGE$q < 0.20) # 409

pdf(paste0("Plots/AGE-mixed-model-p-values-",tissue,".pdf"),height=8,width=12)
hist(reAGE$p,breaks=seq(0,1,.01),las=1,xlab="P",ylab="Number of proteins",main="",col=viridis(5)[4])
dev.off()

# select the results for the interaction
system.time(reDX.AGE<-as.data.frame(rbindlist(mclapply(names(reLMER),selectResults,results=reLMER,effect="DXASD:AGEyr",mc.cores = 4,mc.preschedule = T))))
colnames(reDX.AGE)=c("protein","Estimate","Std.Error","df","t.value","p") # give friendly names to the columns

1-qvalue(reDX.AGE$p)$pi0 # .436
reDX.AGE$q=qvalue(reDX.AGE$p)$qvalues # use qvalue to account for multiple testing


sum(reDX.AGE$q < 0.05) # 0
sum(reDX.AGE$q < 0.10) # 305
sum(reDX.AGE$q < 0.20) # 1109

pdf(paste0("Plots/DX.AGE-mixed-model-p-values-",tissue,".pdf"),height=8,width=12)
hist(reDX.AGE$p,breaks=seq(0,1,.01),las=1,xlab="P",ylab="Number of proteins",main="",col=viridis(5)[2])
dev.off()


###import age results without interaction and save to the workdata
#reAGE.noint=as.data.frame(read_excel("C:/Users/Shelby_macdonaldlab/OneDrive - University of Pittsburgh/MacDonaldLab_Data/ASD_Development_March2023/Paper_files/Supplemental_tables/Results/Protein/Results_protein.xlsx", sheet='Synaptosome_MLM.AGE'))
#reAGE.noint=reAGE.noint[,c(1:7)]

save(reLMER,reDX,reAGE,reAGE.noint,reDX.AGE,meta.samples,meta.pep,norm.prt,selected.model,samples,pools,plexes,plex.border,plex.mid,pchDX,colPLEX,file=paste0("WorkData/LMER-REGRESSION-",tissue,".RData"))






##### Example
feature="Q8TBG9" # "Q9P2U7"
covariates=selected.model[-4]
Y=log2(norm.prt[,samples])
XX=meta.samples[samples,]
random.effect="PAIR"

model.formula<-formula(paste("Y['",feature,"',] ~ ",paste(c(covariates),collapse=" + "),paste0(" + (1|",random.effect,")"),sep=""))
re.noint=lmer(model.formula,data=XX,REML = T,control = lmerControl(optimizer="bobyqa"))

beta=c(re.noint@beta,0)
names(beta)=c(colnames(re.noint@pp$X),"PLEXF29")
y=Y[feature,]-beta[paste0("PLEX",meta.samples[samples,"PLEX"])]

df=data.frame(pair=unique(meta.samples$PAIR))
rownames(df)=df$pair   # 28 pairs
df$ctl.id=rownames(meta.samples)[match(paste(df$pair,"CTL",sep="-"),paste(meta.samples$PAIR,meta.samples$DX,sep="-"))]
df$case.id=rownames(meta.samples)[match(paste(df$pair,"ASD",sep="-"),paste(meta.samples$PAIR,meta.samples$DX,sep="-"))]

df$ctl.plex=meta.samples[df$ctl.id,"PLEX"]
df$case.plex=meta.samples[df$case.id,"PLEX"]

df$ctl.age=meta.samples[df$ctl.id,"AGEyr"]
df$case.age=meta.samples[df$case.id,"AGEyr"]


df$ctl.res=y[df$ctl.id]
df$case.res=y[df$case.id]

plot(meta.samples[samples,"AGEyr"],y[samples])
points(df$case.age,df$case.res,pch=19,col="blue")
points(df$ctl.age,df$ctl.res,pch=19,col="red")
abline(lm(df$case.res~df$case.age)$coef,col="blue",lwd=2)
abline(lm(df$ctl.res~df$ctl.age)$coef,col="red",lwd=2)


covariates=selected.model
model.formula<-formula(paste("Y['",feature,"',] ~ ",paste(c(covariates),collapse=" + "),paste0(" + (1|",random.effect,")"),sep=""))
re.int=lmer(model.formula,data=XX,REML = T,control = lmerControl(optimizer="bobyqa"))

summary(lm(df$case.res~df$case.age))
summary(lm(df$ctl.res~df$ctl.age))



beta=c(re.int@beta)
names(beta)=c(colnames(re.int@pp$X))
x=re.int@pp$X

beta=beta[-c(1:2,8)]
x=x[,-c(1:2,8)]
res=Y[feature,]-x%*%beta
y=-beta[paste0("PLEX",meta.samples[samples,"PLEX"])]
boxplot(res~meta.samples$DX)
