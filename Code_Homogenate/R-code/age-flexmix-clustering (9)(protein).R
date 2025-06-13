rm(list=ls()); gc()
options(stringsAsFactors = F)

# this portion performs a tracjectory analysis for the age effects

require(data.table)           # fast version of read.table
require(readxl)               # used for reading Excel format files
require(viridis)              # colorblind friendly color palette
require(matrixStats)          # to get statistics on matrices
require(doParallel)           # parallel computing
require(limma)
require(sva)
require(flexmix)
require(lmerTest)
require(qvalue)
library(xlsx)


tissue="HOMOGENATE-PROTEIN"

##########################################################################
#setwd
setwd("C:/Users/Shelby_macdonaldlab/OneDrive - University of Pittsburgh/MacDonaldLab_Data/ASD_Development_March2023")
##########################################################################


##########################################################################
### load thenormalized data and the informaation on the model
load("Bert_ASDpipeline_dontalter_2022_09/ASD-Development-Homogenate-Package-Sept2022/WorkData/MODEL-SELECTION-HOMOGENATE-PROTEIN.RData")
##########################################################################
#selected.model=selected.model

##########################################################################
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
##########################################################################




### calcualte the residuals adding back the effect of age
system.time(RESage <- mclapply(rownames(norm.prt),lmerResiduals,Y=log2(norm.prt[,samples]),XX=meta.samples[samples,],
                               covariates=selected.model,random.effect="PAIR",add.back=c("AGEyr"),mc.cores=1,mc.preschedule=T))
RESage=as.data.frame(rbindlist(RESage))
colnames(RESage)=colnames(norm.prt)
rownames(RESage)=rownames(norm.prt)



# load the data to 
#load("WorkData/SURROGATE-VARIABLES-HOMOGENATE.RData")


# set up the dataframe to use for flexmix.
asd=rownames(meta.samples)[meta.samples$DX == "ASD"]
ctl=rownames(meta.samples)[meta.samples$DX == "CTL"]
proteins=rownames(norm.prt)

# use only the the controls
df=data.frame(protein=rep(proteins,each=length(ctl)),sample=rep(ctl,times=length(proteins)))
df[,c("AGEyr")]=meta.samples[df$sample,c("AGEyr")]
df[,"RESage"]=RESage[cbind(df$protein,df$sample)]

# run flexmix to find the mixture distribution, this can take some time
re.flexmix <- flexmix(RESage ~ AGEyr  | protein, k = 3, 
                      data = df, control = list(tolerance = 10^-5),model=FLXMRglmfix(varFix = T))

(n.tbl=table(clusters(re.flexmix))/length(ctl))
#HOM Age_only redo (CTL)
#1(dec)    2(stable)    3 (inc)
#1078 3195  454





df$cluster=clusters(re.flexmix)
prt.clusters=df[!duplicated(df[,c("protein","cluster")]),c("protein","cluster")]


df.param=parameters(re.flexmix)
df.param=df.param[rownames(df.param) != "sigma",]
df.age=df.param[grep("coef.AGE",rownames(df.param)),]
cov.matrix=data.frame(intercept=1,age=seq(min(df$AGEyr),max(df$AGEyr),1))  
y.age=as.matrix(cov.matrix)%*%as.matrix(df.param)  # intercept plus age effect

pdf("Age/Age_only_redo/plots/HOMOGENATE-PROTEIN-AGE-mixture-modeling_04-12-2022.pdf",height=8,width=8)
plot(NULL,NULL,xlim=range(cov.matrix$age),ylim=range(y.age),las=1,xlab="Age",ylab="Effect")
abline(h=0,col="grey50",lty=3)
for(i in 1:ncol(y.age)){
  points(cov.matrix$age,y.age[,i],type="l",col=viridis(9)[c(2,4,6)][i],lwd=2,lty=i)
}

plus.minus=c("-","+","+");names(plus.minus)=c("-1","0","1")
legend("topleft",legend=c(sprintf("abnd = %4.2f %s %4.3f * age; n.prt: %4i",round(df.param[1,1],2),plus.minus[as.character(sign(df.param[2,1]))],abs(df.param[2,1]),n.tbl[1]),
                          sprintf("abnd = %4.2f %s %4.3f * age; n.prt: %4i",round(df.param[1,2],2),plus.minus[as.character(sign(df.param[2,2]))],abs(df.param[2,2]),n.tbl[2]),
                          sprintf("abnd = %4.2f %s %4.3f * age; n.prt: %4i",round(df.param[1,3],2),plus.minus[as.character(sign(df.param[2,3]))],abs(df.param[2,3]),n.tbl[3])),
       lty=1:ncol(y.age),lwd = 2,pch=NULL,col=viridis(9)[c(2,4,6)][1:ncol(y.age)],bty="n")
dev.off()


RE=list(parameters=as.data.frame(df.param),clustering=prt.clusters)

#save as excell sheets
for (i in names(RE)){
  xlsx::write.xlsx(RE[i], file="Age/Age_only_redo/Results/HOMOGENATE-PROTEIN-AGE-mixture-modeling_04-12-2023.xlsx", sheetName=paste(i), append=T)
}

##########################################################################
#perform flexmix clustering for ASD samples as well

# use only the asd
df2=data.frame(protein=rep(proteins,each=length(asd)),sample=rep(asd,times=length(proteins)))
df2[,c("AGEyr")]=meta.samples[df2$sample,c("AGEyr")]
df2[,"RESage"]=RESage[cbind(df2$protein,df2$sample)]

# run flexmix to find the mixture distribution, this can take some time
re.flexmix2 <- flexmix(RESage ~ AGEyr  | protein, k = 3, 
                       data = df2, control = list(tolerance = 10^-5),model=FLXMRglmfix(varFix = T))

(n.tbl=table(clusters(re.flexmix2))/length(asd))
#HOM Age only redo (ASD)
#1 (stable)    2 (dec)   3 (inc)
#3266 1037  424



df2$cluster=clusters(re.flexmix2)
prt.clusters2=df2[!duplicated(df2[,c("protein","cluster")]),c("protein","cluster")]


df.param2=parameters(re.flexmix2)
df.param2=df.param2[rownames(df.param2) != "sigma",]
df.age2=df.param2[grep("coef.AGE",rownames(df.param2)),]
cov.matrix2=data.frame(intercept=1,age=seq(min(df2$AGEyr),max(df2$AGEyr),1))  
y.age2=as.matrix(cov.matrix2)%*%as.matrix(df.param2)  # intercept plus age effect

pdf("Age/Age_only_redo/plots/HOMOGENATE-ASD-PROTEIN-AGE-mixture-modeling_04-12-2023.pdf",height=8,width=8)
plot(NULL,NULL,xlim=range(cov.matrix2$age),ylim=range(y.age2),las=1,xlab="Age",ylab="Effect")
abline(h=0,col="grey50",lty=3)
for(i in 1:ncol(y.age2)){
  points(cov.matrix2$age,y.age2[,i],type="l",col=viridis(9)[c(2,4,6)][i],lwd=2,lty=i)
}

plus.minus=c("-","+","+");names(plus.minus)=c("-1","0","1")
legend("topleft",legend=c(sprintf("abnd = %4.2f %s %4.3f * age; n.prt: %4i",round(df.param2[1,1],2),plus.minus[as.character(sign(df.param2[2,1]))],abs(df.param2[2,1]),n.tbl[1]),
                          sprintf("abnd = %4.2f %s %4.3f * age; n.prt: %4i",round(df.param2[1,2],2),plus.minus[as.character(sign(df.param2[2,2]))],abs(df.param2[2,2]),n.tbl[2]),
                          sprintf("abnd = %4.2f %s %4.3f * age; n.prt: %4i",round(df.param2[1,3],2),plus.minus[as.character(sign(df.param2[2,3]))],abs(df.param2[2,3]),n.tbl[3])),
       lty=1:ncol(y.age2),lwd = 2,pch=NULL,col=viridis(9)[c(2,4,6)][1:ncol(y.age2)],bty="n")
dev.off()


RE2=list(parameters=as.data.frame(df.param2),clustering=prt.clusters2)

#save as excell sheets
for (i in names(RE2)){
  xlsx::write.xlsx(RE[i], file="Age/Age_only_redo/Results/HOMOGENATE-ASD-PROTEIN-AGE-mixture-modeling_04-12-2023.xlsx", sheetName=paste(i), append=T)
}



##########################################################################


save(re.flexmix,meta.samples,meta.pep,norm.prt,selected.model,RESage,RE,RE2,samples,pools,plexes,plex.border,plex.mid,pchDX,colPLEX,file="Age/Age_only_redo/WorkData/TRAJECTORY-HOMOGENATE.RData")

load("Age/Age_only_redo/WorkData/TRAJECTORY-HOMOGENATE.RData")

##### calculate the slope for both asd and ctl

linearRegression<-function(feature,yy,xx){
  lm.formula=as.formula(paste0("yy['",feature,"',]","~","xx"))
  re.lm=summary(lm(lm.formula))
  a=re.lm$coef[1,1]
  b=re.lm$coef[2,1]
  se=re.lm$coef[2,2]
  return(data.frame(feature,a,b,se))
}
lrCTL=as.data.frame(rbindlist(lapply(prt.clusters$protein,linearRegression,yy=as.matrix(RESage[,ctl]),xx=meta.samples[ctl,"AGEyr"])))
lrASD=as.data.frame(rbindlist(lapply(prt.clusters$protein,linearRegression,yy=as.matrix(RESage[,asd]),xx=meta.samples[asd,"AGEyr"])))

prt.clusters=RE$clustering
rownames(prt.clusters)=prt.clusters$protein
prt.clusters[lrASD$feature,c("a.asd","b.asd","se.asd")]=lrASD[,c("a","b","se")]
prt.clusters[lrCTL$feature,c("a.ctl","b.ctl","se.ctl")]=lrCTL[,c("a","b","se")]

########################################################################################
plot(prt.clusters[prt.clusters$cluster == 3,"b.ctl"],prt.clusters[prt.clusters$cluster == 3,"b.asd"],xlim=c(-0.08,0.12),ylim=c(-0.08,0.12))
re.lm=lm(prt.clusters[prt.clusters$cluster == 3,"b.asd"]~prt.clusters[prt.clusters$cluster == 3,"b.ctl"])
abline(a=0,b=1)
#abline(re.lm$coef)

plot(prt.clusters[prt.clusters$cluster == 1,"b.ctl"],prt.clusters[prt.clusters$cluster == 1,"b.asd"],xlim=c(-0.08,0.12),ylim=c(-0.08,0.12))
re.lm=lm(prt.clusters[prt.clusters$cluster == 1,"b.asd"]~prt.clusters[prt.clusters$cluster == 1,"b.ctl"])
abline(a=0,b=1)
#abline(re.lm$coef)

plot(prt.clusters[prt.clusters$cluster == 2,"b.ctl"],prt.clusters[prt.clusters$cluster == 2,"b.asd"],xlim=c(-0.08,0.12),ylim=c(-0.08,0.12))
re.lm=lm(prt.clusters[prt.clusters$cluster == 2,"b.asd"]~prt.clusters[prt.clusters$cluster == 2,"b.ctl"])
abline(a=0,b=1)
#abline(re.lm$coef)




########################################################################################

t.test(prt.clusters[prt.clusters$cluster == 1,"b.asd"],prt.clusters[prt.clusters$cluster == 1,"b.ctl"],paired = T)
t.test(prt.clusters[prt.clusters$cluster == 2,"b.asd"],prt.clusters[prt.clusters$cluster == 2,"b.ctl"],paired = T)
t.test(prt.clusters[prt.clusters$cluster == 3,"b.asd"],prt.clusters[prt.clusters$cluster == 3,"b.ctl"],paired = T)

hist(prt.clusters[prt.clusters$cluster == 2,"b.asd"]-prt.clusters[prt.clusters$cluster == 2,"b.ctl"])


regr.coef=aggregate(cbind(a.asd,b.asd,a.ctl,b.ctl)~cluster,FUN="mean",data=as.data.frame(prt.clusters))

pdf("Age/Age_only_redo/plots/HOMOGENATE-AGE-TRAJ-ASD-NT-2023-04-11.pdf",height=8,width=8)
plot(NULL,NULL,xlim=c(4,32),ylim=c(-.5,.75),las=1,xlab="Age",ylab="Effect")
abline(h=0,col="grey50")
for(i in 1:3){
  abline(a=regr.coef[i,"a.ctl"],b=regr.coef[i,"b.ctl"],lty=2,col=viridis(9)[c(2,5,8)][i],lwd=2)
  abline(a=regr.coef[i,"a.asd"],b=regr.coef[i,"b.asd"],lty=1,col=viridis(9)[c(2,5,8)][i],lwd=2)
}
legend("topleft",legend=c("NT","ASD"),lty=c(2,1),lwd=2,bty="n")
legend("bottomleft",legend=paste0(c("stable","decreasing","increasing"),"; n = ",n.tbl),col=viridis(9)[c(2,5,8)],lty=1,lwd=2,bty="n")
dev.off()


MLM=as.data.frame(read_excel("Bert_ASDpipeline_dontalter_2022_09/ASD-Development-Homogenate-Package-Sept2022/Results/HOMOGENATE-PROTEIN.xlsx",sheet="MLM.DX"))
rownames(MLM)=MLM$protein

########################################################################################
posterior=re.flexmix@posterior$scaled[seq(1,nrow(re.flexmix@posterior$scaled),31),]
prt.clusters[,c("pp1","pp2","pp3")]=posterior
t.test(prt.clusters[prt.clusters$cluster == 2 & prt.clusters$pp2 > .75 & prt.clusters$pp2 <= .90,"b.ctl"],
       prt.clusters[prt.clusters$cluster == 2 & prt.clusters$pp2 > .75 & prt.clusters$pp2 <= .90,"b.asd"],paired = T)
#p=0.001928

t.test(prt.clusters[prt.clusters$cluster == 2 & prt.clusters$pp2 > .90 & prt.clusters$pp2 <= 1,"b.ctl"],
       prt.clusters[prt.clusters$cluster == 2 & prt.clusters$pp2 > .90 & prt.clusters$pp2 <= 1,"b.asd"],paired = T)
#p=0.1628

t.test(prt.clusters[prt.clusters$cluster == 3 & prt.clusters$pp3 > .75 & prt.clusters$pp3 <= .90,"b.ctl"],
       prt.clusters[prt.clusters$cluster == 3 & prt.clusters$pp3 > .75 & prt.clusters$pp3 <= .90,"b.asd"],paired = T)
#p=0.391

t.test(prt.clusters[prt.clusters$cluster == 3 & prt.clusters$pp3 > .90 & prt.clusters$pp3 <= 1,"b.ctl"],
       prt.clusters[prt.clusters$cluster == 3 & prt.clusters$pp3 > .90 & prt.clusters$pp3 <= 1,"b.asd"],paired = T)
#p=0.81

########################################################################################

prt.clusters$delta=prt.clusters$b.asd-prt.clusters$b.ctl
prt.clusters$se.delta=sqrt(prt.clusters$se.asd^2+prt.clusters$se.ctl^2)
prt.clusters$p.delta=2*pt(-abs(prt.clusters$delta/prt.clusters$se.delta),df=30)
prt.clusters$q.delta=qvalue(prt.clusters$p.delta)$qvalue

table(prt.clusters$cluster[prt.clusters$q.delta < 0.05])

prt.clusters[MLM$protein,"dx.fc"]=MLM$Estimate
prt.clusters[MLM$protein,"dx.q"]=MLM$q
prt.clusters[MLM$protein,"gene"]=MLM$gene


#add asd specific age clustering to prt.clusters dataframe
RE.ASD=as.data.frame(RE2$clustering)
names(RE.ASD)[2]="cluster.ASD"
prt.clusters[RE.ASD$protein,"cluster.ASD"]=RE.ASD$cluster.ASD
#reorder
prt.clusters=prt.clusters[,c(1,2,19,3:18)]

#replace cluster numbers with identities
#this may be specific to each separate flex mix run
prt.clusters["cluster"][prt.clusters["cluster"] == 2] <- "stable"
prt.clusters["cluster"][prt.clusters["cluster"] == 3] <- "increasing"
prt.clusters["cluster"][prt.clusters["cluster"] == 1] <- "decreasing"

prt.clusters["cluster.ASD"][prt.clusters["cluster.ASD"] == 1] <- "stable"
prt.clusters["cluster.ASD"][prt.clusters["cluster.ASD"] == 3] <- "increasing"
prt.clusters["cluster.ASD"][prt.clusters["cluster.ASD"] == 2] <- "decreasing"

#HOM Age_only redo (CTL)
#1(dec)    2(stable)    3 (inc)
#1078 3195  454

#HOM Age only redo (ASD)
#1 (stable)    2 (dec)   3 (inc)
#3266 1037  424


xlsx::write.xlsx(prt.clusters,"Age/Age_only_redo/Results/HOMOGENATE-PROTEIN-PRTClusters.xlsx")
save(re.flexmix,meta.samples,meta.pep,norm.prt,selected.model,RESage,RE,RE2,samples,pools,plexes,plex.border,plex.mid,pchDX,colPLEX,prt.clusters, file="Age/Age_only_redo/WorkData/TRAJECTORY-HOMOGENATE.RData")

########################################################################################
