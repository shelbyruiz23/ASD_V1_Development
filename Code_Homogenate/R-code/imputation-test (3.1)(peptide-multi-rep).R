rm(list=ls()); gc()
options(stringsAsFactors = F)

tissue="HOMOGENATE-PEPTIDE"

### this code performs the VAE imputation using softImpute to come up with the initialization values
### the code should not have to be changed other than when one wants to use another initialization process
### the user needs to specify the covariates to use in the VAE
#### libraries used
require(data.table)           # fast version of read.table
require(readxl)               # used for reading Excel format files
require(viridis)              # colorblind friendly color palette
require(matrixStats)          # to get statistics on matrices
require(softImpute)
require(limma)

#### copies of VAEIT and vae-python.R need to be in this directory

############## function
mseImputation <- function(lambda,r,sample.col,pool.col,Cna,thresh=1e-4,maxit=100){
  # function that finds the lambda for a given r using the minimum distance between the observed
  # observed average intensity and the imputed average intensity
  # lambda - lambda to use in the imputation
  # rank - rank to use for the svd decomposition
  # Cna - matrix with missing values    
  re.softimpute<-softImpute(Cna,rank.max = r,lambda=lambda,type="als",trace.it  = F,thresh = thresh, maxit = maxit)
  if(r > 1){
    Cudv=(re.softimpute$u%*%diag(re.softimpute$d))%*%t(re.softimpute$v)
  }else{
    Cudv=(re.softimpute$u*re.softimpute$d)%*%t(re.softimpute$v)
  }
  colnames(Cudv)=colnames(Cna)
  rownames(Cudv)=rownames(Cna)
  
  i.complete=which(rowSums(is.na(Cna)) == 0)
  mse=mean((rowMeans(Cudv[i.complete,sample.col],na.rm=T)-rowMeans(Cna[i.complete,pool.col],na.rm=T))^2)
  
  return(mse)
}


##### load the data
load(paste0("WorkData/QC-",tissue,".RData"))

# start VAE
#####################################
# now perform the VAE step
####################################

##### Instructions to perform setting up an environment
# conda install mamba -n base -c conda-forge
# conda config --add channels conda-forge
# mamba create -n py_aev python==3.8 -y   # run in shell
# conda activate py_aev
# mamba install -n py_env -c conda-forge -c esri tensorflow==2.6 tensorflow-probability==0.14 tensorflow-addons scikit-learn numpy -y




# start conda environment
### library reticulate is needed to run the python code
require(reticulate)
use_condaenv(condaenv = "py_vae")
# load the vae code
source("vae_python.r")

# set up the data
qc.pep=qc.pep[order(rowMeans(qc.pep,na.rm=T)),]


imputationFunction <- function(seed,qc.pep,meta.samples,samples,pools){
#  require(reticulate)
#  use_condaenv(condaenv = "py_vae")
  # load the vae code
 # source("vae_python.r")
  # simulate a dataset with missingness
  # order the peptides by average expression
#  seed=0
  set.seed(seed)
  
  # simulate a dataset with missingness
  Xall=qc.pep      # all data
  Xcomplete=qc.pep[rowSums(is.na(qc.pep)) == 0,] # complete data
  
  # select the pattern to simulate from Xall
  i.xall=sort(sample(nrow(Xall),nrow(Xcomplete)))
  # find the locations with missing values in Xcomplete for the selected peptides
  i.na=which(is.na(Xall[i.xall,]),arr.ind = F)
  # now create a matrix with missing vales based on this
  Xna=Xcomplete                               # start with the complete data
  Xna[i.na]=NA                                # add the missingnes


  #1) find starting values using softImpute
  
  # find the maximum value for lambda. This defines the searched range
  (lam_max = lambda0(log2(Xna), lambda = 0, maxit = 100, trace.it  = FALSE, thresh = 1e-05))  # 11052.99

  # find the optimal lambda, using a rank = 1 (r=1) softImpute
  lambda.sv=NULL
  for(i in 1:10){
    print(i); flush.console()
    re.optimize=optimize(mseImputation,interval = c(0,lam_max),tol = 1e-8,
                       r=1, sample.col=samples, pool.col=pools, Cna=log2(Xna), thresh=1e-5, maxit=100)
    lambda.sv=c(lambda.sv,re.optimize$minimum)
  }
  mean(lambda.sv)
  lambda=round(mean(lambda.sv),2) # 0.35

  # perform softimpute given lambda and r=1
  re.softimpute=softImpute(log2(Xna), rank.max = 1, lambda=lambda, type="als", trace.it = F, thresh = 1e-5, maxit=100, final.svd=T)
  # all imputed/smoothed values
  Xsi=complete(log2(Xna), re.softimpute)
  Xsi=2^Xsi  # put back on original scale

  ############# run vae on both samples and pools
  # set up the model
  # categorical covariates
  batches_cate <- meta.samples[,c("PLEX","DX","RACE","SEX")]
  batches_cate <- as.matrix(batches_cate)
  # continuous covariates
  batches_cont <- meta.samples[,c("AGEyr","AGEyr2","PMIhr")]
  batches_cont <- as.matrix(batches_cont)

  # create the input matrices
  Yinit=as.matrix(log2(Xsi))
  Yna=as.matrix(log2(Xna))

  re.vae<-vae(Xna=Yna, batches_cate=batches_cate, batches_cont=batches_cont, X_init=Yinit, 
              num_epoch=20L, print_every_epoch=1L, X=NULL,
              beta_kl=1., beta_unobs=0.5, p_feat=.2)

  Xvae.sp=2^re.vae$X_imp; rownames(Xvae.sp)=rownames(Xcomplete); colnames(Xvae.sp)=colnames(Xcomplete)
  Xvae.sp=Xvae.sp[,samples]
  
  ############# run vae on samples only
  # set up the model
  # categorical covariates
  batches_cate <- meta.samples[samples,c("PLEX","DX","RACE","SEX")]
  batches_cate <- as.matrix(batches_cate)
  # continuous covariates
  batches_cont <- meta.samples[samples,c("AGEyr","AGEyr2","PMIhr")]
  batches_cont <- as.matrix(batches_cont)
  
  # create the input matrices
  Yinit=as.matrix(log2(Xsi[,samples]))
  Yna=as.matrix(log2(Xna[,samples]))
  
  re.vae<-vae(Xna=Yna, batches_cate=batches_cate, batches_cont=batches_cont, X_init=Yinit, 
              num_epoch=20L, print_every_epoch=1L, X=NULL,
              beta_kl=1., beta_unobs=0.5, p_feat=.2)
  
  Xvae.s=2^re.vae$X_imp; rownames(Xvae.s)=rownames(Xcomplete); colnames(Xvae.s)=colnames(Xcomplete[,samples])
  Xvae.s=Xvae.s[,samples]
  
  ### determine the MSE
  Xna=Xna[,samples]
  Xcomplete=Xcomplete[,samples]
  meanSE.imp.SP=mean((log2(Xvae.sp[is.na(Xna)])-log2(Xcomplete[is.na(Xna)]))^2)  # 1.7828
  meanSE.imp.S=mean((log2(Xvae.s[is.na(Xna)])-log2(Xcomplete[is.na(Xna)]))^2)  # 1.7828
  
  meanSE.smth.SP=mean((log2(Xvae.sp[!is.na(Xna)])-log2(Xcomplete[!is.na(Xna)]))^2)  # 1.7828
  meanSE.smth.S=mean((log2(Xvae.s[!is.na(Xna)])-log2(Xcomplete[!is.na(Xna)]))^2)  # 1.7828
 
  medianSE.imp.SP=median((log2(Xvae.sp[is.na(Xna)])-log2(Xcomplete[is.na(Xna)]))^2)  # 1.7828
  medianSE.imp.S=median((log2(Xvae.s[is.na(Xna)])-log2(Xcomplete[is.na(Xna)]))^2)  # 1.7828
  
  medianSE.smth.SP=median((log2(Xvae.sp[!is.na(Xna)])-log2(Xcomplete[!is.na(Xna)]))^2)  # 1.7828
  medianSE.smth.S=median((log2(Xvae.s[!is.na(Xna)])-log2(Xcomplete[!is.na(Xna)]))^2)  # 1.7828
  
  return(data.frame(meanSE.imp.SP,meanSE.imp.S,meanSE.smth.SP,meanSE.smth.S,
                    medianSE.imp.SP,medianSE.imp.S,medianSE.smth.SP,medianSE.smth.S))
}

seeds=sample(5e5,10)
re=lapply(seeds,imputationFunction,qc.pep=qc.pep,meta.samples = meta.samples,samples=samples,pools=pools)

#meanSE.imp.SP meanSE.imp.S meanSE.smth.SP meanSE.smth.S medianSE.imp.SP medianSE.imp.S medianSE.smth.SP medianSE.smth.S
#1:      1.799776     1.739035    0.004908383  0.0013675874       0.6168122      0.5903945     0.0001591407    2.735184e-04
#2:      1.841767     1.761281    0.004445963  0.0012859904       0.6245714      0.5748666     0.0001393146    2.461982e-04
#3:      1.824643     1.811379    0.005134660  0.0008427234       0.6043966      0.5905693     0.0002012263    1.318215e-04
#4:      1.804660     1.770573    0.004833635  0.0008457552       0.6189333      0.5961791     0.0001361225    7.181267e-05
#5:      1.773094     1.724937    0.005120990  0.0015488794       0.6006656      0.5678326     0.0001788250    3.537149e-04
#6:      1.775385     1.701819    0.004664478  0.0017646212       0.5967683      0.5570720     0.0001809394    4.675940e-04
#7:      1.806951     1.745841    0.004967394  0.0007980583       0.6041405      0.5736315     0.0001422103    7.698094e-05
#8:      1.765454     1.708668    0.004548055  0.0008853806       0.5967869      0.5693675     0.0001526394    8.480849e-05
#9:      1.785440     1.729216    0.005418587  0.0013162254       0.6269714      0.5818873     0.0002851707    2.769097e-04
#10:     1.815106     1.743219    0.004991982  0.0010245658       0.6173474      0.5808113     0.0001548503    9.811163e-05

> colMeans(re)
meanSE.imp.SP     meanSE.imp.S   meanSE.smth.SP    meanSE.smth.S  medianSE.imp.SP   medianSE.imp.S medianSE.smth.SP  medianSE.smth.S 
1.7992274961     1.7435968097     0.0049034126     0.0011679787     0.6107393502     0.5782611544     0.0001730439     0.0002081470 
> colSds(as.matrix(re))
[1] 2.444806e-02 3.193976e-02 2.935757e-04 3.372609e-04 1.143878e-02 1.206097e-02 4.454044e-05 1.365812e-04

> t.test(re$meanSE.imp.SP,re$meanSE.imp.S,paired = T)

Paired t-test

data:  re$meanSE.imp.SP and re$meanSE.imp.S
t = 8.8248, df = 9, p-value = 1.002e-05
alternative hypothesis: true mean difference is not equal to 0
95 percent confidence interval:
  0.04137034 0.06989104
sample estimates:
  mean difference 
0.05563069 

> t.test(re$meanSE.smth.SP,re$meanSE.smth.S,paired = T)

Paired t-test

data:  re$meanSE.smth.SP and re$meanSE.smth.S
t = 26.059, df = 9, p-value = 8.708e-10
alternative hypothesis: true mean difference is not equal to 0
95 percent confidence interval:
  0.003411161 0.004059707
sample estimates:
  mean difference 
0.003735434 

> t.test(re$medianSE.smth.SP,re$medianSE.smth.S,paired = T)

Paired t-test

data:  re$medianSE.smth.SP and re$medianSE.smth.S
t = -0.87143, df = 9, p-value = 0.4062
alternative hypothesis: true mean difference is not equal to 0
95 percent confidence interval:
  -1.262274e-04  5.602112e-05
sample estimates:
  mean difference 
-3.510312e-05 

> t.test(re$medianSE.imp.SP,re$medianSE.imp.S,paired = T)

Paired t-test

data:  re$medianSE.imp.SP and re$medianSE.imp.S
t = 9.5867, df = 9, p-value = 5.079e-06
alternative hypothesis: true mean difference is not equal to 0
95 percent confidence interval:
  0.02481438 0.04014201
sample estimates:
  mean difference 
0.0324782 
