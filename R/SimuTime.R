setwd("C:/dev/FMPR")

library(FMPR)
library(glmnet)
######################################################################################
# simulation cases from MTLasso2G
######################################################################################
nDataSet<-30

FMPR_time <- rep(0,nDataSet)
GFLasso_time <- rep(0,nDataSet)
LassoM_time <- rep(0,nDataSet)
RidgeM_time <- rep(0,nDataSet)
elasticM_time <- rep(0,nDataSet)

LassoMP_time <- rep(0,nDataSet)
RidgeMP_time <- rep(0,nDataSet)
elasticMP_time <- rep(0,nDataSet)

nfolds <- 5


K<-10
#K<-20
J<-100
#J<-500

for (i in 1:30) {
  
  print(i)
  
  xFile <- paste("./data/simu/X",i,"N50K",K,"J",J,".txt", sep ="")  
  yFile <- paste("./data/simu/Y",i,"N50K",K,"J",J,".txt", sep ="")  
  bFile <- paste("./data/simu/B",i,"N50K",K,"J",J,".txt", sep ="") 


  X<-read.table(file=xFile,sep=",")

  Y<-read.table(file=yFile,sep=",")

  B<-read.table(file=bFile,sep=",")


  X<-as.matrix(X)
  Y<-as.matrix(Y)
  B<-as.matrix(B)


  X<-scale(X)
  Y<-scale(X %*% B + rnorm(nrow(X) * ncol(Y)))

################################################################################
# same lambda and gamma for all algorithms
################################################################################
  l <- max(maxlambda1(X, Y))


  ngrid <- 25
  lambda <- 2^seq(-0.01, -10, length=ngrid) * l
  gamma <- c(0, 10^seq(-5, 5, length=ngrid))

###############################################################################
# cross validation
###############################################################################

  
  # FMPR
  FMPR_time[i]<-system.time({
    opt.f <- optim.fmpr(X=X, Y=Y, cortype=2,
                        lambda=lambda, gamma=gamma, nfolds=nfolds)
  })
  
  
  # GFLasso
  GFLasso_time[i]<-system.time({
    opt.s <- optim.spg(X=X, Y=Y, cortype=2,
                       lambda=lambda, gamma=gamma, nfolds=nfolds)
  })
  
  # LassoM - lasso Multiple Gaussian model
  LassoM_time[i]<-system.time({
  cvfit = cv.glmnet(x=X, y=Y, family = "mgaussian", standardize.response=FALSE, nfolds=5,lambda=lambda)
  
  })

  # Ridge
  RidgeM_time[i]<-system.time({
  cvfit_Ridge = cv.glmnet(x=X, y=Y, family = "mgaussian", standardize.response=FALSE, nfolds=5,lambda=lambda, alpha=0)
  })  
  
  # elastic 
  elasticM_time[i]<-system.time({
  cvfit_elastic = cv.glmnet(x=X, y=Y, family = "mgaussian", standardize.response=FALSE, nfolds=5,lambda=lambda, alpha=0.5)
  })
  
  
  p.fac = rep(1, ncol(X))
  
  index<-sample(x=which(rowMeans(B)>0),size=1)
  
  p.fac[index] = 0
  
  LassoMP_time[i]<-system.time({  
  cvfitK = cv.glmnet(x=X, y=Y, family = "mgaussian", standardize.response=FALSE, nfolds=5,lambda=lambda, penalty.factor=p.fac)
  })  

  RidgeMP_time[i]<-system.time({  
  cvfitK1 = cv.glmnet(x=X, y=Y, family = "mgaussian", standardize.response=FALSE, nfolds=5,lambda=lambda, penalty.factor=p.fac, alpha=0)
  })

  elasticMP_time[i]<-system.time({
  cvfitK2 = cv.glmnet(x=X, y=Y, family = "mgaussian", standardize.response=FALSE, nfolds=5,lambda=lambda, penalty.factor=p.fac, alpha=0.5)
  })

}

print(mean(FMPR_time))

print(mean(GFLasso_time))

print(mean(RidgeM_time))

print(mean(RidgeMP_time))

print(mean(elasticM_time))

print(mean(elasticMP_time))

print(mean(LassoM_time))

print(mean(LassoMP_time))


