setwd("C:/dev")

#FMRP installation
#library(devtools)
#install_github("gabraham/FMPR/FMPR")


library(FMPR)
library(glmnet)

library(ROCR)

######################################################################################
# simulation cases from MTLasso2G
######################################################################################

K<-10
#K<-20
J<-100
#J<-500

for ( i in 1:30) {
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

  ################################################################################
  # normaliztion 
  ################################################################################

  X<-scale(X)
  Y<-scale(X %*% B + rnorm(nrow(X) * ncol(Y)))


  # FMPR GFLasso
  C <- gennetwork(Y, cortype=2)

################################################################################
  l <- max(maxlambda1(X, Y))


  ngrid <- 25
  lambda <- 2^seq(-0.01, -10, length=ngrid) * l

  gamma <- c(0, 10^seq(-5, 5, length=ngrid))


  nfolds <- 5

###############################################################################
# cross validation
###############################################################################


  # FMPR

  opt.f <- optim.fmpr(X=X, Y=Y, cortype=2,
                      lambda=lambda, gamma=gamma, nfolds=nfolds)

  # GFLasso
  opt.s <- optim.spg(X=X, Y=Y, cortype=2,
                       lambda=lambda, gamma=gamma, nfolds=nfolds)

  # LassoM -- lasso Multiple Gaussian model
  cvfit = cv.glmnet(x=X, y=Y, family = "mgaussian", standardize.response=FALSE, nfolds=5,lambda=lambda)


  # RidgeM -- Ridge Multiple Gaussian model
  cvfit_Ridge = cv.glmnet(x=X, y=Y, family = "mgaussian", standardize.response=FALSE, nfolds=5,lambda=lambda, alpha=0)


  # elasticM -- elastic Multiple Gaussian model
  cvfit_elastic = cv.glmnet(x=X, y=Y, family = "mgaussian", standardize.response=FALSE, nfolds=5,lambda=lambda, alpha=0.5)



  p.fac = rep(1, ncol(X))

  index<-sample(x=which(rowMeans(B)>0),size=1)

  p.fac[index] = 0

  # LassoMP
  cvfitK = cv.glmnet(x=X, y=Y, family = "mgaussian", standardize.response=FALSE, nfolds=5,lambda=lambda, penalty.factor =p.fac)

  # RidgeMP
  cvfitK1 = cv.glmnet(x=X, y=Y, family = "mgaussian", standardize.response=FALSE, nfolds=5,lambda=lambda, penalty.factor =p.fac, alpha=0)

  # elasticMP
  cvfitK2 = cv.glmnet(x=X, y=Y, family = "mgaussian", standardize.response=FALSE, nfolds=5,lambda=lambda, penalty.factor =p.fac, alpha=0.5)


  simuFile <- paste("./data/simu/simu",i,"N50K",K,"J",J,".RData", sep ="") 

  save(X,Y,B,C, p.fac,opt.f, opt.s,cvfit,cvfit_Ridge,cvfit_elastic,cvfitK,cvfitK1,cvfitK2,file=simuFile)

}

###################################################################################################
# evaluate the performance of FMPR, GFLasso, LassoM, LassoMP, RidgeM, RidgeMP, elasticM, elasticMP
###################################################################################################


AUC_matrix<- matrix(0,nrow=30,ncol=8)

MSE_matrix<- matrix(0,nrow=30,ncol=8)

DF_matrix<- matrix(0,nrow=30,ncol=8)

for ( i in 1:30) {
  print(i)
  
  simuFile <- paste("./data/simu/simu",i,"N50K",K,"J",J,".RData", sep ="") 
  
  load(simuFile)
  ###################################################################################################
  # FMPR, GFLasso, LassoM, LassoMP, RidgeM, RidgeMP, elasticM, elasticMP using the optimal parameters
  ###################################################################################################

  #FMPR
  f <- fmpr(X=X, Y=Y, C=C, lambda=opt.f$opt["lambda"],
          gamma=opt.f$opt["gamma"], simplify=TRUE)

  #GFLasso

  s <- spg(X=X, Y=Y, C=C, lambda=opt.s$opt["lambda"],
         gamma=opt.s$opt["gamma"], simplify=TRUE)

  #LasspM
  lassoM<- glmnet(x=X, y=Y, family = "mgaussian", standardize.response=FALSE,lambda=cvfit$lambda.min)
  lm<-sapply(lassoM$beta,function(x) as.matrix(x))

  #RidgeM
  RidgeM<- glmnet(x=X, y=Y, family = "mgaussian", standardize.response=FALSE,lambda=cvfit_Ridge$lambda.min, alpha=0)
  rm<-sapply(RidgeM$beta,function(x) as.matrix(x))

  #elasticM
  elasticM<- glmnet(x=X, y=Y, family = "mgaussian", standardize.response=FALSE,lambda=cvfit_elastic$lambda.min, alpha=0.5)
  em<-sapply(elasticM$beta,function(x) as.matrix(x))

  # LassoMP
  lassoMP<- glmnet(x=X, y=Y, family = "mgaussian", standardize.response=FALSE,lambda=cvfitK$lambda.min, penalty.factor =p.fac)
  lmK<-sapply(lassoMP$beta,function(x) as.matrix(x))


  RidgeMP<- glmnet(x=X, y=Y, family = "mgaussian", standardize.response=FALSE,lambda=cvfitK1$lambda.min, penalty.factor =p.fac, alpha=0)
  lmK1<-sapply(RidgeMP$beta,function(x) as.matrix(x))

  elasticMP<- glmnet(x=X, y=Y, family = "mgaussian", standardize.response=FALSE,lambda=cvfitK2$lambda.min, penalty.factor =p.fac, alpha=0.5)
  lmK2<-sapply(elasticMP$beta,function(x) as.matrix(x))



#######################################################################
# AUC
#######################################################################



resAUC <- lapply(list(FMPR=f, GFlasso=s, RidgeM=rm, RidgeMP=lmK1, elasticM=em, elasticMP=lmK2, LassoM=lm, lassoMP=lmK), function(B2) {
    performance(prediction(
      labels=as.numeric(B != 0),
      predictions=as.numeric(abs(B2))), "auc")
})

AUC<-lapply(resAUC, function(x) unlist(slot(x, "y.values")))

AUC<-unlist(AUC)

#######################################################################
# MSE
#######################################################################


resMSE <- lapply(list(FMPR=f, GFlasso=s, RidgeM=rm, RidgeMP=lmK1, elasticM=em, elasticMP=lmK2,LassoM=lm, lassoMP=lmK), function(B2) {
  performance(prediction(
    labels=as.numeric(B != 0),
    predictions=as.numeric(abs(B2))), "rmse")
})

nSample = nrow(X)

MSE<-lapply(resMSE, function(x) unlist(slot(x, "y.values")))

MSE<-unlist(MSE)

#######################################################################
# DF
#######################################################################

df <- lapply(list(FMPR=f, GFlasso=s,  RidgeM=rm, RidgeMP=lmK1, elasticM=em, elasticMP=lmK2, LassoM=lm, lassoMP=lmK), function(B2) {
  m=sum(abs(rowMeans(B2))>0)
})

df<-unlist(df)


AUC_matrix[i,]<- AUC

MSE_matrix[i,]<- MSE

DF_matrix[i,]<- df


}

##############################################################################################################
# save the results
##############################################################################################################
simuSumFile <- paste("./data/simu/simu","N50K",K,"J",J,".RData", sep ="") 


save(AUC_matrix, MSE_matrix,DF_matrix,MSEm_matrix,R2_matrix, aR2_matrix, file = simuSumFile)

##############################################################################################################
# mean of 30 data sets
##############################################################################################################

MSE_matrix<-t(MSE_matrix)

AUC_matrix<-t(AUC_matrix)

DF_matrix<-t(DF_matrix)


RMSE<-rowMeans(MSE_matrix)

AUC<-rowMeans(AUC_matrix)

DF<-rowMeans(DF_matrix)

##############################################################################################################
# plot RMSE, AUC and DF
##############################################################################################################

dev.off()
par(mfrow = c(3,1))
par(mar=c(5,6,4,3))
barplot(RMSE,col= c(3:8,"darkblue",2), ylab="RMSE",cex.lab=1.5,cex.axis=1.5,
        names.arg=c("FMPR", "GFLasso","RidgeM", "RidgeMP", "elasticM", "elasticMP", "LassoM", "LassoMP"),cex.names=1.5)
barplot(AUC, ylim=c(0,1), col= c(3:8,"darkblue",2),  ylab="AUC",cex.lab=1.5,cex.axis=1.5,
        names.arg=c("FMPR", "GFLasso","RidgeM", "RidgeMP", "elasticM", "elasticMP", "LassoM", "LassoMP"),cex.names=1.5)
barplot(DF, col= c(3:8,"darkblue",2), ylab="DF", cex.lab=1.5,cex.axis=1.5,
        names.arg=c("FMPR", "GFLasso","RidgeM", "RidgeMP", "elasticM", "elasticMP", "LassoM", "LassoMP"),cex.names=1.5)





