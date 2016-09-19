setwd("C:/dev")

library(glmnet)
library(ROCR)
######################################################################################
# simulation cases from MTLasso2G
######################################################################################
# choose one dataset from each experimental setting
K<-10
#K<-20
J<-100
#J<-500

# choose one of datasets
i<-1

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
# same lambda for LassoMP, RidgeMP and elasticMP
################################################################################
l <- max(maxlambda1(X, Y))

ngrid <- 25

lambda <- 2^seq(-0.01, -10, length=ngrid) * l

nfolds <- 5


###############################################################################
# cross validation
###############################################################################

AUC_matrix<- matrix(0,nrow=11,ncol=3)

MSE_matrix<- matrix(0,nrow=11,ncol=3)

DF_matrix<- matrix(0,nrow=11,ncol=3)

MSEm_matrix<- matrix(0,nrow=11,ncol=3)

R2_matrix<- matrix(0,nrow=11,ncol=3)

aR2_matrix<- matrix(0,nrow=11,ncol=3)

cvm_matrix<- matrix(0,nrow=11,ncol=3)

p.fac = rep(1, ncol(X))

set.seed(123)
index<-sample(x=which(rowMeans(B)>0),size=1)

for ( i in 1:11) {
  print(i)
  # set the pentalty factor 
  p.fac[index] = (i-1)*0.1
  
  print(p.fac[index])

  # cross validation -- LassoMP
  cvLassoMP = cv.glmnet(x=X, y=Y, family = "mgaussian", standardize.response=FALSE, nfolds=5,lambda=lambda, penalty.factor =p.fac)

  # cross validation -- RidgeMP
  cvRidgeMP = cv.glmnet(x=X, y=Y, family = "mgaussian", standardize.response=FALSE, nfolds=5,lambda=lambda, penalty.factor =p.fac, alpha=0)

  # cross validation -- elasticMP
  cvElasticMP = cv.glmnet(x=X, y=Y, family = "mgaussian", standardize.response=FALSE, nfolds=5,lambda=lambda, penalty.factor =p.fac, alpha=0.5)


###############################################################################
# LassoMP, RidgeMP, elasticMP using the optimal lambda
###############################################################################


  lassoMP<- glmnet(x=X, y=Y, family = "mgaussian", standardize.response=FALSE,lambda=cvLassoMP$lambda.min, penalty.factor=p.fac)
  BlassoMP<-sapply(lassoMP$beta,function(x) as.matrix(x))

  RidgeMP<- glmnet(x=X, y=Y, family = "mgaussian", standardize.response=FALSE,lambda=cvRidgeMP$lambda.min, penalty.factor=p.fac, alpha=0)
  BRidgeMP<-sapply(RidgeMP$beta,function(x) as.matrix(x))

  elasticMP<- glmnet(x=X, y=Y, family = "mgaussian", standardize.response=FALSE,lambda=cvElasticMP$lambda.min, penalty.factor=p.fac, alpha=0.5)
  BelasticMP<-sapply(elasticMP$beta,function(x) as.matrix(x))

#######################################################################
# AUC
#######################################################################

  resAUC <- lapply(list(lassoMP=BlassoMP,RidgeMP=BRidgeMP, elasticMP=BelasticMP), function(B2) {
    performance(prediction(
      labels=as.numeric(B != 0),
      predictions=as.numeric(abs(B2))), "auc")
  })

  AUC<-lapply(resAUC, function(x) unlist(slot(x, "y.values")))

#print(AUC)
  AUC<-unlist(AUC)

#######################################################################
# RMSE
#######################################################################


  resMSE <- lapply(list(lassoMP=BlassoMP,RidgeMP=BRidgeMP, elasticMP=BelasticMP), function(B2) {
    performance(prediction(
      labels=as.numeric(B != 0),
      predictions=as.numeric(abs(B2))), "rmse")
  })

  nSample = nrow(X)

  MSE<-lapply(resMSE, function(x) unlist(slot(x, "y.values")))

  MSE<-unlist(MSE)

#######################################################################
# DF -- degree of freedom
#######################################################################

  df <- lapply(list(lassoMP=BlassoMP, RidgeMP=BRidgeMP, elasticMP=BelasticMP), function(B2) {
    m=sum(abs(rowMeans(B2))>0)
  })

  df<-unlist(df)


  AUC_matrix[i,]<- AUC

  MSE_matrix[i,]<- MSE

  DF_matrix[i,]<- df

}


##################################################################################################################
# plot AUC
##################################################################################################################

par(mfrow = c(1,1))
s<-seq(0,1,0.1)

AUC<-as.data.frame(AUC_matrix)

colnames(AUC)<-c("LassoMP","RidgeMP","elasticMP")


par(mar=c(5,6,4,2))
matplot(s, cbind(AUC$LassoMP,AUC$RidgeMP, AUC$elasticMP), pch=19, type="b", col=c(2,4,1),xlab = "Penalty factor", ylab="AUC",lwd=5,cex=2,cex.lab=2,cex.axis=2)
legend('bottomright', c("LassoMP", "RidgeMP", "elasticMP"), 
       pch = 19, col=c(2,4,1), bty='n',cex=2, pt.cex = 2)

##########################################################################################################
# plot MSE
##################################################################################################################

MSE<-as.data.frame(MSE_matrix)

colnames(MSE)<-c("LassoMP","RidgeMP","elasticMP")

par(mar=c(5,6,4,2))
matplot(s, cbind(MSE$LassoMP,MSE$RidgeMP, MSE$elasticMP), pch=19, type="b", col=c(2,4,1),xlab = "Penalty factor", ylab="RMSE",lwd=5,cex=2,cex.lab=2,cex.axis=2)
legend('bottomright', c("LassoMP", "RidgeMP", "elasticMP"), 
       pch = 19, col=c(2,4,1), bty='n',cex=2, pt.cex = 2)

##########################################################################################################
# plot DF
##################################################################################################################

DF<-as.data.frame(DF_matrix)

colnames(DF)<-c("LassoMP","RidgeMP","elasticMP")

par(mar=c(5,6,4,2))
matplot(s, cbind(DF$LassoMP,DF$RidgeMP, DF$elasticMP), pch=19, type="b", col=c(2,4,1),xlab = "Penalty factor", ylab="DF",lwd=5,cex=2,cex.lab=2,cex.axis=2)
legend('right', c("LassoMP", "RidgeMP", "elasticMP"), 
       pch = 19, col=c(2,4,1), bty='n',cex=2, pt.cex = 2)

