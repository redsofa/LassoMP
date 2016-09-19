setwd("C:/dev/FMPR")

library(glmnet)
library(Matrix)

##################################################################################################
# cis eQTLs
##################################################################################################

cis_eQTL <- function(corMatrix, snpLabel,probesPos, cisThreshold=15){
  
  corMatrix<-sapply(corMatrix, function(x) as.matrix(x))
  
  corMatrix1<-Matrix(corMatrix)
  
  corMatrix2<-as.matrix(summary(corMatrix1))
  
  
  # get eQTLs from index Matrix 
  eQTLs <- data.frame(Probe=rownames(ge)[corMatrix2[,2]],SNP=colnames(snp)[corMatrix2[,1]], Coef=corMatrix2[,3])
  
  # merge eQTLs with SNP and Probe position information
  merged <- Reduce(function(x, y) merge(x, y), list(eQTLs, probePos, snpLabel))
  
  eQTLsAll<-merged[complete.cases(merged),]
  
  # threshold for cis eQTLs
  cisThreshold <- 15
  
  cis<-with(eQTLsAll, (ProbeChromosome==SNPchr) & (abs((as.numeric(ProbecM)-as.numeric(SNPpos)))<cisThreshold))
  
  
  eQTLsAll<-cbind(eQTLsAll,cis)
  
  # pcercent of cis eQTLs
  perciseQTL<-sum(cis+0)/nrow(eQTLsAll)
  
  return(list(eQTLsAll,perciseQTL))
}

#############################################################################
# # read SNP and gene expression data
#############################################################################

snpFile = "./data/StMx/BarleySNP.txt"

snp <- read.table(snpFile, header=TRUE, sep="\t", row.names=1,as.is=TRUE, check.names=FALSE)

snp<-as.matrix(snp)

# read gene expression  data as Y

geneExpressionFile= "./data/StMx/plum1.txt"  

#geneExpressionFile= "./data/StMx/skyblue.txt"  

#geneExpressionFile= "./data/StMx/saddlebrown.txt"  


ge <- read.table(geneExpressionFile, header=TRUE, sep="\t", row.names=1)

# read snp label data

snpLabelFile = "./data/StMx/snpAll.txt"

snpLabel <- read.table(snpLabelFile, header=TRUE, sep="\t")

snpLabel<-snpLabel[,1:3]

# read probe position data

probePosFile = "./data/StMx/ProbePos.txt"

probePos <- read.table(probePosFile, header=TRUE, sep="\t")


probePos<-subset(probePos,probePos$ProbeChromosome != 'unknown')

#############################################################################
# normliazation
#############################################################################

X<-snp

# glmnet handle normalization standardize.response=TRUE
Y<-t(ge)

##################################################################
# cross validation on LassoM
##################################################################

# LassoM  
set.seed(123)
cvfit = cv.glmnet(x=X, y=Y, family = "mgaussian", standardize.response=TRUE, nfolds=5)


#############################################################################
# LassoM -- MSE, DF, %DEV, and Percent of cis eQTLs
#############################################################################
#LassoM is the solution $lambda.min

LassoM<- glmnet(x=X, y=Y, family = "mgaussian", standardize.response=TRUE,lambda=cvfit$lambda.min)


lambdaMinIndex = which(cvfit$glmnet.fit$lambda==cvfit$lambda.min)

# Mean Sqaured Error
LassoM_MSE<-cvfit$cvm[lambdaMinIndex]

# Number of predictors
LassoM_DF<-LassoM$df

#Proportion of Error explained
LassoM_PDE<-LassoM$dev.ratio

LassoM_cis<-cis_eQTL(LassoM$beta,snpLabel,probesPos)
# Percent of cis eQTLs
print(LassoM_cis[[2]])

#############################################################################
# LassoMs -- MSE, DF, %DEV, and Percent of cis eQTLs
#############################################################################
#LassoM is the solution $lambda.1se

LassoMs<- glmnet(x=X, y=Y, family = "mgaussian", standardize.response=TRUE,lambda=cvfit$lambda.1se)

lambda1seIndex = which(cvfit$glmnet.fit$lambda==cvfit$lambda.1se)

# Mean Sqaured Error
LassoMs_MSE<-cvfit$cvm[lambda1seIndex]

# Number of predictors
LassoMs_DF<-LassoMs$df

#Proportion of Error explained
LassoMs_PDE<-LassoMs$dev.ratio

LassoMs_cis<-cis_eQTL(LassoMs$beta,snpLabel,probesPos)
# Percent of cis eQTLs
print(LassoMs_cis[[2]])

##################################################################
# cross validation on LassoMP
##################################################################


#LassoMP
p.fac = rep(1, 413)
# prior knowledge on plum1 and skyblue
p.fac[c(353)] = 0
# prior knowledge on saddlebrown
#p.fac[c(367,370:371)] = 0

set.seed(123)
cvfitK = cv.glmnet(x=X, y=Y, family = "mgaussian", standardize.response=TRUE, nfolds=5, penalty.factor=p.fac)

#############################################################################
# LassoMP -- MSE, DF, %DEV, and Percent of cis eQTLs
#############################################################################
#LassoMP is the solution $lambda.min + Penalty factor
LassoMP<- glmnet(x=X, y=Y, family = "mgaussian", standardize.response=TRUE,lambda=cvfitK$lambda.min, penalty.factor=p.fac)

lambdaMinIndex = which(cvfitK$glmnet.fit$lambda==cvfitK$lambda.min)

# Mean Sqaured Error
LassoMP_MSE<-cvfitK$cvm[lambdaMinIndex]


# Number of predictors
LassoMP_DF<-LassoMP$df

#Proportion of Error explained
LassoMP_PDE<-LassoMP$dev.ratio

LassoMP_cis<-cis_eQTL(LassoMP$beta,snpLabel,probesPos)
# Percent of cis eQTLs
print(LassoMP_cis[[2]])

#############################################################################
# LassoMPs -- MSE, DF, %DEV, and Percent of cis eQTLs
#############################################################################
#LassoMPs is the solution $lambda.1se + Penalty factor

LassoMPs<- glmnet(x=X, y=Y, family = "mgaussian", standardize.response=TRUE,lambda=cvfitK$lambda.1se, penalty.factor=p.fac)


lambda1seIndex = which(cvfitK$glmnet.fit$lambda==cvfitK$lambda.1se)

# Mean Sqaured Error
LassoMPs_MSE<-cvfitK$cvm[lambda1seIndex]

# Number of predictors
LassoMPs_DF<-LassoMPs$df

#Proportion of Error explained
LassoMPs_PDE<-LassoMPs$dev.ratio


LassoMPs_cis<-cis_eQTL(LassoMPs$beta,snpLabel,probesPos)

# Percent of cis eQTLs
print(LassoMPs_cis[[2]])

