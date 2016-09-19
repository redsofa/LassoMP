setwd("C:/dev/FMPR")

library(glmnet)
#############################################################################
# # read genotype data , phenotype  and gene expression data
#############################################################################
# read SNP data  as X

snpFile = "./data/StMx/BarleySNP.txt"

snp <- read.table(snpFile, header=TRUE, sep="\t", row.names=1,as.is=TRUE, check.names=FALSE)

snp<-as.matrix(snp)

# read gene expression  data as Y

geneExpressionFile= "./data/StMx/plum1.txt"  

geneExpressionFile= "./data/StMx/skyblue.txt"  

geneExpressionFile= "./data/StMx/saddlebrown.txt"  


ge <- read.table(geneExpressionFile, header=TRUE, sep="\t", row.names=1)


ge<-t(ge)
###################################################################
# normalization
###################################################################

X<-snp

# glmnet handle normalization standardize.response=TRUE
Y<-ge

##############################################################################
# LassoM
# lasso Multiple Gaussian model
##############################################################################

# lasso Multiple Gaussian model
set.seed(123)
lasso.err = rep(0,30)

for (i in 1:30) {
  print(i)

  cvfit = cv.glmnet(x=X, y=Y, family = "mgaussian", standardize.response=TRUE, nfolds=5)
  lasso.err[i] <- sort(cvfit$cvm)[1]
}

lasso.err

##################################################################################
# LassoMP
# lasso Multiple Gaussian model with prior knowledge
##################################################################################
set.seed(123)
p.fac = rep(1, 413)
# plum1 
p.fac[c(353)] = 0

#skyblue
#p.fac[353] = 0

# saddlebrown
#p.fac[c(367,370:371)] = 0

lassoK.err = rep(0,30)

for (i in 1:30) {

  print(i)
  
  cvfitK = cv.glmnet(x=X, y=Y, family = "mgaussian", standardize.response=TRUE, nfolds=5,penalty.factor=p.fac)
  lassoK.err[i] <- sort(cvfitK$cvm)[1]
}

############################################################################################
# MSE plot from glmnet lasso cvm
############################################################################################

MSE<- data.frame(LassoM=lasso.err, LassoMP = lassoK.err)


boxplot(MSE, xlab="", ylab="Mean Squared Error",notch=TRUE, col=(c("red", "blue")))


