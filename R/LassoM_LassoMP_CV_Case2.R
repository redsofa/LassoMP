
setwd("C:/dev/FMPR")

library(glmnet)

#############################################################################
# # read genotype data 
#############################################################################

snp <- read.csv(file= "./data/QSM/Geno_QSM_Impute.csv", header=TRUE, stringsAsFactors = FALSE)
#head(snp)[1:5,1:5]

snpData<-snp[,10:ncol(snp)]

snpData[snpData=="A"]<-0
snpData[snpData=="B"]<-2


snpData<-as.matrix(snpData)

snpDataNew<-matrix(as.numeric(snpData),nrow=nrow(snpData),ncol=ncol(snpData))

snpDataNew<-t(snpDataNew)

###############################################################################
# read gene expression data
###############################################################################
# read gene expression  as Y

geneExpressionFile = "./data/QSM/saddlebrownINCO-MOCK.csv"

#geneExpressionFile = "./data/QSM/darkgreyINCO-MOCK.csv"

#geneExpressionFile = "./data/QSM/darkgreenINCO-MOCK.csv"

#geneExpressionFile = "./data/QSM/yellowINCO-MOCK.csv"


ge <- read.csv(geneExpressionFile, header=TRUE, row.names=1)
ge<-t(ge)


###################################################################
###################################################################

X <- snpDataNew

Y <- ge

##############################################################################
# compare LassoM, LassoMP on normalized Y
##############################################################################

# lasso Multiple Gaussian model


set.seed(123)
lasso.err = rep(0,30)

for (i in 1:30) {
  print(i)
  
  cvfit = cv.glmnet(x=X, y=Y, family = "mgaussian", standardize.response=TRUE, nfolds=5)
  lasso.err[i] <- sort(cvfit$cvm)[1]
}


##################################################################################
# lasso Multiple Gaussian model with prior knowledge
##################################################################################
set.seed(123)


# saddlebrown
p.fac = rep(1, 378)
p.fac[330] = 0

# darkgrey
#p.fac = rep(1, 378)
#p.fac[255] = 0

# darkgreen
#p.fac = rep(1, 378)
#p.fac[351] = 0

# yellow
#p.fac = rep(1, 378)
#p.fac[99] = 0


lassoK.err = rep(0,30)

for (i in 1:30) {
  
  print(i)
  
  
  cvfitK = cv.glmnet(x=X, y=Y, family = "mgaussian", standardize.response=TRUE, nfolds=5, penalty.factor =p.fac)
  lassoK.err[i] <- sort(cvfitK$cvm)[1]
}

############################################################################################
# MSE
############################################################################################


MSE<- data.frame(LassoM=lasso.err, LassoMP = lassoK.err)

boxplot(MSE, xlab="", ylab="Mean Squared Error",notch=TRUE, col=(c("red", "blue")))

