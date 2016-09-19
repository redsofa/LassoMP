# LassoMP

Software Description
====================


Source code used for simulation and case studies are described in the followings. 

1. Simulation

  The simulation data was generated using the Matlab code in MMT (http://www.ece.ubc.ca/~xiaohuic/code/multi-task_lasso/multi-task_lasso.htm). Thirty datasets in each of four experimental settings: the number of predictors X =100 and X=500, and the number of responses Y=10 and Y-20 were generated using simuData in MMT. The number of samples is set as 50. In total, 120 datasets(X, Y, and the covariate matrix B) used in the simulation can be found in ./data/simu. 
  
  SimuLoop - compare eight multivariate models
  
  glmnet(https://cran.r-project.org/web/packages/glmnet/index.html) and FMPR R packages(https://github.com/gabraham/FMPR) were used in simulation. LassoM, LassoMP, RidgeM, RidgeMP, elasticM, and elasticMP are implemented in glmnet. FMPR and GFLasso are implemented in FMPR.  The comparison of these eight algorithms is provided in SimuLoop.R. Figure 1 is also generated.
  
  
  SimuTime - compare the computational time of eight multivariate models
  
  SimuTime was used for comparing the average computational time of eight algorithms (FMPR, GFLasso, RidgeM, RidgeMP, elasticM, elasticMP, LassoM and LassoMP) in four experimental settings.
  
  
  
  SimuPenalty is used for comparing the influence of the penalty factor on LassoMP,  RidgeMP, and elasticMP. 



2.	Case studies

  Case study 1: Genotype and gene expression data were downloaded from http://www.genenetwork.org/. Genotype data and three gene modules data are available from ./data/StMx folder. 
  
  Case study 2: The genotype and phenotype data were downloaded from Moscou 2011, and the gene expression data were downloaded from GEO http://www.ncbi.nlm.nih.gov/geo/ (GSE20416). Genotype data and four gene modules are available from ./data/QSM folder. 
  
  
  
  LassoM_LassoMP_CV_Case1.R and LassoM_LassoMP_CV_Case2.R are used to compare LassoM and LassoMP in case study 1 and 2, and Figure 2 and 3 are generated respectively. 
  
  
  
  LassoM_LassoMPs.R is used to compare LassoM, LassoMs, LassoMP, LassoMPs in case study 1. Table 4 in the main text and additional files are generated.  



References:
===========
Moscou, M.J., Lauter, N., Ste_enson, B., Wise, R.P.: Quantitative and qualitative stem rust resistance factors in barley are associated with transcriptional suppression of defense regulons. PLoS Genet 7(7), 1-17 (2011) (http://journals.plos.org/plosgenetics/article?id=10.1371%2Fjournal.pgen.1002208)
