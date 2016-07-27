## (1) GWAS with BGData


**(1.1) Loading the data**
```R
 rm(list=ls())
 library(BGData)
 setwd('../output')
 load('../data/example_sire_wg.RData')
```

**(1.2) Creating a BGData object***
```R
 DATA=BGData(geno=X, pheno=data.frame(y=blup,pc1=PC[,1],pc2=PC[,2]),map=data.frame())
```
**(1.3) Running GWAS using BGData**
```R
 B<-GWAS(y~pc1+pc2,data=DATA,method='lm',mc.cores=10,verbose=T) # 6.9 seq with 10 cores
 plot(-log10(B[,4]),type='o',col=4,cex=.5)
```

## (2) BGLR

**(2.1) Gaussian Prior (`BRR`)**
 ```R
  nIter=55000
  burnIn=5000
  X=X/sqrt(ncol(X))
  fmBRR=BGLR(y=blup,ETA=list(list(X=X,model='BRR')),nIter=nIter,burnIn=burnIn,saveAt='BRR_') # ~.05 sec/iteration
  save(fmBRR,file='fmBRR.RData')
 ```
 
 **(2.1) Scaled-t prior (`BayesA`)**
 
 ```R
  fmBA=BGLR(y=blup,ETA=list(list(X=X,model='BayesA')),nIter=nIter,burnIn=burnIn,saveAt='BA_') 
  save(fmBA,file='fmBA.RData')
 ```
 
 **(2.2) Spike-Slab-1 (`BayesB`)**
 
 ```R
  fmBB=BGLR(y=blup,ETA=list(list(X=X,model='BayesB')),nIter=nIter,burnIn=burnIn,saveAt='BB_') 
  save(fmBB,file='fmBB.RData')
 ```
 

 **(2.3) BayesB (`BayesB`)**
 
 ```R
  fm=BGLR(y=blup,ETA=list(list(X=X,model='BRR')),nIter=55000,burnIn=5000)
 ```
 
 **(2.4) BayesC (`BayesC`)**
 
 ```R
  fm=BGLR(y=blup,ETA=list(list(X=X,model='BRR')),nIter=55000,burnIn=5000)
 ```
 
 **(2.5) Bayesian-LASSO (`BL`)**
 
 ```R
  fm=BGLR(y=blup,ETA=list(list(X=X,model='BRR')),nIter=55000,burnIn=5000)
 ```
 
## (3) 5-fold Cross-validation

```R
N=nrow(X)
folds=rep(1:5,each=ceiling(N/5))[1:N]
folds=folds[order(runif(N))]

yHatCV_BRR=rep(NA,N)
yHatCV_BB=rep(NA,N)

for(i in 1:5){
 yNA=blup
 tst=which(folds==i)
 yNA[tst]=NA
 fm=BGLR(y=yNA,ETA=list(list(X=X,model='BRR')),nIter=3500,burnIn=500,verbose=F)
  yHatCV_BRR[tst]=fm$yHat[tst]
 fm=BGLR(y=yNA,ETA=list(list(X=X,model='BayesB')),nIter=3500,burnIn=500,verbose=F)
  yHatCV_BB[tst]=fm$yHat[tst]
  print(i)
}

```
