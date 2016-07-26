## (1) GWAS with BGData


** (1.1) Loading the data**
```R
 rm(list=ls())
 library(BGData)
 load('~/PIC/example_sire_wg.RData')
```

** (1.2) Creating a BGData object***
```R
 DATA=BGData(geno=X, pheno=data.frame(y=blup,pc1=PC[,1],pc2=PC[,2]),map=data.frame())
```
**(1.3) Running GWAS using BGData**
```R
 B=GWAS(y~pc1+pc2,data=DATA,method='lm',mc.cores=4,verbose=T) 
 plot(-log10(B[,4]),type='o',col=4,cex=.5)
```


## (2) BGLR

** (2.1) Gaussian Prior (`BRR`)
 
 ```R
  fm=BGLR(y=blup,ETA=list(list(X=X,model='BRR')),nIter=55000,burnIn=5000)
 ```
 
 ** (2.1) Gaussian Prior (`BRR`)
 
 ```R
  fm=BGLR(y=blup,ETA=list(list(X=X,model='BRR')),nIter=55000,burnIn=5000)
 ```
 
 ** (2.2) BayesA (`BayesA`)
 
 ```R
  fm=BGLR(y=blup,ETA=list(list(X=X,model='BRR')),nIter=55000,burnIn=5000)
 ```
 

 ** (2.3) BayesB (`BayesB`)
 
 ```R
  fm=BGLR(y=blup,ETA=list(list(X=X,model='BRR')),nIter=55000,burnIn=5000)
 ```
 
 ** (2.4) BayesC (`BayesC`)
 
 ```R
  fm=BGLR(y=blup,ETA=list(list(X=X,model='BRR')),nIter=55000,burnIn=5000)
 ```
 
 ** (2.5) Bayesian-LASSO (`BL`)
 
 ```R
  fm=BGLR(y=blup,ETA=list(list(X=X,model='BRR')),nIter=55000,burnIn=5000)
 ```
 
