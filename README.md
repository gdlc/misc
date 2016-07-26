## (1) GWAS with BGData


**Loading the data**
```R
 load('~/PIC/example_sire_wg.RData')
 PC=svd(X,nu=3,nv=0)$u
 DATA=BGData(geno=X, pheno=data.frame(y=blup,pc1=PC[,1],pc2=PC[,2]),map=data.frame())
 B=GWAS(y~PC,data=DATA,model='lm',mc.cores=4) 
 plot(-log10(B[,4],type='o',col=4,cex=.5))
```
