#-------------------------------------------------------------------
# Update Date: 06/06/2020
# Create Date: 06/06/2020
# 
# Author: Jin Jin
# Email: jjin31@jhu.edu
#-------------------------------------------------------------------
## directly run on cluster
rm(list=ls())
temp <- commandArgs(TRUE)
chr = as.numeric(temp[1])
plinkout_path = '/dcl01/chatterj/data/jin/UKB/GWAS/egfr/'
setwd(plinkout_path)


# PLINK code for running simple linear regression:
slmcode=paste(
  paste0('/dcl01/chatterj/data/tools/plink2'),
  paste0('--threads 1'),
  paste0('--bfile chr',chr),
  paste0('--snps-only'),
  paste0('--keep gwasid.txt'),
  paste0('--pheno pheno_adj.txt'),
  paste0('--glm'),
  paste0('--out ',plinkout_path,'GWAS/chr',chr)
)
system(slmcode)

#--glm hide-covar --vif 10000 --covar-variance-standardize \
#--out /dcl01/chatterj/data/jzhang2/pwas/pipeline/Results/White/pQTL/gwas_result/allchr/${n}/chr${n}





## validate:
temp = read_plink('/dcl01/chatterj/data/jin/UKB/GWAS/egfr/chr22')
names(temp)=c('bed','bim','fam')
#dim(temp$bed)
#[1]  17435 289432

X = temp$bed['rs7287144',]
nomissing.id = names(which(!is.na(X)))
X = X[nomissing.id]
Y = read.table('/dcl01/chatterj/data/jin/UKB/GWAS/egfr/pheno_adj.txt',header=T)
rownames(Y) = Y$IID
Y = Y[nomissing.id,'pheno']

fit1 = lm(Y~X)