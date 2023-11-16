#----------------------
# Created 02/27/2020
# Updated 02/27/2020
# calculate the PRS by summing the 22 PRSs from 22 chromosomes
#----------------------
rm(list=ls())
library(readr)
library(data.table)
library(mvtnorm)
library(devtools)
library(lavaan)
library(gdata)
library(xtable)
library(MASS) # for the ginv
library(data.table)
library(corpcor) #for pseudoinverse
library(parallel)
library(MendelianRandomization) # for mr_ivw
library(dplyr)
library(R.utils) # for gzip
library(stringr) # for str_detect
library(genio) # a package to facilitate reading and writing genetics data. The focus of this vignette is processing plink BED/BIM/FAM files.
library(data.table)
library(pROC)

traitlist= c("eGFRcr-LDPRED")
maindirec = "/dcl01/chatterj/data/yzhang/"
jindir = "/dcl01/chatterj/data/jin/"

temp <- commandArgs(TRUE)
trait_name =  traitlist[as.numeric(temp[1])]
trait_name = 'eGFRcr-LDPRED'
trait = 'eGFRcr'
traitcolname = 'egfr'
prsmethod = 'p+t'
refdat ='1K'  #'5000UKB'

population = 'AA'
#K=5 # the number of genetic PCs, not necessary for eGFRcr
#---------------------------------------#---------------------------------------
# input variables !!!!!!!!!!!!!!!!!!!
#---------------------------------------#---------------------------------------
for (population in c('EA', 'AA')){
  validate.dir = paste0('/dcl01/chatterj/data/zyu/')
  output_path = paste0('/dcs04/nilanjan/data/jjin/UKB/transethnic/aric/',population,'/')
  writescore_path = paste0('/dcs04/nilanjan/data/jjin/UKB/transethnic/aric-prs/',population,'/')
  if (!dir.exists(writescore_path)){dir.create(writescore_path)}
  
  #if ((prsmethod == 'ldpred')&(refdat == '1K'))
  
  fractions = 0.3 #c(1,0.3,0.1,0.03,0.01,0.003,0.001)
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # freeze below
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  for(pthr in fractions){
    if(pthr==1) pthrname = "_p1.0000e+00"
    if(pthr==0) pthrname = "-inf"
    if(pthr==0.1) pthrname = "_p1.0000e-01"
    if(pthr==0.3) pthrname = "_p3.0000e-01"
    if(pthr==0.01) pthrname = "_p1.0000e-02"
    if(pthr==0.03) pthrname = "_p3.0000e-02"
    if(pthr==4e-3) pthrname = "_p4.0000e-03"
    if(pthr==1e-3) pthrname = "_p1.0000e-03"
    if(pthr==3e-3) pthrname = "_p3.0000e-03"
    if(pthr==1e-4) pthrname = "_p1.0000e-04"
    if(pthr==3e-4) pthrname = "_p3.0000e-04"
    #---------------------------------------#---------------------------------------
    y = 0;
    for(chrnum in c(1:22)){
      temfile = paste0(output_path,"prs_ldpred_chr",chrnum, "_LDpred", pthrname, ".txt")
      if(file.exists(temfile)){
        dftem = read.table(temfile,sep=",",header=T)
        y = y + dftem[,'PRS']
        print(paste0('Chr ', chrnum,' Completed'))
      }
      if(!file.exists(temfile)) print(chrnum)
    }
    df.prs = data.frame(cbind(dftem[,'IID'], y))
    colnames(df.prs) = c("id", "prs")
    write.table(df.prs, paste0(writescore_path,"PRS",pthrname,".txt"), row.names = F,col.names = T, quote = FALSE, sep = "\t" )
  }
  
  ####### all patients, adjusting for age, sex and prs
  if (population == 'EA'){
    load(paste0(validate.dir,"egfr_pc_withgwasid.RData"))
    names.cols = c('gwasid','white','female','age_v1','egfrcr_v1',paste0('PC',1:10))
    validatetable = ophe_withiid[complete.cases(ophe_withiid[,names.cols]),names.cols]
    
    colnames(validatetable) = c('id','white','female','age',trait,paste0('PC',1:10))
    validatetable[['id']] = as.character(validatetable[['id']])
    validatetable[['female']] = as.factor(validatetable[['female']])
    R2.type = c('R2 Unadjusted','R2 Adjusted')
    output = matrix(NA,length(fractions),length(R2.type)+2)
    colnames(output) = c(R2.type,'Regression Coeff','P-value')
    rownames(output) = paste0(fractions,' Causal')  
  }
  
  if (population == 'AA'){
    ophe_withiid = read.csv(paste0(validate.dir,"egfr_pc_AA_withgwasid.csv"), header=T)
    names.cols = c('gwasid','white', 'female',as.vector(sapply(1,function(x){paste0(c('age_v','egfrcr_v'),x)})), 
                   paste0('PC',1:10))
    validatetable = ophe_withiid[complete.cases(ophe_withiid[,names.cols]),names.cols]
    
    colnames(validatetable) = c('id','white','female','age',trait,paste0('PC',1:10))
    validatetable[['id']] = as.character(validatetable[['id']])
    validatetable[['female']] = as.factor(validatetable[['female']])
    R2.type = c('R2 Unadjusted','R2 Adjusted')
    output = matrix(NA,length(fractions),length(R2.type)+2)
    colnames(output) = c(R2.type,'Regression Coeff','P-value')
    rownames(output) = paste0(fractions,' Causal')
  }
  
  for (i in 1:length(fractions)){
    pthr = fractions[i]
    if(pthr==1) pthrname = "_p1.0000e+00"
    if(pthr==0) pthrname = "-inf"
    if(pthr==0.1) pthrname = "_p1.0000e-01"
    if(pthr==0.3) pthrname = "_p3.0000e-01"
    if(pthr==0.01) pthrname = "_p1.0000e-02"
    if(pthr==0.03) pthrname = "_p3.0000e-02"
    if(pthr==4e-3) pthrname = "_p4.0000e-03"
    if(pthr==1e-3) pthrname = "_p1.0000e-03"
    if(pthr==3e-3) pthrname = "_p3.0000e-03"
    if(pthr==1e-4) pthrname = "_p1.0000e-04"
    if(pthr==3e-4) pthrname = "_p3.0000e-04"
    
    tem = read.table(paste0(writescore_path,"PRS",pthrname,".txt"), header = T)
    tem$id = as.character(tem$id)
    prstable = validatetable %>% inner_join(tem, by = 'id')
    prstable = prstable[complete.cases(prstable),]
    prstable[,trait] = log(prstable[,trait])
    prstable$prs = scale(prstable$prs,center=T,scale=T)
    prstable[,trait] = scale(prstable[,trait],center=T,scale=T)
    ##### unadjusted
    output[i,'R2 Unadjusted'] = (cor(prstable[,trait],prstable[,'prs']))^2
    formula.prs=formula(paste0(paste(trait, paste(c('prs','age','female',paste0('PC',1:10)),collapse="+"), sep='~')))
    fit = lm(formula.prs, data=prstable)
    #rsq.model = summary(fit)$r.squared
    output[i,'Regression Coeff'] = coefficients(fit)['prs']
    output[i,'R2 Adjusted'] = (coefficients(fit)['prs'])^2/var(na.omit(prstable[,trait]))
    output[i,'P-value'] = summary(fit)$coefficients['prs','Pr(>|t|)']
  }
  
  # p=0.3
  save(output,file=paste0('/dcl01/chatterj/data/jin/collaboration/CKDPRS/validation.ldpred.aric.',population,'.trans.RData'))
  
}
  # EUR
  # R2 Unadjusted R2 Adjusted Regression Coeff       P-value
  # 0.3 Causal    0.06751751  0.07099333        0.2664457 4.603932e-174
  # wrong 
  # 0.3 Causal    0.06656798  0.07341369       -0.2709496 1.35798e-170
  
  # AFR
  # R2 Unadjusted R2 Adjusted Regression Coeff       P-value
  # 0.3 Causal    0.01618915  0.01808404        0.1344769 8.571614e-12 # UKB GWAS without PC
  # 0.3 Causal    0.01702559  0.01650662        0.1284781 9.640903e-12
  # wrong 
  # 0.3 Causal    0.03553319   0.0410682        0.2026529 1.194893e-25
  
  
  
for (population in c('EA', 'AA')){
  validate.dir = paste0('/dcl01/chatterj/data/zyu/')
  output_path = paste0('/dcs04/nilanjan/data/jjin/UKB/transethnic/aric/',population,'/')
  writescore_path = paste0('/dcs04/nilanjan/data/jjin/UKB/transethnic/aric-prs/',population,'/')
  if (!dir.exists(writescore_path)){dir.create(writescore_path)}
  #---------------------------------------#---------------------------------------
  # input variables !!!!!!!!!!!!!!!!!!!
  #---------------------------------------#---------------------------------------
  Pvalthr = c(0.05,5e-8) #c(0.5, 5e-8)#c(1,5e-1,5e-2,5e-4,5e-6,5e-8)
  R2 = 0.1 #c(0.1,0.2,0.4,0.6,0.8)
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # freeze below
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  for(pvalthr in Pvalthr){
    if(pvalthr==1) pthrname = "_p1.0000e+00"
    if(pvalthr==5e-1) pthrname = "_p5.0000e-01"
    if(pvalthr==5e-2) pthrname = "_p5.0000e-02"
    if(pvalthr==5e-4) pthrname = "_p5.0000e-04"
    if(pvalthr==5e-6) pthrname = "_p5.0000e-06"
    if(pvalthr==5e-8) pthrname = "_p5.0000e-08"
    for (r2 in R2){
      if(r2==0.1) r = "_r0.10"
      if(r2==0.2) r = "_r0.20"
      if(r2==0.4) r = "_r0.40"
      if(r2==0.6) r = "_r0.60"
      if(r2==0.8) r = "_r0.80"
      #---------------------------------------#---------------------------------------
      y = 0;
      for(chrnum in c(1:22)){
        temfile = paste0(output_path,"prs_ldpred_chr",chrnum, "_P+T",r,pthrname, ".txt")
        if(file.exists(temfile)){
          dftem = read.table(temfile,sep=",",header=T)
          y = y + dftem[,'PRS']
          print(paste0('Chr ', chrnum,' Completed'))
        }
        if(!file.exists(temfile)) print(chrnum)
      }
      df.prs = data.frame(cbind(dftem[,'IID'], y))
      colnames(df.prs) = c("id", "prs")
      write.table(df.prs, paste0(writescore_path,"PRS_pt_",r,pthrname,".txt"), row.names = F,col.names = T, quote = FALSE, sep = "\t" )
    }
  }
  ####### all patients, adjusting for age, sex and prs
  if (population == 'EA'){
    load(paste0(validate.dir,"egfr_pc_withgwasid.RData"))
    names.cols = c('gwasid','white','female','age_v1','egfrcr_v1',paste0('PC',1:10))
    validatetable = ophe_withiid[complete.cases(ophe_withiid[,names.cols]),names.cols]
    
    colnames(validatetable) = c('id','white','female','age',trait,paste0('PC',1:10))
    validatetable[['id']] = as.character(validatetable[['id']])
    validatetable[['female']] = as.factor(validatetable[['female']])
    R2.type = c('R2 Unadjusted','R2 Adjusted')
    output = matrix(NA,length(Pvalthr)*length(R2),length(R2.type)+2)
    colnames(output) = c(R2.type,'Regression Coeff','P-value')
    settings = expand.grid(pvalthr=Pvalthr,r2=R2)
    rownames(output) = sapply(1:(length(Pvalthr)*length(R2)),function(x){paste0('pthr=',settings[x,1],' r2=', settings[x,2])})
  }
  
  if (population == 'AA'){
    ophe_withiid = read.csv(paste0(validate.dir,"egfr_pc_AA_withgwasid.csv"), header=T)
    names.cols = c('gwasid','white', 'female',as.vector(sapply(1,function(x){paste0(c('age_v','egfrcr_v'),x)})), 
                   paste0('PC',1:10))
    validatetable = ophe_withiid[complete.cases(ophe_withiid[,names.cols]),names.cols]
    
    colnames(validatetable) = c('id','white','female','age',trait,paste0('PC',1:10))
    validatetable[['id']] = as.character(validatetable[['id']])
    validatetable[['female']] = as.factor(validatetable[['female']])
    R2.type = c('R2 Unadjusted','R2 Adjusted')
    output = matrix(NA,length(Pvalthr)*length(R2),length(R2.type)+2)
    colnames(output) = c(R2.type,'Regression Coeff','P-value')
    settings = expand.grid(pvalthr=Pvalthr,r2=R2)
    rownames(output) = sapply(1:(length(Pvalthr)*length(R2)),function(x){paste0('pthr=',settings[x,1],' r2=', settings[x,2])})
  }
  
  
  for (i in 1:nrow(settings)){
    pvalthr = settings[i,'pvalthr']
    r2 = settings[i,'r2']
    if(pvalthr==1) pthrname = "_p1.0000e+00"
    if(pvalthr==5e-1) pthrname = "_p5.0000e-01"
    if(pvalthr==5e-2) pthrname = "_p5.0000e-02"
    if(pvalthr==5e-4) pthrname = "_p5.0000e-04"
    if(pvalthr==5e-6) pthrname = "_p5.0000e-06"
    if(pvalthr==5e-8) pthrname = "_p5.0000e-08"
    if(r2==0.2) r = "_r0.20"
    if(r2==0.4) r = "_r0.40"
    if(r2==0.6) r = "_r0.60"
    if(r2==0.8) r = "_r0.80"
    
    tem = read.table(paste0(writescore_path,"PRS_pt_",r,pthrname,".txt"), header = T)
    tem$id = as.character(tem$id)
    prstable = validatetable %>% inner_join(tem, by = 'id')
    prstable = prstable[complete.cases(prstable),]
    prstable[,trait] = log(prstable[,trait])
    prstable$prs = scale(prstable$prs,center=T,scale=T)
    prstable[,trait] = scale(prstable[,trait],center=T,scale=T)
    
    ##### unadjusted
    output[i,'R2 Unadjusted'] = (cor(prstable[,trait],prstable[,'prs']))^2
    ##### adjust for age and sex:
    ##### Currently the handling of covariates is not great in LDpred. 
    ##### the author suggests that we rather load the resulting scores 
    ##### and covariates in R, and get final variance explained using glm or lm in R.
    formula.prs=formula(paste0(paste(trait, paste(c('prs','age','female',paste0('PC',1:10)),collapse="+"), sep='~')))
    fit = lm(formula.prs, data=prstable)
    #rsq.model = summary(fit)$r.squared
    output[i,'Regression Coeff'] = coefficients(fit)['prs']
    output[i,'R2 Adjusted'] = (coefficients(fit)['prs'])^2/var(na.omit(prstable[,trait]))
    output[i,'P-value'] = summary(fit)$coefficients['prs','Pr(>|t|)']
  }
  
  # 0.5 0.1
  save(output,file=paste0('/dcl01/chatterj/data/jin/collaboration/CKDPRS/validation.pt.aric.',population,'.trans.RData'))
}
  # R2 Unadjusted R2 Adjusted Regression Coeff       P-value
  # EA
  # pthr=0.05 r2=0.1     0.05624631  0.05978826        0.2445164 3.652551e-146
  # pthr=5e-08 r2=0.1    0.04758893  0.04846052        0.2201375 7.495700e-119
  # wrong EA
  # pthr=0.5 r2=0.1      0.05388124  0.06186027       -0.2487172 4.685379e-141
  # pthr=5e-08 r2=0.1    0.04800085  0.04853500       -0.2203066 7.404133e-119
  
  # AA
  #                   R2 Unadjusted R2 Adjusted Regression Coeff      P-value
  # pthr=0.05 r2=0.1    0.008763627 0.009267574       0.09626824 3.703234e-07
  # pthr=5e-08 r2=0.1   0.003542102 0.005087775       0.07132864 9.920039e-05
  # wrong AA
  # pthr=0.5 r2=0.1     0.020954625 0.027778424       0.16666861 1.102042e-16
  # pthr=5e-08 r2=0.1   0.003221319 0.003994378       0.06320109 5.665135e-04


