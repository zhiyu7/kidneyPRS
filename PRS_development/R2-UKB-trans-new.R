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

#K=5 # the number of genetic PCs, not necessary for eGFRcr
#---------------------------------------#---------------------------------------
# input variables !!!!!!!!!!!!!!!!!!!
#---------------------------------------#---------------------------------------

# ldpred
validate.dir = paste0('/dcs04/nilanjan/data/jjin/UKB/transethnic/validation/pheno/')
# note: the output was save in the wrong path..
output_path = paste0('/dcs04/nilanjan/data/jjin/UKB/transethnic/pt/')
writescore_path = paste0('/dcs04/nilanjan/data/jjin/UKB/transethnic/prs/')
if (!dir.exists(writescore_path)){dir.create(writescore_path)}

  fractions = c(1,0.3,0.1,0.03,0.01,0.003,0.001)#,GWASalphathr
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
  
  ####### all patients, adjusting for age and sex AND 10 PC
  validatetable = bigreadr::fread2(paste0('/dcs04/nilanjan/data/jjin/UKB/transethnic/pheno/valpheno.txt'))
  validatetable = validatetable[,c('eid','age','sex','race',traitcolname,paste0('pc',1:40))]
  # remove the UKB training patients 
  # validationid = read.table(paste0(output.dir,'UKBvalidation.idlist.txt'),header=F)
  # validationid = validationid[,1]
  # validatetable = validatetable[(validatetable[,'eid'] %in% validationid),] # 32159
  trait = 'eGFRcr'
  colnames(validatetable) = c('id','age','sex','race',trait,paste0('pc',1:40))
  R2.type = c('R2 Unadjusted','R2 Adjusted')
  output = matrix(NA,length(fractions),length(R2.type)+2)
  colnames(output) = c(R2.type,'Regression Coeff','P-value')
  rownames(output) = paste0(fractions,' Causal')
  
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
    tem$prs = scale(tem$prs,center=T,scale=T)
    prstable = validatetable %>% inner_join(tem, by = 'id')
    prstable = prstable[complete.cases(prstable),]
    prstable[,trait] = log(prstable[,trait])
    
    ##### unadjusted
    output[i,'R2 Unadjusted'] = (cor(prstable[,trait],prstable[,'prs']))^2
    ##### adjust for age and sex:
    ##### Currently the handling of covariates is not great in LDpred. 
    ##### the author suggests that we rather load the resulting scores 
    ##### and covariates in R, and get final variance explained using glm or lm in R.
    formula.prs=formula(paste0(paste(trait, paste(c('prs','age','sex','race',paste0('pc',1:40)),collapse="+"), sep='~')))
    fit = lm(formula.prs, data=prstable)
    #rsq.model = summary(fit)$r.squared
    output[i,'Regression Coeff'] = coefficients(fit)['prs']
    output[i,'R2 Adjusted'] = (coefficients(fit)['prs'])^2/var(na.omit(prstable[,trait]))
    output[i,'P-value'] = summary(fit)$coefficients['prs','Pr(>|t|)']
  }
  # p = 0.3
  save(output,file='/dcl01/chatterj/data/jin/collaboration/CKDPRS/validation.ldpred.ukbtrans.ckdgentrans.RData')

  # wrong
  # R2 Unadjusted R2 Adjusted Regression Coeff       P-value
  # 1 Causal       0.052499679 0.062720208     -0.010409752  0.000000e+00
  # 0.3 Causal     0.055273080 0.065552988     -0.010642236  0.000000e+00
  # 0.1 Causal     0.054170287 0.060908662     -0.010258318  0.000000e+00
  # 0.03 Causal    0.030563225 0.030516135     -0.007261083 3.158906e-250
  # 0.01 Causal    0.012988607 0.012506758     -0.004648459 6.410941e-105
  # 0.003 Causal   0.005762951 0.005374099     -0.003047119  4.496619e-46
  # 0.001 Causal   0.003680273 0.003296068     -0.002386353  6.904683e-29

  # correct
  # R2 Unadjusted R2 Adjusted Regression Coeff       P-value
  # 1 Causal       0.052454071 0.062861568      0.010421476  0.000000e+00
  # 0.3 Causal     0.055068049 0.065654174      0.010650446  0.000000e+00
  # 0.1 Causal     0.054975606 0.062574040      0.010397615  0.000000e+00
  # 0.03 Causal    0.032250031 0.032373688      0.007478814 1.831096e-264
  # 0.01 Causal    0.012439101 0.011990055      0.004551423 1.228898e-100
  # 0.003 Causal   0.005150372 0.004803787      0.002880902  2.326890e-41
  # 0.001 Causal   0.003342578 0.002973002      0.002266388  3.258497e-26
  
  
  
  
  
if (prsmethod == 'p+t')
{
  validate.dir = paste0('/dcs04/nilanjan/data/jjin/UKB/transethnic/validation/pheno/')
  output_path = paste0('/dcs04/nilanjan/data/jjin/UKB/transethnic/ldpred/')
  writescore_path = paste0('/dcs04/nilanjan/data/jjin/UKB/transethnic/prs/')
  if (!dir.exists(writescore_path)){dir.create(writescore_path)}

  Pvalthr = c(1,5e-1,5e-2,5e-4,5e-6,5e-8)
  R2 = c(0.1,0.2,0.4,0.6,0.8)
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
        #prs_ldpred_chr9_P+T_r0.40_p1.0000e+00.txt
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
    
  
  ####### all patients, adjusting for age and sex
  ####### all patients, adjusting for age and sex AND 10 PC
  validatetable = bigreadr::fread2(paste0('/dcs04/nilanjan/data/jjin/UKB/transethnic/pheno/valpheno.txt'))
  validatetable = validatetable[,c('eid','age','sex','race',traitcolname,paste0('pc',1:40))]
  # remove the UKB training patients 
  # validationid = read.table(paste0(output.dir,'UKBvalidation.idlist.txt'),header=F)
  # validationid = validationid[,1]
  # validatetable = validatetable[(validatetable[,'eid'] %in% validationid),] # 32159
  trait = 'eGFRcr'
  colnames(validatetable) = c('id','age','sex','race',trait,paste0('pc',1:40))
  R2.type = c('R2 Unadjusted','R2 Adjusted')
  output = matrix(NA,length(Pvalthr)*length(R2),length(R2.type)+2)
  colnames(output) = c(R2.type,'Regression Coeff','P-value')
  settings = expand.grid(pvalthr=Pvalthr,r2=R2)
  rownames(output) = sapply(1:(length(Pvalthr)*length(R2)),function(x){paste0('pthr=',settings[x,1],' r2=', settings[x,2])})
  
  for (i in 1:nrow(settings)){
    pvalthr = settings[i,'pvalthr']
    r2 = settings[i,'r2']
    if(pvalthr==1) pthrname = "_p1.0000e+00"
    if(pvalthr==5e-1) pthrname = "_p5.0000e-01"
    if(pvalthr==5e-2) pthrname = "_p5.0000e-02"
    if(pvalthr==5e-4) pthrname = "_p5.0000e-04"
    if(pvalthr==5e-6) pthrname = "_p5.0000e-06"
    if(pvalthr==5e-8) pthrname = "_p5.0000e-08"
    if(r2==0.1) r = "_r0.10"
    if(r2==0.2) r = "_r0.20"
    if(r2==0.4) r = "_r0.40"
    if(r2==0.6) r = "_r0.60"
    if(r2==0.8) r = "_r0.80"
    
    ##### unadjusted
    tem = read.table(paste0(writescore_path,"PRS_pt_",r,pthrname,".txt"), header = T)
    prstable = validatetable %>% inner_join(tem, by = 'id')
    prstable = prstable[complete.cases(prstable),]
    prstable[,trait] = log(prstable[,trait])
    prstable$prs = scale(prstable$prs,center=T,scale=T)
    
    ##### unadjusted
    output[i,'R2 Unadjusted'] = (cor(prstable[,trait],prstable[,'prs']))^2
    ##### adjust for age and sex:
    ##### Currently the handling of covariates is not great in LDpred. 
    ##### the author suggests that we rather load the resulting scores 
    ##### and covariates in R, and get final variance explained using glm or lm in R.
    formula.prs=formula(paste0(paste(trait, paste(c('prs','age','sex','race',paste0('pc',1:40)),collapse="+"), sep='~')))
    fit = lm(formula.prs, data=prstable)
    output[i,'Regression Coeff'] = coefficients(fit)['prs']
    output[i,'R2 Adjusted'] = (coefficients(fit)['prs'])^2/var(na.omit(prstable[,trait]))
    output[i,'P-value'] = summary(fit)$coefficients['prs','Pr(>|t|)']
  }
  # 0.05, 0.1
  save(output,file='/dcl01/chatterj/data/jin/collaboration/CKDPRS/ukb.validation.pt.transukb.transckdgen.RData')
}
