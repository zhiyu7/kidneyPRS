rm(list=ls())
library(data.table)
library(dplyr)

traitlist= c("eGFRcr-LDPRED")
Sys.setenv(HDF5_USE_FILE_LOCKING='FALSE')
Sys.setenv(OPENBLAS_NUM_THREADS="1")

temp <- commandArgs(TRUE)
trait_name =  traitlist[as.numeric(temp[1])]
chr  = as.numeric(temp[2])
#trait_name = 'eGFRcr-LDPRED'
#chr=1
ssf.format = 'LDPRED'

population = 'EA'
if (population == 'EA') pop = 'EUR'
if (population == 'AA') pop = 'AFR'
#---------------------------------------#---------------------------------------


#---------------------------------------#---------------------------------------
# input variables !!!!!!!!!!!!!!!!!!!
#---------------------------------------#---------------------------------------
ld_radius = 400 # Jin: (1million / 3000)
#---------------------------------------#---------------------------------------

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# freeze below
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

maindirec = "/dcl01/chatterj/data/yzhang/"
jindir = "/dcl01/chatterj/data/jin/"

plinkpath = paste0(maindirec,"software/plink")
output_path = paste0("/dcs04/nilanjan/data/jjin/UKB/transethnic/ldpred/")
if (!dir.exists(output_path)){dir.create(output_path)}

LDrefpanel_plink_path = paste0('/dcs04/nilanjan/data/jjin/UKB/transethnic/ref/chr',chr)
validationUKBiobank_path = paste0('/dcs04/nilanjan/data/jjin/UKB/transethnic/validation/chr',chr)
sumstatsLDpred_path = paste0('/dcs04/nilanjan/data/jjin/UKB/transethnic/ldfile/')

tmp_output_path = paste0('/dcs04/nilanjan/data/jjin/UKB/transethnic/tmp/')
if (!dir.exists(tmp_output_path)){dir.create(tmp_output_path)}

#-------------------------------------------------------------------
linuxcode = paste(paste0("if [ -e ", tmp_output_path,"/coord_chr", chr, ".hdf5 ];"),
                  paste0("then rm ", tmp_output_path,"/coord_chr", chr, ".hdf5;"),
                  paste0("echo FILE FOUND AND REMOVED;"),
                  "fi")
system(linuxcode)
linuxcode = paste(paste0("if [ -e ", sumstatsLDpred_path,"/ldfilename* ];"),
                  paste0("then rm ", sumstatsLDpred_path,"/ldfilename*;"),
                  paste0("echo FILE FOUND AND REMOVED;"),
                  "fi")
system(linuxcode)



if (ssf.format == 'LDPRED'){
  # Run the coord file to generate hdf5 file
  # hdf5 file create is only possible in /tmp directory
  coordcode = paste(
    "/users/jjin/.local/bin/ldpred",
    paste0("--debug"),
    paste0("coord"),
    paste0("--vbim ",validationUKBiobank_path,".bim"),
    paste0("--gf=",LDrefpanel_plink_path),
    # only select the SNPs that are among the 1.2M HapMap3 SNPs
    #"--only-hm3",
    #paste0("--ssf=/dcl01/chatterj/data/zyu/summary_meta_CKDGen_UKB_1.2m"), 
    #paste0("--ssf=/dcl01/chatterj/data/zyu/summary_meta_CKDGen_UKB_1.2m_highprecision_10digit.txt"),
    #paste0("--ssf=/dcl01/chatterj/data/zyu/multiethnic_ckdgen/eGFRcr/fromJin/meta-ukbtrans-ckdgentrans.txt"),
    paste0('--ssf=/dcs04/nilanjan/data/jjin/UKB/transethnic/sumdata/mega-meta-ukbtrans-ckdgentrans.txt'),
    "--ssf-format=LDPRED",
    paste0("--out=",tmp_output_path, "/coord_chr",chr,".hdf5")
  )
  system(coordcode)
}


# Let's make sure this actually ran
system("echo Coord ran successfully")
print("-----------------")

Pvalthr = c(1,0.5,0.05,0.0005,0.000005,0.00000005)
#Pvalthr = 0.00000005 #c(0.05,0.0005,0.000005,0.00000005)
#R2 = 0.1 #c(0.2,0.4,0.6,0.8)
R2 = c(0.1,0.2,0.4,0.6,0.8)
for (pvalthr in Pvalthr){
  for (r2 in R2){
    #-------------------------------------------------------------------
    # Using hdf5 file generated from above, run ldpred
    # Remember to copy all result files to my directory.
    ldpredcode = paste(
      "/users/jjin/.local/bin/ldpred p+t",
      paste0("--cf=",tmp_output_path, "/coord_chr",chr,".hdf5"),
      paste0("--ldr=", ld_radius), 
      paste0("--p=", pvalthr),
      paste0("--r2=", r2),
      #paste0("--ldf=",sumstatsLDpred_path, "/ldfilename_chr",chr),
      paste0("--out=",tmp_output_path, "/ldpredout_chr",chr)
    )
    system(ldpredcode)
    
    # Let's make sure this actually ran
    system("echo LDpred ran successfully")
    print("-----------------")
    
    #-------------------------------------------------------------------
    # Copy the files to my directory
    linuxcode = paste("mv",
                      paste0(tmp_output_path, "/ldpredout_chr",chr,"*"), 
                      paste0(output_path,"/"))
    system(linuxcode)
    
    #-------------------------------------------------------------------
    # validate step to calculate PRS
    validatecode = paste(
      #"/jhpce/shared/jhpce/core/python/2.7.9/bin/python -m",
      "/users/jjin/.local/bin/ldpred",
      paste0("--debug"),
      paste0("score"),
      paste0("--gf=", validationUKBiobank_path), 
      paste0("--rf=",  output_path, "/ldpredout_chr",chr), 
      paste0("--rf-format=P+T"),
      paste0("--pf-format=STANDARD"),
      paste0("--p=", pvalthr),
      paste0("--r2=", r2),
      paste0("--only-score"),
      paste0("--out=", output_path,"/prs_ldpred_","chr",chr)
    )
    system(validatecode)
    
    #system(paste0("rm ", tmp_output_path,"/prs_ldpred_","chr",chr,"*.hdf5"))
    system("echo validate ran successfully")
    print("-----------------")
  }
}

