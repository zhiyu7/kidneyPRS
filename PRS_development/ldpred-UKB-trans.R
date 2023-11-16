rm(list=ls())
library(data.table)
library(dplyr)

traitlist= c("eGFRcr-LDPRED")
Sys.setenv(HDF5_USE_FILE_LOCKING='FALSE')
Sys.setenv(OPENBLAS_NUM_THREADS="1")

temp <- commandArgs(TRUE)
trait_name =  traitlist[as.numeric(temp[1])]
chr  = as.numeric(temp[2])
ssf.format = 'LDPRED'
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
output_path = paste0("/dcs04/nilanjan/data/jjin/UKB/transethnic/pt/")
if (!dir.exists(output_path)){dir.create(output_path)}

LDrefpanel_plink_path = paste0('/dcs04/nilanjan/data/jjin/UKB/transethnic/ref/chr',chr)
validationUKBiobank_path = paste0('/dcs04/nilanjan/data/jjin/UKB/transethnic/validation/chr',chr)
sumstatsLDpred_path = paste0('/dcs04/nilanjan/data/jjin/UKB/transethnic/ldfileldpred/')

fractions = c(1,0.3,0.1,0.03,0.01,0.003,0.001)

tmp_output_path = paste0('/dcs04/nilanjan/data/jjin/UKB/transethnic/tmpldpred')
if (!dir.exists(tmp_output_path)){dir.create(tmp_output_path)}

#-------------------------------------------------------------------
linuxcode = paste(paste0("if [ -e ", tmp_output_path,"/coord_chr", chr, ".hdf5 ];"),
                  paste0("then rm ", tmp_output_path,"/coord_chr", chr, ".hdf5;"),
                  paste0("echo FILE FOUND AND REMOVED;"),
                  "fi")
system(linuxcode)
linuxcode = paste(paste0("if [ -e ", sumstatsLDpred_path,"/ldfilename_chr",chr,"_ldradius400.pkl.gz ];"),
                  paste0("then rm ", sumstatsLDpred_path,"/ldfilename_chr",chr,"_ldradius400.pkl.gz;"),
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
    paste0("--gf=",LDrefpanel_plink_path),
    paste0("--vbim ",validationUKBiobank_path,".bim"),
    #paste0("--ssf=/dcl01/chatterj/data/zyu/multiethnic_ckdgen/eGFRcr/fromJin/meta_ckdgen_trans_ukb_wholegenome.txt"),
    #paste0("--ssf=/dcl01/chatterj/data/zyu/multiethnic_ckdgen/eGFRcr/fromJin/mega-meta-ukbtrans-ckdgentrans.txt"),
    paste0('--ssf=/dcs04/nilanjan/data/jjin/UKB/transethnic/sumdata/mega-meta-ukbtrans-ckdgentrans.txt'),
    "--ssf-format=LDPRED",
    paste0("--out=",tmp_output_path, "/coord_chr",chr,".hdf5")
  )
  system(coordcode)
}

# Let's make sure this actually ran
system("echo Coord ran successfully")
print("-----------------")

#-------------------------------------------------------------------
# Using hdf5 file generated from above, run ldpred
# Remember to copy all result files to my directory.
ldpredcode = paste(
  "/users/jjin/.local/bin/ldpred fast",
  paste0("--cf=",tmp_output_path, "/coord_chr",chr,".hdf5"),
  paste0("--ldr=", ld_radius), 
  paste0("--ldf=",sumstatsLDpred_path, "/ldfilename_chr",chr),
  #paste0("--f=0.004"), 
  #paste0("--N=",samplesize),
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



for (i in 1:length(fractions))
{
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
  
  temfile = paste0(output_path,"/ldpredout_chr",chr, "_LDpred_fast", pthrname, ".txt")
  if(file.exists(temfile)){
    linuxcode = paste("mv",
                      paste0(output_path,"/ldpredout_chr",chr, "_LDpred_fast", pthrname, ".txt"), 
                      paste0(output_path,"/ldpredout_chr",chr, "_LDpred", pthrname, ".txt"))
    system(linuxcode)
  }
}

#-------------------------------------------------------------------
# validate step to calculate PRS
validatecode = paste(
  #"/jhpce/shared/jhpce/core/python/2.7.9/bin/python -m",
  "/users/jjin/.local/bin/ldpred",
  paste0("--debug"),
  paste0("score"),
  paste0("--gf=", validationUKBiobank_path), 
  paste0("--rf=",  output_path, "/ldpredout_chr",chr), 
  #paste0("--f=0.004"), 
  paste0("--rf-format=LDPRED"),
  paste0("--pf-format=STANDARD"),
  paste0("--only-score"),
  paste0("--out=", output_path,"/prs_ldpred_","chr",chr)
)
system(validatecode)
#summary_dict[0]={'name':'Validation genotype file (prefix):','value':p_dict['gf']}
#summary_dict[0.1]={'name':'Input weight file(s) (prefix):','value':p_dict['rf']}
#summary_dict[0.2]={'name':'Output scores file(s) (prefix):','value':p_dict['out']}

system("echo validate ran successfully")
print("-----------------")



