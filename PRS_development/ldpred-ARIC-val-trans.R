rm(list=ls())
library(data.table)
library(dplyr)
traitlist= c("eGFRcr-LDPRED")
trait_name =  traitlist[1]
trait = 'eGFRcr'
Sys.setenv(HDF5_USE_FILE_LOCKING='FALSE')
Sys.setenv(OPENBLAS_NUM_THREADS="1")
temp <- commandArgs(TRUE)
pop = c('EA','AA')
population = pop[as.numeric(temp[2])]
#population = 'EA'
if (population == 'EA') pop = 'EUR'
if (population == 'AA') pop = 'AFR'
chrnum  = as.numeric(temp[1])
#trait_name = 'eGFRcr-LDPRED'
#chrnum=1
ssf.format = 'LDPRED'
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
#plinkpath = paste0(maindirec,"software/plink-1.07-x86_64/")
output_path = paste0('/dcs04/nilanjan/data/jjin/UKB/transethnic/aric/',population,'/')
if (!dir.exists(output_path)){dir.create(output_path)}
validation_file = paste0("/dcl01/chatterj/data/jin/prs/realdata/ARIC/eGFRcr/",pop,"/geno/mega/chr.qc",chrnum)
### the selected tuning parameter
fraction = 0.3

LDrefpanel_plink_path = paste0('/dcs04/nilanjan/data/jjin/UKB/transethnic/ref/chr')
sumstatsLDpred_path = paste0('/dcs04/nilanjan/data/jjin/UKB/transethnic/ldfile-aricldpred/')

tmp_output_path = paste0('/dcs04/nilanjan/data/jjin/UKB/transethnic/tmp-aricldpred/')
if (!dir.exists(tmp_output_path)){dir.create(tmp_output_path)}

#-------------------------------------------------------------------
linuxcode = paste(paste0("if [ -e ", tmp_output_path,"/coord_chr", chrnum, ".hdf5 ];"),
                  paste0("then rm ", tmp_output_path,"/coord_chr", chrnum, ".hdf5;"),
                  paste0("echo FILE FOUND AND REMOVED;"),
                  "fi")
system(linuxcode)
linuxcode = paste(paste0("if [ -e ", sumstatsLDpred_path,"/ldfilename_chr",chrnum,"_ldradius400.pkl.gz ];"),
                  paste0("then rm ", sumstatsLDpred_path,"/ldfilename_chr",chrnum,"_ldradius400.pkl.gz;"),
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
    paste0("--vbim ",validation_file,".bim"),
    paste0("--gf=",LDrefpanel_plink_path,chrnum),
    paste0('--ssf=/dcs04/nilanjan/data/jjin/UKB/transethnic/sumdata/mega-meta-ukbtrans-ckdgentrans.txt'),
    # only select the SNPs that are among the 1.2M HapMap3 SNPs
    "--ssf-format=LDPRED",
    paste0("--out=",tmp_output_path, "/coord_chr",chrnum,".hdf5")
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
  paste0("--cf=",tmp_output_path, "/coord_chr",chrnum,".hdf5"),
  paste0("--ldr=", ld_radius), 
  paste0("--ldf=",sumstatsLDpred_path, "/ldfilename_chr",chrnum),
  paste0("--f=",fraction), 
  #paste0("--N=",samplesize),
  paste0("--out=",tmp_output_path, "/ldpredout_chr",chrnum)
)
system(ldpredcode)

# Let's make sure this actually ran
system("echo LDpred ran successfully")
print("-----------------")


#-------------------------------------------------------------------
# Copy the files to my directory
linuxcode = paste("mv",
                  paste0(tmp_output_path, "/ldpredout_chr",chrnum,"*"), 
                  paste0(output_path,"/"))
system(linuxcode)


for (i in 1:length(fraction)){
  pthr = fraction[i]
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
  
  temfile = paste0(output_path,"/ldpredout_chr",chrnum, "_LDpred_fast", pthrname, ".txt")
  if(file.exists(temfile)){
    linuxcode = paste("mv",
                      paste0(output_path,"/ldpredout_chr",chrnum, "_LDpred_fast", pthrname, ".txt"), 
                      paste0(output_path,"/ldpredout_chr",chrnum, "_LDpred", pthrname, ".txt"))
    system(linuxcode)
  }
}
#ldpredout_chr18_LDpred_fast_p4.0000e-03.txt




#-------------------------------------------------------------------
# validate step to calculate PRS
validatecode = paste(
  "/users/jjin/.local/bin/ldpred",
  paste0("--debug"),
  paste0("score"),
  paste0("--gf=", validation_file), 
  paste0("--rf=",  output_path, "/ldpredout_chr",chrnum), 
  paste0("--f=",fraction), 
  paste0("--rf-format=LDPRED"),
  paste0("--pf-format=STANDARD"),
  paste0("--only-score"),
  paste0("--out=", output_path,"/prs_ldpred_","chr",chrnum)
)
system(validatecode)

system("echo validate ran successfully")
print("-----------------")

