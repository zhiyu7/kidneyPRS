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
library(bigreadr)
library(MendelianRandomization) # for mr_ivw
library(dplyr)
library(R.utils) # for gzip
library(stringr) # for str_detect
library(genio) # a package to facilitate reading and writing genetics data. The focus of this vignette is processing plink BED/BIM/FAM files.

trait = 'eGFRcr-LDPRED'
output.dir = paste0('/dcl01/chatterj/data/jin/collaboration/CKDPRS/data/validation/',trait,'/')
#rawdata = readRDS(paste0(output.dir,trait,".rds")) # 145634 x 15

ukbdat = bigreadr::fread2(paste0('~/Dropbox/JHU/PRS-Yu/code/gwasdata_raw_for_lm.txt'))
#ind.eur = which(ukbdat$race == 'EUR')
#ind.others = which(ukbdat$race != 'EUR')
#set.seed(2021)
#keep.eur = sample(ind.eur, ceiling(length(ind.eur)/5))
#ukbdat = ukbdat[c(keep.eur, ind.others),]
#ukbdat$pc1 = scale(ukbdat$pc1); ukbdat$pc2 = scale(ukbdat$pc2)
colnames(ukbdat)[which(colnames(ukbdat) == 'pc1')] = 'PC1'
colnames(ukbdat)[which(colnames(ukbdat) == 'pc2')] = 'PC2'


library(RColorBrewer)
library(ggplot2)
library(hrbrthemes)
library(ggExtra)
library(gtable)
library(gbm)

mycols = c(brewer.pal(8,"Reds")[8], brewer.pal(8,"YlOrBr")[4], 
           brewer.pal(8,"Paired")[2], 'darkolivegreen4')

# p <- ggplot(ukbdat, aes(x=pc1, y=pc2, color = race)) +
#   geom_point() + 
#   scale_color_manual(values = mycols) + 
#   theme_ipsum() +
#   #geom_point(aes(colour = factor(race))) + #, size = 4) +
#   theme(text = element_text(size=14, face="bold"), # axis.text.x = element_text(size=14, face="bold"), 
#         panel.spacing = unit(3, "lines"), 
#         legend.position = "right", legend.title = element_blank(), 
#         panel.background = element_rect(fill = "white", colour = "white", size = 0.2, linetype = "solid"),
#         panel.grid = element_line(size = 0.5, linetype = 'solid', colour = 'black'),
#         panel.border = element_blank()) +
#   removeGrid(x = TRUE, y = TRUE)



p <- ggplot(ukbdat, aes(x=PC1, y=PC2)) +
  geom_point(pch=21, size = 2.3, color = 'black', aes(fill = race)) + 
  scale_fill_manual(values = mycols) + 
  #theme_ipsum() +
  #geom_point(aes(colour = factor(race))) + #, size = 4) +
  theme(text = element_text(size=13), # axis.text.x = element_text(size=14, face="bold"), 
        axis.title.x = element_text(size=11, face="bold"), axis.title.y = element_text(size=11, face="bold"),
        #panel.spacing = unit(3, "lines"), 
        legend.position = "right", legend.title = element_blank(), 
        legend.text = element_text(size = 12, face = 'bold'),
        panel.background = element_rect(fill = "white", colour = "white", size = 0.2, linetype = "solid"),
        panel.grid = element_line(size = 1, linetype = 'solid', colour = 'black'), 
        panel.border = element_rect(colour = 'black', fill=NA, size=1)) +
  removeGrid(x = TRUE, y = TRUE)

png(paste0('~/Dropbox/JHU/PRS-Yu/code/PCAplots.png'),width = 2000, height = 1600, res = 315)
p
# grid.arrange(grobs = list(p1,p2,leg),
#              widths = c(4), heights = c(1,4,4), #nrow=3,ncol=2, 
#              layout_matrix = rbind(3,1,2))
dev.off()

