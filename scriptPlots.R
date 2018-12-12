###############
### Load
###############
library (RColorBrewer)
#load ('~/Dropbox/BED_files_miRSNPS/databases/navyData.RData')
load ('~/Dropbox/BED_files_miRSNPS/databases/navyResultsData.RData')
load ('~/Dropbox/BED_files_miRSNPS/TCGA/datasets/BRCA_gxp_samps_mirids.RData')
source ('~/Dropbox/BED_files_miRSNPS/TCGA/script/Sfunctions.R')

###############
### n miRNAs
### by percent
###############

##----without intersect
barPlotNas (gexp = gxmir)

##----with intersect
tablemiR <- read.csv (file = '~/Dropbox/BED_files_miRSNPS/databases/SNPmiRGene.csv', 
                      header = T, sep = '\t', stringsAsFactors = F)
##----intersection miRs
miRNavy <- tablemiR$miRNA
barPlotNas (gexp = gxmir, miRNavy = miRNavy)

###############
### Heatmap
### cor
###############
###############
### Select miRs
###by %NAs
###############
tgxmir <- gxmir
gxmir <- naFilter (gexp = gxmir, maxPercentNas = 15, byNa = TRUE)
mirids <- rownames (gxmir)
miRNavy <- tablemiR$miRNA
miRInter <- intersect (mirids, miRNavy)

##----intersection Genes
gensInter <- unique (tablemiR$Predicted.Target.Gene)
tp <- gensInter[which(gensInter %in% geneidsComplete$SYMBOL)]
gensInter <- geneidsComplete[which(geneidsComplete$SYMBOL %in% tp),"ENTREZ"]
names (gensInter) <- tp

##----gexp intersections
ids <- colnames (gxmir)
gxrna <- gxrna[,ids]
gexp <- rbind (gxrna[gensInter, ], gxmir[miRInter, ])

##----correlation matrix
gexp <- t (gexp)
cormat <- cor (gexp[, miRInter], gexp[, gensInter], method = "spearman", use =  "complete.obs")
colnames (cormat) <- names (gensInter)
heatmap.cormat(cormat)
