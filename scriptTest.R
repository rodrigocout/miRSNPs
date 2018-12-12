##############################################
####Script Test
##############################################
###############
### Load
###############
library (RColorBrewer)
#load ('~/Dropbox/BED_files_miRSNPS/databases/navyData.RData')
load ('~/Dropbox/BED_files_miRSNPS/databases/navyResultsData.RData')
load ('~/Dropbox/BED_files_miRSNPS/TCGA/datasets/BRCA_gxp_samps_mirids.RData')
source ('~/Dropbox/BED_files_miRSNPS/TCGA/script/Sfunctions.R')

###############
### Select %NAs
###############

gxmir <- naFilter (gexp = gxmir, maxPercentNas = 15, byNa = TRUE)
mirids <- rownames (gxmir)
# ###############
# ### Re-order
# ###############
# gxmirOrd <- expOrder (gexp = gxmir, decreasing = TRUE, bySample = FALSE)

###############
### miRs characteristics
###sd, var, 
###############
tableChar <- gxmirChar (gexp = gxmir)

# ###############
# ### heatMap gexp
# ###############
# heatmap.gexp (gexp = gxmirOrd)

###############
### miR GWAS
###############
#load (file = "~/Dropbox/BED_files_miRSNPS/TCGA/datasets/extraGenes.RData")
tablemiR <- read.csv (file = '~/Dropbox/BED_files_miRSNPS/databases/SNPmiRGene.csv', 
                      header = T, sep = '\t', stringsAsFactors = F)
##----intersection miRs
miRNavy <- tablemiR$miRNA
miRInter <- intersect (mirids, miRGWAS)

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

##----heatmap cormat
heatmap.cormat (cormat)
