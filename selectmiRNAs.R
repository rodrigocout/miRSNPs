############################################
### Load
############################################
#load ('~/Dropbox/BED_files_miRSNPS/databases/navyData.RData')
load ('~/Dropbox/BED_files_miRSNPS/databases/navyResultsData.RData')
load ('~/Dropbox/BED_files_miRSNPS/TCGA/datasets/BRCA_gxp_samps_mirids.RData')
library (RColorBrewer)
library (RTN)
library(gplots)
tablemiR <- read.csv (file = '~/Dropbox/BED_files_miRSNPS/databases/SNPmiRGene.csv', header = T, sep = '\t', stringsAsFactors = F)

############################################
###Select genes
############################################
vecGwasS <- c('ROPN1L', 'TERC', 'BCL2L15', 'AP4B1', 'DCLRE1B', 'HIPK1', 'PTPN22', 'BCL2L15', 'AP4B1', 'DCLRE1B',
              'HIPK1', 'PTPN22', 'MAP3K1', 'PRKG1', 'CSTF2T', 'STAM2', 'ESR1', 'C6orf97', 'ECHDC1', 'RNF146',
              'ATP2B4', 'MDM4', 'SSBP4', 'ISYNA1', 'ELL', 'COX11', 'SLC4A7', 'ANKLE1', 'C19orf62', 'ABHD8', 
              'PRKACB') 
vecGwasS <- unique (vecGwasS)
symbolPredicted <- unique (tablemiR$Predicted.Target.Gene)

genes <- unique (c(vecGwasS, symbolPredicted))
genes <- genes[which(genes %in% geneidsComplete$SYMBOL)]
##---
idx <- which (geneidsComplete$SYMBOL %in% genes)
entrz <- geneidsComplete[idx, "ENTREZ"]
names (entrz) <- genes

############################################
###Select miRNAs
############################################
mirs <- mirids[which(mirids %in% unique(tablemiR$miRNA))]
##---
