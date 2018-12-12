############################################
### Load
############################################
library (RTN)
load ('~/Dropbox/BED_files_miRSNPS/TCGA/datasets/BRCA_gxp_samps_mirids.RData')
names (mirids) <- NULL

rtniBRCA_TCGA <- new ("TNI", gexp = gexp, transcriptionFactors = mirids)
save (rtniBRCA_TCGA, file = '~/Dropbox/BED_files_miRSNPS/TCGA/datasets/BRCAtni_may2016.RData')

rtniBRCA_TCGA <- tni.preprocess(rtniBRCA_TCGA, gexpIDs = geneidsComplete)
save (rtniBRCA_TCGA, file = '~/Dropbox/BED_files_miRSNPS/TCGA/datasets/BRCAtni_may2016.RData')

rtniBRCA_TCGA <- tni.permutation(rtniBRCA_TCGA, pValueCutoff = 1e-7)
save (rtniBRCA_TCGA, file = '~/Dropbox/BED_files_miRSNPS/TCGA/datasets/BRCAtni_may2016.RData')

rtniBRCA_TCGA <- tni.bootstrap (rtniBRCA_TCGA, nBootstraps = 200)
save (rtniBRCA_TCGA, file = '~/Dropbox/BED_files_miRSNPS/TCGA/datasets/BRCAtni_may2016.RData')

rtniBRCA_TCGA <- tni.dpi.filter(rtniBRCA_TCGA)
save (rtniBRCA_TCGA, file = '~/Dropbox/BED_files_miRSNPS/TCGA/datasets/BRCAtni_may2016.RData')
