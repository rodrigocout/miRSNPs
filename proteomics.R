#########################
####load
#########################
load (file="~/Dropbox/BED_files_miRSNPS/TCGA/datasets/BRCA_gxp_samps_mirids.RData")
#----readTableExpression
protTable <- read.table(file="~/Dropbox/BED_files_miRSNPS/TCGA/dados_TCGA/gdac.broadinstitute.org_BRCA.Merge_protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data.Level_3.2016012800.0.0/BRCA.protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data.data.txt", sep="\t",  header=T,  stringsAsFactors=F)

##----pre-process
##---rownames
protTable <- protTable[-1,]
rownames(protTable) <- protTable[, 1]
protTable <- protTable[, -1]
##---colnames
colnames(protTable) <- gsub(".","-",colnames(protTable), fixed=TRUE)
tp <- colnames(protTable)
tpp <- substr(tp, 1, 16)
idx <- which(tpp%in%rownames(dataFreeze))
protTable <- protTable[, idx]
tpp <- tpp[idx]
sampsIds <- colnames(protTable);names(sampsIds) <- tpp
##---update dataFreeze
sampsIds <- sampsIds[rownames(dataFreeze)]
dataFreeze <- cbind(dataFreeze, sampsIds)
colnames(dataFreeze)[3] <- "Proteomic"

##---IdsGenes
refIdsGenes <- read.table(file="~/Dropbox/BED_files_miRSNPS/TCGA/dados_TCGA/gdac.broadinstitute.org_BRCA.RPPA_AnnotateWithGene.Level_3.2016012800.0.0/BRCA.antibody_annotation.txt", sep="\t",  stringsAsFactors=F,  header=T)

#########################
####naFilter
#########################
source(file = '~/Dropbox/BED_files_miRSNPS/TCGA/script/Sfunctions.R')
#----rodrigoTable
tablemiR <- read.csv (file =
                          '~/Dropbox/BED_files_miRSNPS/databases/SNPmiRGene.csv',
                      sep = "\t",  stringsAsFactors = F)

#----filterbyNAs
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
