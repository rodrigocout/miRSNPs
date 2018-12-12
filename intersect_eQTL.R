###Integration with eQTL data public available#####
###Blood eQTL####
##cis eQTL###
cis_eqtl <- read.delim('~/Dropbox/BED_files_miRSNPS/databases/Blood_eQTL_CisAssociationsProbeLevelFDR0.5.txt', header = T)
head(cis_eqtl)
colnames(cis_eqtl)
length(cis_eqtl$SNPName)

load(file = "~/Dropbox/BED_files_miRSNPS/databases/navyResultsData.RData")
colnames(navyD)
length(navyD$SNP)

#Intersect
inter_eqtl <- intersect(cis_eqtl$SNPName,navyD$SNP)
length(inter_eqtl)
inter_eqtl <- as.data.frame(inter_eqtl)
colnames(inter_eqtl) = 'SNPs_eqtl'

save(inter_eqtl, file = "~/Dropbox/BED_files_miRSNPS/databases/inter_Blood_cis_eQTL.RData")

###trans-eqtl#######
trans_eqtl <- read.delim('~/Dropbox/BED_files_miRSNPS/databases/Blood_eQTL_TransEQTLsFDR0.5.txt', header = T)
head(trans_eqtl)
colnames(trans_eqtl)
length(trans_eqtl$SNPName)

#Intersect
inter_trans_eqtl <- intersect(trans_eqtl$SNPName,navyD$SNP)
length(inter_trans_eqtl)
inter_trans_eqtl <- as.data.frame(inter_trans_eqtl)
colnames(inter_trans_eqtl) = 'SNPs_trans_eqtl'

save(inter_eqtl, file = "~/Dropbox/BED_files_miRSNPS/databases/inter_Blood_trans_eQTL.RData")

#Intersect GWAS with eQTLs (Blood eQTL)
gwas_ciseqtl <- intersect(inter_gwas$SNPs, inter_eqtl$SNPs_eqtl)
length(gwas_ciseqtl)
gwas_ciseqtl <- as.data.frame(gwas_ciseqtl)
head(gwas_ciseqtl)
save(gwas_ciseqtl, file = "~/Dropbox/BED_files_miRSNPS/databases/GWAS_Blood_cis_eQTL.RData")

gwas_transeqtl <- intersect(inter_gwas$SNPs, inter_trans_eqtl$SNPs_trans_eqtl)
length(gwas_transeqtl)
gwas_transeqtl <- as.data.frame(gwas_transeqtl)
save(gwas_transeqtl, file = "~/Dropbox/BED_files_miRSNPS/databases/GWAS_Blood_trans_eQTL.RData")
