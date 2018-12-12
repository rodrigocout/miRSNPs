####Analysis miRSNP Naive Bayes and Autoimmune########
setwd('/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/summary_statistics')

##Summary statistics from Immunobase###########
cro <- read.table('hg38_gwas_cro_franke_4_19_1.tab', header = T)
class(cro)
dim(cro)

##Subset to MH plot##
library(qqman)
cro_sub <- subset(cro, select = c( 'Chr' , 'Marker', 'Position', 'PValue' ))
colnames(cro_sub) <- c('CHR', 'SNP', 'BP', 'P')
head(cro_sub, n =10)
class(cro_sub)
#Delete row CHR X
cro_sub <- cro_sub[-c(941368), ]

manhattan(cro_sub) ##Error in manhattan(cro_sub) : 
                   ##CHR column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again.
qq(gwasResults$P)

#############################################################################

atd <- read.table('hg38_gwas_ic_atd_cooper_4_19_1.tab', header = T)

cel <- read.table('hg38_gwas_ic_cel_trynka_4_19_1.tab', header = T)

jia <- read.table('hg38_gwas_ic_jia_hinks_UK_4_19_1.tab', header = T)

ms <- read.table('hg38_gwas_ic_ms_imsgc_4_19_2.tab', header = T)

nar <- read.table('hg38_gwas_ic_nar_faraco_4_19_1.tab', header = T)

pbc <- read.table('hg38_gwas_ic_pbc_liu_4_19_1.tab', header = T)

pso <- read.table('hg38_gwas_ic_pso_tsoi_4_19_1.tab', header = T)

ra <- read.table('hg38_gwas_ic_ra_eyre_4_19_1.tab', header = T)

t1d <- read.table('hg38_gwas_ic_t1d_onengut_cc_4_19_1.tab', header = T)

t1d_meta <- read.table('summary_statistics/hg38_gwas_ic_t1d_onengut_meta_4_19_1.tab', header = T)

sle <- read.table('hg38_gwas_sle_bentham_4_20_0.tab', header = T)

uc <- read.table('hg38_gwas_uc_anderson_4_19_2.tab', header = T)

#################Associate miRSNPs intersected in Haploreg###############
mirsnps_haplo <- read.delim("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/mirsnps_shared_haploreg.txt", header = T)
colnames(mirsnps_haplo)
length(unique(mirsnps_haplo$rsID))

mirsnps_df <- subset(mirsnps_haplo, select = c("chr", "pos_hg38", "rsID" ,"query_snp_rsid","GENCODE_name", "eQTL"))

###intersect with Navy Bayes###########################
inter_atd <- intersect(mirsnps_haplo$rsID,navyD$SNP)
length(inter_mirsnps)
inter_mirsnps <- as.data.frame(inter_mirsnps)
colnames(inter_mirsnps) = 'SNP'

mirsnps_navyscore <- (merge(inter_mirsnps, navyD, by = 'SNP'))
colnames(mirsnps_navyscore) = c("rsID", "miRNA", "combined_score")
colnames(mirsnps_navyscore)
mirsnps_navy_haplo <- merge(mirsnps_navyscore, mirsnps_haplo, by = 'rsID')
length(unique(mirsnps_navy_haplo$rsID))
colnames(mirsnps_navy_haplo)

mirsnps_df_navy <- subset(mirsnps_navy_haplo, select = c("chr", "pos_hg38", "rsID" ,"combined_score", "miRNA", "GENCODE_name","ref", "alt", "eQTL"))

save(mirsnps_df_navy, file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/mirsnps_navyscore_eqtls.RData")
write.table(mirsnps_df_navy, "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/mirsnps_navyscore_eqtls.txt", sep = '\t',row.names = FALSE, quote = FALSE)




#############Intersect with miRSNP Navy Bayes ##############################
load(file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Navy_miRSNP_data/navyResultsData.RData")
colnames(navyD)
length(unique(navyD$SNP))
length(unique(navyD$miRNA))

inter_ige <- intersect(ige_df$rsID,navyD$SNP)

alop_navyscore <- (merge(inter_alop, navyD, by = 'SNP'))
colnames(alop_navyscore) = c("rsID", "miRNA", "combined_score")

crohn_navy_haplo <- merge(crohn_navyscore, crohn_df, by = 'rsID')