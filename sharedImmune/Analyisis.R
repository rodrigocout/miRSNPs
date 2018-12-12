####Analysis miRSNP Naive Bayes and Autoimmune########
setwd('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/')

####Extracting shared SNPs in each AID#####################
#####Shared SNPs in AIT ##################################
ait <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/AIT_shared.txt', header  = T)
colnames(ait)
dim(ait)
class(ait)
head(ait, n = 10)
ait_SNPs <- unique(ait$Rs.Id)
ait_SNPs <-as.data.frame(ait_SNPs)
ait_SNPs[,1]
write.table(ait_SNPs[,1], "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/AIT_SNPs.txt", sep ='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

###Integrate AIT shared SNPs with Haploreg############
ait_haploreg <- read.delim("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/AIT_haploreg.txt", header = T)
head(ait_haploreg)
colnames(ait_haploreg)
length(ait_haploreg$eQTL)
length(unique(ait_haploreg$rsID))
ait_proxies <- unique(ait_haploreg$rsID)
class(ait_proxies)

ait_df <- subset(ait_haploreg, select = c("chr", "pos_hg38","query_snp_rsid", "rsID","r2", "ref", "alt","EUR", "gwas", "GENCODE_name", "dbSNP_functional_annotation", "eQTL"))

#####Shared SNPs in Alopecia#############
alop <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Alopecia_shared.txt', header = T)
colnames(alop)
dim(alop)
class(alop)
head(alop, n = 10)
alop_SNPs <- unique(alop$Rs.Id)
alop_SNPs <-as.data.frame(alop_SNPs)
alop_SNPs[,1]
write.table(alop_SNPs[,1], "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Alop_SNPs.txt", sep ='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

############Integrate with Haploreg v4.1 results ###############
alop_haploreg <- read.delim("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Alop_haploreg.txt", header = T)
head(alop_haploreg)
colnames(alop_haploreg)
length(alop_haploreg$eQTL)
length(unique(alop_haploreg$rsID))
alop_proxies <- unique(alop_haploreg$rsID)
length(alop_proxies)

alop_df <- subset(alop_haploreg, select = c("chr", "pos_hg38","query_snp_rsid", "rsID","r2", "ref", "alt","EUR","eQTL", "gwas", "GENCODE_name", "dbSNP_functional_annotation"))

########Shared SNPs in Anky######
ank <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Ankylosing_Spondylitis_shared.txt', header = T)
ank_SNPs <- unique(ank$Rs.Id)
ank_SNPs <-as.data.frame(ank_SNPs)
ank_SNPs[,1]
write.table(ank_SNPs[,1], "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/Ank_SNPs.txt", sep ='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

#####Integrate with Haploreg v4.1 results ##############
ank_haploreg <- read.delim("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Ank_haploreg.txt", header = T)
head(ank_haploreg)
colnames(ank_haploreg)
length(unique(ank_haploreg$rsID))
ank_proxies <- unique(ank_haploreg$rsID)
length(ank_proxies)

ank_df <- subset(ank_haploreg, select = c("chr", "pos_hg38","query_snp_rsid", "rsID","r2", "ref", "alt","EUR","eQTL", "gwas", "GENCODE_name", "dbSNP_functional_annotation"))

#######Shared SNPs in Cel#########
cel <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Celiac_shared.txt', header = T)
cel_SNPs <- unique(cel$Rs.Id)
cel_SNPs <-as.data.frame(cel_SNPs)
cel_SNPs[,1]
length(cel_SNPs[,1])
write.table(cel_SNPs[,1], "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Cel_SNPs.txt", sep ='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

#####Integrate with Haploreg v4.1 results ##############
cel_haploreg <- read.delim("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Cel_haploreg.txt", header = T)
head(cel_haploreg)
colnames(cel_haploreg)
length(unique(cel_haploreg$rsID))
cel_proxies <- unique(cel_haploreg$rsID)
length(cel_proxies)

cel_df <- subset(cel_haploreg, select = c("chr", "pos_hg38","query_snp_rsid", "rsID","r2", "ref", "alt","EUR","eQTL", "gwas", "GENCODE_name", "dbSNP_functional_annotation"))

cel_haplo_sub <- subset(cel_haploreg, select = c("chr", "pos_hg38","query_snp_rsid", "rsID","r2", "ref", "alt","EUR", "GENCODE_name", "eQTL"))
colnames(cel_haplo_sub) <- c("chr", "pos_hg38","query_snp_rsid", "SNP","r2", "ref", "alt","EUR", "GENCODE_name", "eQTL")

########Shared SNPs in Crohn######
crohn <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Crohns_shared.txt', header = T)
crohn_SNPs <- unique(crohn$Rs.Id)
crohn_SNPs <-as.data.frame(crohn_SNPs)
crohn_SNPs[,1]
length(crohn_SNPs[,1])
write.table(crohn_SNPs[,1], "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Crohn_SNPs.txt", sep ='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

#####Integrate with Haploreg v4.1 results ##############
crohn_haploreg <- read.delim("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Crohn_haploreg.txt", header = T)
head(crohn_haploreg)
colnames(crohn_haploreg)
length(unique(crohn_haploreg$rsID))
crohn_proxies <- unique(crohn_haploreg$rsID)
length(crohn_proxies)

crohn_df <- subset(crohn_haploreg, select = c("chr", "pos_hg38","query_snp_rsid", "rsID","r2", "ref", "alt","EUR","eQTL", "gwas", "GENCODE_name", "dbSNP_functional_annotation"))

cro_haplo_sub <- subset(crohn_haploreg, select = c("chr", "pos_hg38","query_snp_rsid", "rsID","r2", "ref", "alt","EUR", "GENCODE_name", "eQTL"))
colnames(cro_haplo_sub) <- c("chr", "pos_hg38","query_snp_rsid", "SNP","r2", "ref", "alt","EUR", "GENCODE_name", "eQTL")


########Shared SNPs in IBD###################
ibd <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/IBD_shared.txt', header = T)
ibd_SNPs <- unique(ibd$Rs.Id)
ibd_SNPs <-as.data.frame(ibd_SNPs)
ibd_SNPs[,1]
length(ibd_SNPs[,1])
write.table(ibd_SNPs[,1], "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/IBD_SNPs.txt", sep ='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

#####Integrate with Haploreg v4.1 results ##############
ibd_haploreg <- read.delim("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/IBD_haploreg.txt", header = T)
colnames(ibd_haploreg)
length(unique(ibd_haploreg$rsID))
ibd_proxies <- unique(ibd_haploreg$rsID)
length(ibd_proxies)

ibd_df <- subset(ibd_haploreg, select = c("chr", "pos_hg38","query_snp_rsid", "rsID","r2", "ref", "alt","EUR","eQTL", "gwas", "GENCODE_name", "dbSNP_functional_annotation"))

##########Shared IGE##############
ige <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/IGE_Allergic.txt', header = T)
ige_SNPs <- unique(ige$Rs.Id)
ige_SNPs <-as.data.frame(ige_SNPs)
ige_SNPs[,1]
length(ige_SNPs[,1])
write.table(ige_SNPs[,1], "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/IGE_SNPs.txt", sep ='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

#####Integrate with Haploreg v4.1 results ##############
ige_haploreg <- read.delim("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/IGE_haploreg.txt", header = T)
colnames(ige_haploreg)
length(unique(ige_haploreg$rsID))
ige_proxies <- unique(ige_haploreg$rsID)
length(ige_proxies)

ige_df <- subset(ige_haploreg, select = c("chr", "pos_hg38","query_snp_rsid", "rsID","r2", "ref", "alt","EUR","eQTL", "gwas", "GENCODE_name", "dbSNP_functional_annotation"))

########Shared Juvenile#############
juv <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Juvenile_Arthistis_shared.txt', header = T)
juv_SNPs <- unique(juv$Rs.Id)
juv_SNPs <-as.data.frame(juv_SNPs)
juv_SNPs[,1]
length(juv_SNPs[,1])
write.table(juv_SNPs[,1], "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Juv_SNPs.txt", sep ='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

#####Integrate with Haploreg v4.1 results ##############
juv_haploreg <- read.delim("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Juv_haploreg.txt", header = T)
colnames(juv_haploreg)
length(unique(juv_haploreg$rsID))
juv_proxies <- unique(juv_haploreg$rsID)
length(juv_proxies)

juv_df <- subset(juv_haploreg, select = c("chr", "pos_hg38","query_snp_rsid", "rsID","r2", "ref", "alt","EUR","eQTL", "gwas", "GENCODE_name", "dbSNP_functional_annotation"))

jia_haplo_sub <- subset(juv_haploreg, select = c("chr", "pos_hg38","query_snp_rsid", "rsID","r2", "ref", "alt","EUR", "GENCODE_name", "eQTL"))
colnames(jia_haplo_sub) <- c("chr", "pos_hg38","query_snp_rsid", "SNP","r2", "ref", "alt","EUR", "GENCODE_name", "eQTL")

########Shared Lupus##############
lups <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Lupus_shared.txt' , header = T)
lups_SNPs <- unique(lups$Rs.Id)
lups_SNPs <-as.data.frame(lups_SNPs)
lups_SNPs[,1]
length(lups_SNPs[,1])
write.table(lups_SNPs[,1], "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Lupus_SNPs.txt", sep ='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

#####Integrate with Haploreg v4.1 results ##############
lups_haploreg <- read.delim("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Lupus_haploreg.txt", header = T)
colnames(lups_haploreg)
length(unique(lups_haploreg$rsID))
lups_proxies <- unique(lups_haploreg$rsID)
length(lups_proxies)

lups_df <- subset(lups_haploreg, select = c("chr", "pos_hg38","query_snp_rsid", "rsID","r2", "ref", "alt","EUR","eQTL", "gwas", "GENCODE_name", "dbSNP_functional_annotation"))

sle_haplo_sub <- subset(lups_haploreg, select = c("chr", "pos_hg38","query_snp_rsid", "rsID","r2", "ref", "alt","EUR", "GENCODE_name", "eQTL"))
colnames(sle_haplo_sub) <- c("chr", "pos_hg38","query_snp_rsid", "SNP","r2", "ref", "alt","EUR", "GENCODE_name", "eQTL")


##########Shared MS#####################
ms <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/MS_shared.txt', header = T)
ms_SNPs <- unique(ms$Rs.Id)
ms_SNPs <-as.data.frame(ms_SNPs)
ms_SNPs[,1]
length(ms_SNPs[,1])
write.table(ms_SNPs[,1], "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/MS_SNPs.txt", sep ='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

#####Integrate with Haploreg v4.1 results ##############
ms_haploreg <- read.delim("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/MS_haploreg.txt", header = T)
colnames(ms_haploreg)
length(unique(ms_haploreg$rsID))
ms_proxies <- unique(ms_haploreg$rsID)
length(ms_proxies)

ms_df <- subset(ms_haploreg, select = c("chr", "pos_hg38","query_snp_rsid", "rsID","r2", "ref", "alt","EUR","eQTL", "gwas", "GENCODE_name", "dbSNP_functional_annotation"))

ms_haplo_sub <- subset(ms_haploreg, select = c("chr", "pos_hg38","query_snp_rsid", "rsID","r2", "ref", "alt","EUR", "GENCODE_name", "eQTL"))
colnames(ms_haplo_sub) <- c("chr", "pos_hg38","query_snp_rsid", "SNP","r2", "ref", "alt","EUR", "GENCODE_name", "eQTL")

##########Shared Narcolepsia###############
narc <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Narcolepsy_shared.txt', header = T)
narc_SNPs <- unique(narc$Rs.Id)
narc_SNPs <-as.data.frame(narc_SNPs)
narc_SNPs[,1]
length(narc_SNPs[,1])
write.table(narc_SNPs[,1], "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Narc_SNPs.txt", sep ='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

#####Integrate with Haploreg v4.1 results ##############
narc_haploreg <- read.delim("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Narc_haploreg.txt", header = T)
colnames(narc_haploreg)
length(unique(narc_haploreg$rsID))
narc_proxies <- unique(narc_haploreg$rsID)
length(narc_proxies)

narc_df <- subset(narc_haploreg, select = c("chr", "pos_hg38","query_snp_rsid", "rsID","r2", "ref", "alt","EUR","eQTL", "gwas", "GENCODE_name", "dbSNP_functional_annotation"))

nar_haplo_sub <- subset(narc_haploreg, select = c("chr", "pos_hg38","query_snp_rsid", "rsID","r2", "ref", "alt","EUR", "GENCODE_name", "eQTL"))
colnames(nar_haplo_sub) <- c("chr", "pos_hg38","query_snp_rsid", "SNP","r2", "ref", "alt","EUR", "GENCODE_name", "eQTL")

########Shared Primary Biliary############
prim <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Primary_Biliary_shared.txt', header =T)
prim_SNPs <- unique(prim$Rs.Id)
prim_SNPs <-as.data.frame(prim_SNPs)
prim_SNPs[,1]
length(prim_SNPs[,1])
write.table(prim_SNPs[,1], "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Prim_SNPs.txt", sep ='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

#####Integrate with Haploreg v4.1 results ##############
prim_haploreg <- read.delim("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Prim_haploreg.txt", header = T)
colnames(prim_haploreg)
length(unique(prim_haploreg$rsID))
prim_proxies <- unique(prim_haploreg$rsID)
length(prim_proxies)

prim_df <- subset(prim_haploreg, select = c("chr", "pos_hg38","query_snp_rsid", "rsID","r2", "ref", "alt","EUR","eQTL", "gwas", "GENCODE_name", "dbSNP_functional_annotation"))

pbc_haplo_sub <- subset(prim_haploreg, select = c("chr", "pos_hg38","query_snp_rsid", "rsID","r2", "ref", "alt","EUR", "GENCODE_name", "eQTL"))
colnames(pbc_haplo_sub) <- c("chr", "pos_hg38","query_snp_rsid", "SNP","r2", "ref", "alt","EUR", "GENCODE_name", "eQTL")

#######Shared Psoriasis#################
psori <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Psoriasis_shared.txt', header = T)
psori_SNPs <- unique(psori$Rs.Id)
psori_SNPs <-as.data.frame(psori_SNPs)
psori_SNPs[,1]
length(psori_SNPs[,1])
write.table(psori_SNPs[,1], "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Psori_SNPs.txt", sep ='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

#####Integrate with Haploreg v4.1 results ##############
psori_haploreg <- read.delim("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Psori_haploreg.txt", header = T)
colnames(psori_haploreg)
length(unique(psori_haploreg$rsID))
psori_proxies <- unique(psori_haploreg$rsID)
length(psori_proxies)

psori_df <- subset(psori_haploreg, select = c("chr", "pos_hg38","query_snp_rsid", "rsID","r2", "ref", "alt","EUR","eQTL", "gwas", "GENCODE_name", "dbSNP_functional_annotation"))

pso_haplo_sub <- subset(psori_haploreg, select = c("chr", "pos_hg38","query_snp_rsid", "rsID","r2", "ref", "alt","EUR", "GENCODE_name", "eQTL"))
colnames(pso_haplo_sub) <- c("chr", "pos_hg38","query_snp_rsid", "SNP","r2", "ref", "alt","EUR", "GENCODE_name", "eQTL")

##########Shared RA####################
ra <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/RA_shared.txt', header = T)
ra_SNPs <- unique(ra$Rs.Id)
ra_SNPs <-as.data.frame(ra_SNPs)
ra_SNPs[,1]
length(ra_SNPs[,1])
write.table(ra_SNPs[,1], "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/RA_SNPs.txt", sep ='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

#####Integrate with Haploreg v4.1 results ##############
ra_haploreg <- read.delim("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/RA_haploreg.txt", header = T)
colnames(ra_haploreg)
length(unique(ra_haploreg$rsID))
ra_proxies <- unique(ra_haploreg$rsID)
length(ra_proxies)

ra_df <- subset(ra_haploreg, select = c("chr", "pos_hg38","query_snp_rsid", "rsID","r2", "ref", "alt","EUR","eQTL", "gwas", "GENCODE_name", "dbSNP_functional_annotation"))

ra_haplo_sub <- subset(ra_haploreg, select = c("chr", "pos_hg38","query_snp_rsid", "rsID","r2", "ref", "alt","EUR", "GENCODE_name", "eQTL"))
colnames(ra_haplo_sub) <- c("chr", "pos_hg38","query_snp_rsid", "SNP","r2", "ref", "alt","EUR", "GENCODE_name", "eQTL")

#######Shared SCLE #####################
scle <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Scleroderma_shared.txt', header = T)
scle_SNPs <- unique(scle$Rs.Id)
scle_SNPs <-as.data.frame(scle_SNPs)
scle_SNPs[,1]
length(scle_SNPs[,1])
write.table(scle_SNPs[,1], "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Sclero_SNPs.txt", sep ='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

#####Integrate with Haploreg v4.1 results ##############
scle_haploreg <- read.delim("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Sclero_haploreg.txt", header = T)
colnames(scle_haploreg)
length(unique(scle_haploreg$rsID))
scle_proxies <- unique(scle_haploreg$rsID)
length(scle_proxies)

scle_df <- subset(scle_haploreg, select = c("chr", "pos_hg38","query_snp_rsid", "rsID","r2", "ref", "alt","EUR","eQTL", "gwas", "GENCODE_name", "dbSNP_functional_annotation"))

##########Shared Sjogren ######################
sjogren <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Sjogren_SD_shared.txt', header = T)
sj_SNPs <- unique(sjogren$Rs.Id)
sj_SNPs <-as.data.frame(sj_SNPs)
sj_SNPs[,1]
length(sj_SNPs[,1])
write.table(sj_SNPs[,1], "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Sjogren_SNPs.txt", sep ='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

#####Integrate with Haploreg v4.1 results ##############
sj_haploreg <- read.delim("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Sjogren_haploreg.txt", header = T)
colnames(sj_haploreg)
length(unique(sj_haploreg$rsID))
sj_proxies <- unique(sj_haploreg$rsID)
length(sj_proxies)

sj_df <- subset(sj_haploreg, select = c("chr", "pos_hg38","query_snp_rsid", "rsID","r2", "ref", "alt","EUR","eQTL", "gwas", "GENCODE_name", "dbSNP_functional_annotation"))

###########Shared T1D #############################
t1d <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/T1D_shared.txt', header = T)
t1d_SNPs <- unique(t1d$Rs.Id)
t1d_SNPs <-as.data.frame(t1d_SNPs)
t1d_SNPs[,1]
length(t1d_SNPs[,1])
write.table(t1d_SNPs[,1], "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/T1D_SNPs.txt", sep ='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

#####Integrate with Haploreg v4.1 results ##############
t1d_haploreg <- read.delim("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/T1D_haploreg.txt", header = T)
colnames(t1d_haploreg)
length(unique(t1d_haploreg$rsID))
t1d_proxies <- unique(t1d_haploreg$rsID)
length(t1d_proxies)

t1d_df <- subset(t1d_haploreg, select = c("chr", "pos_hg38","query_snp_rsid", "rsID","r2", "ref", "alt","EUR","eQTL", "gwas", "GENCODE_name", "dbSNP_functional_annotation"))

t1d_haplo_sub <- subset(t1d_haploreg, select = c("chr", "pos_hg38","query_snp_rsid", "rsID","r2", "ref", "alt","EUR", "GENCODE_name", "eQTL"))
colnames(t1d_haplo_sub) <- c("chr", "pos_hg38","query_snp_rsid", "SNP","r2", "ref", "alt","EUR", "GENCODE_name", "eQTL")

##########Shared UC ############################
uc <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/UC_shared.txt', header =T)
uc_SNPs <- unique(uc$Rs.Id)
uc_SNPs <-as.data.frame(uc_SNPs)
uc_SNPs[,1]
length(uc_SNPs[,1])
write.table(uc_SNPs[,1], "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/UC_SNPs.txt", sep ='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

#####Integrate with Haploreg v4.1 results ##############
uc_haploreg <- read.delim("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/UC_haploreg.txt", header = T)
colnames(uc_haploreg)
length(unique(uc_haploreg$rsID))
uc_proxies <- unique(uc_haploreg$rsID)
length(uc_proxies)

uc_df <- subset(uc_haploreg, select = c("chr", "pos_hg38","query_snp_rsid", "rsID","r2", "ref", "alt","EUR","eQTL", "gwas", "GENCODE_name", "dbSNP_functional_annotation"))

uc_haplo_sub <- subset(uc_haploreg, select = c("chr", "pos_hg38","query_snp_rsid", "rsID","r2", "ref", "alt","EUR", "GENCODE_name", "eQTL"))
colnames(uc_haplo_sub) <- c("chr", "pos_hg38","query_snp_rsid", "SNP","r2", "ref", "alt","EUR", "GENCODE_name", "eQTL")

#########Shared Vitiligo ###################
vit <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Vitiligo_shared.txt', header = T)
vit_SNPs <- unique(vit$Rs.Id)
vit_SNPs <-as.data.frame(vit_SNPs)
vit_SNPs[,1]
length(vit_SNPs[,1])
write.table(vit_SNPs[,1], "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Vit_SNPs.txt", sep ='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

#####Integrate with Haploreg v4.1 results ##############
vit_haploreg <- read.delim("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/Vit_haploreg.txt", header = T)
colnames(vit_haploreg)
length(unique(vit_haploreg$rsID))
vit_proxies <- unique(vit_haploreg$rsID)
length(vit_proxies)

vit_df <- subset(vit_haploreg, select = c("chr", "pos_hg38","query_snp_rsid", "rsID","r2", "ref", "alt","EUR","eQTL", "gwas", "GENCODE_name", "dbSNP_functional_annotation"))

##########Shared SNPs with proxies#############
snps_proxies <- c(as.character(alop_haploreg$rsID), as.character(ait_haploreg$rsID), as.character(ank_haploreg$rsID), 
                  as.character(cel_haploreg$rsID), as.character(crohn_haploreg$rsID), as.character(ibd_haploreg$rsID),
                  as.character(ige_haploreg$rsID), as.character(juv_haploreg$rsID), as.character(lups_haploreg$rsID),
                  as.character(ms_haploreg$rsID), as.character(narc_haploreg$rsID), as.character(prim_haploreg$rsID),
                  as.character(psori_haploreg$rsID), as.character(ra_haploreg$rsID), as.character(scle_haploreg$rsID),
                  as.character(sj_haploreg$rsID), as.character(t1d_haploreg$rsID), as.character(uc_haploreg$rsID),
                  as.character(vit_haploreg$rsID))
length(unique(snps_proxies))
write.table(snps_proxies, 'SNPs_proxies_shared_autoimmune.txt', sep = '\t', quote = F, row.names = F)

#######Summary statistics from Immunobase###########
cro <- read.table('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/summary_statistics/hg38_gwas_cro_franke_4_19_1.tab', header = T)
class(cro)
dim(cro)
colnames(cro)
head(cro, n =10)
colnames(cro) <- c("SNP", "Chr" , "Position", "PValue",          
                   "OR.MinAllele.", "LowerOR", "UpperOR", "Alleles.Maj.Min.")

atd <- read.table('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/summary_statistics/hg38_gwas_ic_atd_cooper_4_19_1.tab', header = T)

cel <- read.table('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/summary_statistics/hg38_gwas_ic_cel_trynka_4_19_1.tab', header = T)
colnames(cel) <- c("SNP", "Chr" , "Position", "PValue",          
                   "OR.MinAllele.", "LowerOR", "UpperOR", "Alleles.Maj.Min.")


jia <- read.table('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/summary_statistics/hg38_gwas_ic_jia_hinks_UK_4_19_1.tab', header = T)
colnames(jia) <- c("SNP", "Chr" , "Position", "PValue",          
                   "OR.MinAllele.", "LowerOR", "UpperOR", "Alleles.Maj.Min.")

ms <- read.table('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/summary_statistics/hg38_gwas_ic_ms_imsgc_4_19_2.tab', header = T)
colnames(ms) <- c("SNP", "Chr" , "Position", "PValue",          
                  "OR.MinAllele.", "LowerOR", "UpperOR", "Alleles.Maj.Min.")

nar <- read.table('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/summary_statistics/hg38_gwas_ic_nar_faraco_4_19_1.tab', header = T)
colnames(nar) <- c("SNP", "Chr" , "Position", "PValue",          
                   "OR.MinAllele.", "LowerOR", "UpperOR", "Alleles.Maj.Min.")

pbc <- read.table('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/summary_statistics/hg38_gwas_ic_pbc_liu_4_19_1.tab', header = T)
colnames(pbc) <- c("SNP", "Chr" , "Position", "PValue",          
                   "OR.MinAllele.", "LowerOR", "UpperOR", "Alleles.Maj.Min.")

pso <- read.table('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/summary_statistics/hg38_gwas_ic_pso_tsoi_4_19_1.tab', header = T)
colnames(pso) <- c("SNP", "Chr" , "Position", "PValue",          
                   "OR.MinAllele.", "LowerOR", "UpperOR", "Alleles.Maj.Min.")

ra <- read.table('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/summary_statistics/hg38_gwas_ic_ra_eyre_4_19_1.tab', header = T)
colnames(ra) <- c("SNP", "Chr" , "Position", "PValue",          
                  "OR.MinAllele.", "LowerOR", "UpperOR", "Alleles.Maj.Min.")

t1d <- read.table('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/summary_statistics/hg38_gwas_ic_t1d_onengut_cc_4_19_1.tab', header = T)
colnames(t1d) <- c("SNP", "Chr" , "Position", "PValue",          
                   "OR.MinAllele.", "LowerOR", "UpperOR", "Alleles.Maj.Min.")

t1d_meta <- read.table('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/summary_statistics/hg38_gwas_ic_t1d_onengut_meta_4_19_1.tab', header = T)
colnames(t1d_meta) <- c("SNP", "Chr" , "Position", "PValue",          
                        "OR.MinAllele.", "LowerOR", "UpperOR", "Alleles.Maj.Min.")

sle <- read.table('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/summary_statistics/hg38_gwas_sle_bentham_4_20_0.tab', header = T)
colnames(sle) <- c("SNP", "Chr" , "Position", "PValue",          
                   "OR.MinAllele.", "LowerOR", "UpperOR", "Alleles.Maj.Min.")

uc <- read.table('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/summary_statistics/hg38_gwas_uc_anderson_4_19_2.tab', header = T)
colnames(uc) <- c("SNP", "Chr" , "Position", "PValue",          
                  "OR.MinAllele.", "LowerOR", "UpperOR", "Alleles.Maj.Min.")

#############Intersect with miRSNP Navy Bayes ##############################
load(file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Navy_miRSNP_data/navyResultsData.RData")
colnames(navyD)
length(unique(navyD$SNP))
length(unique(navyD$miRNA))

##AIT###
colnames(ait_df)
inter_ait <- intersect(ait_df$rsID,navyD$SNP)
length(inter_ait)
inter_ait

##Alopecia###
inter_alop <- intersect(alop_df$rsID,navyD$SNP)
length(inter_ait)
inter_alop
inter_alop <- as.data.frame(inter_alop)
colnames(inter_alop) = 'SNP'

alop_navyscore <- (merge(inter_alop, navyD, by = 'SNP'))
colnames(alop_navyscore) = c("rsID", "miRNA", "combined_score")
colnames(alop_navyscore)
alop_navy_haplo <- merge(alop_navyscore, alop_df, by = 'rsID')
length(unique(alop_navy_haplo$rsID))

save(alop_navyscore, file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/alop_navyscore.RData")
write.table(alop_navy_haplo, "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/alop_navyscore.txt", sep = '\t',row.names = FALSE, quote = FALSE)


#######Anky######
inter_anky <- intersect(ank_df$rsID,navyD$SNP)
length(inter_anky)
inter_anky
inter_anky <- as.data.frame(inter_anky)
colnames(inter_anky) = 'SNP'

anky_navyscore <- (merge(inter_anky, navyD, by = 'SNP'))

save(anky_navyscore, file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/anky_navyscore.RData")
anky_navyscore
colnames(anky_navyscore) = c("rsID", "miRNA", "combined_score")
colnames(anky_navyscore)
anky_navy_haplo <- merge(anky_navyscore, ank_df, by = 'rsID')
length(unique(anky_navy_haplo$rsID))

write.table(anky_navy_haplo, "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/anky_navyscore.txt", sep = '\t',row.names = FALSE, quote = FALSE)

######Cel##############################
inter_cel <- intersect(cel_df$rsID,navyD$SNP)
length(inter_cel)
inter_cel
inter_cel <- as.data.frame(inter_cel)
colnames(inter_cel) = 'SNP'

cel_navyscore <- (merge(inter_cel, navyD, by = 'SNP'))
colnames(cel_navyscore)

save(cel_navyscore, file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/cel_navyscore.RData")
cel_navyscore
colnames(cel_navyscore) = c("rsID", "miRNA", "combined_score")
colnames(cel_navyscore)
cel_navy_haplo <- merge(cel_navyscore, cel_df, by = 'rsID')
cel_navy_haplo
length(unique(cel_navy_haplo$rsID))

write.table(cel_navy_haplo, "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/cel_navyscore.txt", sep = '\t', row.names = FALSE, quote = FALSE)

cel_immunobase_navy <- merge(cel_navyscore, cel, by = 'SNP')
length(unique(cel_immunobase_navy$SNP))

cel_mirsnps <- merge(cel_immunobase_navy, cel_haplo_sub, by = 'SNP')
write.table(cel_mirsnps, 'celiac_mirsnps.txt', sep = '\t', quote = F, row.names = F)

####Crohn#########################################
inter_crohn <- intersect(crohn_df$rsID,navyD$SNP)
length(inter_crohn)
inter_crohn
inter_crohn <- as.data.frame(inter_crohn)
colnames(inter_crohn) = 'SNP'

crohn_navyscore <- (merge(inter_crohn, navyD, by = 'SNP'))

save(crohn_navyscore, file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/crohn_navyscore.RData")
colnames(crohn_navyscore) = c("rsID", "miRNA", "combined_score")
colnames(crohn_navyscore)
crohn_navy_haplo <- merge(crohn_navyscore, crohn_df, by = 'rsID')
length(unique(crohn_navy_haplo$rsID))

write.table(crohn_navy_haplo, "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/crohn_navyscore.txt", sep = '\t', row.names = FALSE, quote = FALSE)

cro_immunobase_navy <- merge(crohn_navyscore, cro, by = 'SNP')
cro_mirsnps <- merge(cro_immunobase_navy, cro_haplo_sub, by = 'SNP')
write.table(cro_mirsnps, 'cro_mirsnps.txt', sep = '\t', quote = F, row.names = F)

#####IBD##########
inter_ibd <- intersect(ibd_df$rsID,navyD$SNP)
length(inter_ibd)
inter_ibd
inter_ibd <- as.data.frame(inter_ibd)
colnames(inter_ibd) = 'SNP'

ibd_navyscore <- (merge(inter_ibd, navyD, by = 'SNP'))

save(ibd_navyscore, file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/ibd_navyscore.RData")

colnames(ibd_navyscore) = c("rsID", "miRNA", "combined_score")
colnames(ibd_navyscore)
ibd_navy_haplo <- merge(ibd_navyscore, ibd_df, by = 'rsID')
length(unique(ibd_navy_haplo$rsID))

write.table(ibd_navy_haplo, "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/ibd_navyscore.txt", sep = '\t', row.names = FALSE, quote = FALSE)

####IGE#########
inter_ige <- intersect(ige_df$rsID,navyD$SNP)
length(inter_ige)
inter_ige
inter_ige <- as.data.frame(inter_ige)
colnames(inter_ige) = 'SNP'

ige_navyscore <- (merge(inter_ige, navyD, by = 'SNP'))

save(ige_navyscore, file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/ige_navyscore.RData")

colnames(ige_navyscore) = c("rsID", "miRNA", "combined_score")
colnames(ige_navyscore)
ige_navy_haplo <- merge(ige_navyscore, ige_df, by = 'rsID')
length(unique(ige_navy_haplo$rsID))

write.table(ige_navy_haplo, "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/ige_navyscore.txt", sep = '\t', row.names = FALSE, quote = FALSE)

##########Juv########
inter_juv <- intersect(juv_df$rsID,navyD$SNP)
length(inter_juv)
inter_juv
inter_juv <- as.data.frame(inter_juv)
colnames(inter_juv) = 'SNP'

juv_navyscore <- (merge(inter_juv, navyD, by = 'SNP'))

save(juv_navyscore, file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/juv_navyscore.RData")

colnames(juv_navyscore) = c("rsID", "miRNA", "combined_score")
colnames(juv_navyscore)
juv_navy_haplo <- merge(juv_navyscore, juv_df, by = 'rsID')
length(unique(juv_navy_haplo$rsID))

write.table(juv_navy_haplo, "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/juv_navyscore.txt", sep = '\t', row.names = FALSE, quote = FALSE)

jia_immunobase_navy <- merge(juv_navyscore, jia, by = 'SNP')

jia_mirsnps <- merge(jia_immunobase_navy, jia_haplo_sub, by = 'SNP')
write.table(jia_mirsnps, 'jia_mirsnps.txt', sep = '\t', quote = F, row.names = F)

########Lupus######
inter_lupus <- intersect(lups_df$rsID,navyD$SNP)
length(inter_lupus)
inter_lupus
inter_lupus <- as.data.frame(inter_lupus)
colnames(inter_lupus) = 'SNP'

lupus_navyscore <- (merge(inter_lupus, navyD, by = 'SNP'))

save(lupus_navyscore, file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/lupus_navyscore.RData")
colnames(lupus_navyscore) = c("rsID", "miRNA", "combined_score")
colnames(lupus_navyscore)
lup_navy_haplo <- merge(lupus_navyscore, lups_df, by = 'rsID')
length(unique(lup_navy_haplo$rsID))

write.table(lup_navy_haplo, "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/lupus_navyscore.txt", sep = '\t', row.names = FALSE, quote = FALSE)

sle_immunobase_navy <- merge(lupus_navyscore, sle, by = 'SNP')
sle_mirsnps <- merge(sle_immunobase_navy, sle_haplo_sub, by = 'SNP')
write.table(sle_mirsnps, 'sle_mirsnps.txt', sep = '\t', quote = F, row.names = F)

#######MS###########
inter_ms <- intersect(ms_df$rsID,navyD$SNP)
length(inter_ms)
inter_ms
inter_ms <- as.data.frame(inter_ms)
colnames(inter_ms) = 'SNP'

ms_navyscore <- (merge(inter_ms, navyD, by = 'SNP'))

save(ms_navyscore, file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/ms_navyscore.RData")

colnames(ms_navyscore) = c("rsID", "miRNA", "combined_score")
colnames(ms_navyscore)
ms_navy_haplo <- merge(ms_navyscore, ms_df, by = 'rsID')
length(unique(ms_navy_haplo$rsID))
write.table(ms_navy_haplo, "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/ms_navyscore.txt", sep = '\t', row.names = FALSE, quote = FALSE)

ms_immunobase_navy <- merge(ms_navyscore, ms, by = 'SNP')
ms_mirsnps <- merge(ms_immunobase_navy, ms_haplo_sub, by = 'SNP')
write.table(ms_mirsnps, 'ms_mirsnps.txt', sep = '\t', quote = F, row.names = F)

####Narcolepsy#######
inter_narc <- intersect(narc_df$rsID,navyD$SNP)
length(inter_narc)
inter_narc

########Primary B#########
inter_prim <- intersect(prim_df$rsID,navyD$SNP)
length(inter_prim)
inter_prim
inter_prim <- as.data.frame(inter_prim)
colnames(inter_prim) = 'SNP'

prim_navyscore <- (merge(inter_prim, navyD, by = 'SNP'))

save(prim_navyscore, file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/prim_navyscore.RData")
colnames(prim_navyscore) = c("rsID", "miRNA", "combined_score")
colnames(prim_navyscore)
prim_navy_haplo <- merge(prim_navyscore, prim_df, by = 'rsID')
length(unique(prim_navy_haplo$rsID))
write.table(prim_navy_haplo, "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/prim_navyscore.txt", sep = '\t', row.names = FALSE, quote = FALSE)

pbc_immunobase_navy <- merge(prim_navyscore, pbc, by = 'SNP')
pbc_mirsnps <- merge(pbc_immunobase_navy, pbc_haplo_sub, by = 'SNP')
write.table(pbc_mirsnps, 'pbc_mirsnps.txt', sep = '\t', quote = F, row.names = F)

###Psoriasis#######
inter_psori <- intersect(psori_df$rsID,navyD$SNP)
length(inter_psori)
inter_psori
inter_psori <- as.data.frame(inter_psori)
colnames(inter_psori) = 'SNP'

psori_navyscore <- (merge(inter_psori, navyD, by = 'SNP'))

save(psori_navyscore, file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/psori_navyscore.RData")
colnames(psori_navyscore) = c("rsID", "miRNA", "combined_score")
colnames(psori_navyscore)
psori_navy_haplo <- merge(psori_navyscore, psori_df, by = 'rsID')
length(unique(psori_navy_haplo$rsID))
write.table(psori_navy_haplo, "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/psori_navyscore.txt", sep = '\t', row.names = FALSE, quote = FALSE)

pso_immunobase_navy <- merge(psori_navyscore, pso, by = 'SNP')
pso_mirsnps <- merge(pso_immunobase_navy, pso_haplo_sub, by = 'SNP')
write.table(pso_mirsnps, 'pso_mirsnps.txt', sep = '\t', quote = F, row.names = F)

####RA#######
inter_ra <- intersect(ra_df$rsID,navyD$SNP)
length(inter_ra)
inter_ra
inter_ra <- as.data.frame(inter_ra)
colnames(inter_ra) = 'SNP'

ra_navyscore <- (merge(inter_ra, navyD, by = 'SNP'))

save(ra_navyscore, file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/ra_navyscore.RData")
colnames(ra_navyscore) = c("rsID", "miRNA", "combined_score")
colnames(ra_navyscore)
ra_navy_haplo <- merge(ra_navyscore, ra_df, by = 'rsID')
length(unique(ra_navy_haplo$rsID))
write.table(ra_navy_haplo, "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/ra_navyscore.txt", sep = '\t', row.names = FALSE, quote = FALSE)

ra_immunobase_navy <- merge(ra_navyscore, ra, by = 'SNP')
ra_mirsnps <- merge(ra_immunobase_navy, ra_haplo_sub, by = 'SNP')
write.table(ra_mirsnps, 'ra_mirsnps.txt', sep = '\t', quote = F, row.names = F)

#####Scleroderma######
inter_scle <- intersect(scle_df$rsID,navyD$SNP)
length(inter_scle)
inter_scle

####Sjogren######
inter_sj <- intersect(sj_df$rsID,navyD$SNP)
length(inter_sj)
inter_sj

###T1D#######
inter_t1d <- intersect(t1d_df$rsID,navyD$SNP)
length(inter_t1d)
inter_t1d
inter_t1d <- as.data.frame(inter_t1d)
colnames(inter_t1d) = 'SNP'

t1d_navyscore <- (merge(inter_t1d, navyD, by = 'SNP'))

save(t1d_navyscore, file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/t1d_navyscore.RData")
colnames(t1d_navyscore) = c("rsID", "miRNA", "combined_score")
colnames(t1d_navyscore)
t1d_navy_haplo <- merge(t1d_navyscore, t1d_df, by = 'rsID')
length(unique(t1d_navy_haplo$rsID))
write.table(t1d_navy_haplo, "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/t1d_navyscore.txt", sep = '\t', row.names = FALSE, quote = FALSE)

t1d_immunobase_navy <- merge(t1d_navyscore, t1d, by = 'SNP')
t1d_mirsnps <- merge(t1d_immunobase_navy, t1d_haplo_sub, by = 'SNP')
write.table(t1d_mirsnps, 't1d_mirsnps.txt', sep = '\t', quote = F, row.names = F)

t1d_meta_immunobase_navy <- merge(t1d_navyscore, t1d_meta, by = 'SNP')
t1d_meta_mirsnps <- merge(t1d_meta_immunobase_navy, t1d_haplo_sub, by = 'SNP')
write.table(t1d_meta_mirsnps, 't1d_meta_mirsnps.txt', sep = '\t', quote = F, row.names = F)

####UC#######
inter_uc <- intersect(uc_df$rsID,navyD$SNP)
length(inter_uc)
inter_uc
inter_uc <- as.data.frame(inter_uc)
colnames(inter_uc) = 'SNP'

uc_navyscore <- (merge(inter_uc, navyD, by = 'SNP'))

save(uc_navyscore, file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/uc_navyscore.RData")

colnames(uc_navyscore) = c("rsID", "miRNA", "combined_score")
colnames(uc_navyscore)
uc_navy_haplo <- merge(uc_navyscore, uc_df, by = 'rsID')
length(unique(uc_navy_haplo$rsID))
write.table(uc_navy_haplo, "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/uc_navyscore.txt", sep = '\t', row.names = FALSE, quote = FALSE)

uc_immunobase_navy <- merge(uc_navyscore, uc, by = 'SNP')
uc_mirsnps <- merge(uc_immunobase_navy, uc_haplo_sub, by = 'SNP')
write.table(uc_mirsnps, 'uc_mirsnps.txt', sep = '\t', quote = F, row.names = F)

####Vitiligo####
inter_vit <- intersect(vit_df$rsID,navyD$SNP)
length(inter_vit)
inter_vit
inter_vit <- as.data.frame(inter_vit)
colnames(inter_vit) = 'SNP'

vit_navyscore <- (merge(inter_vit, navyD, by = 'SNP'))

save(vit_navyscore, file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/vit_navyscore.RData")
colnames(vit_navyscore) = c("rsID", "miRNA", "combined_score")
colnames(vit_navyscore)
vit_navy_haplo <- merge(vit_navyscore, vit_df, by = 'rsID')
length(unique(vit_navy_haplo$rsID))
write.table(vit_navy_haplo, "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/vit_navyscore.txt", sep = '\t', row.names = FALSE, quote = FALSE)

#################Associate miRSNPs intersected in Haploreg###############
mirsnps_haplo <- read.delim("mirsnps_shared_haploreg.txt", header = T)
colnames(mirsnps_haplo)
length(unique(mirsnps_haplo$rsID))

mirsnps_df <- subset(mirsnps_haplo, select = c("chr", "pos_hg38", "rsID" ,"query_snp_rsid","GENCODE_name"))
write.table(mirsnps_df, 'mirsnp_shared_haploreg_genesannot.txt', sep = '\t', quote = F, row.names = F)
head(mirsnps_df, n = 20)

###intersect with Navy Bayes####
inter_mirsnps <- intersect(mirsnps_haplo$rsID,navyD$SNP)
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


###########Make plots#############################
###Piechart with number of SNPs per disease########

slices <- c(length(unique(alop_navy_haplo$rsID)), length(unique(anky_navy_haplo$rsID)),
            length(unique(crohn_navy_haplo$rsID)),length(unique(ibd_navy_haplo$rsID)),
            length(unique(uc_navy_haplo$rsID)),length(unique(cel_navy_haplo$rsID)),
            length(unique(t1d_navy_haplo$rsID)),length(unique(ige_navy_haplo$rsID)),
            length(unique(juv_navy_haplo$rsID)),length(unique(lup_navy_haplo$rsID)),
            length(unique(ms_navy_haplo$rsID)),length(unique(prim_navy_haplo$rsID)),
            length(unique(psori_navy_haplo$rsID)),length(unique(ra_navy_haplo$rsID)),
            length(unique(vit_navy_haplo$rsID)))


lbls <- c("Alopecia","Anky","Crohn", "IBD", "UC","Celiac","T1D","IgE", 
          "Juvenile", "Lupus","MS", "Primary B", "Psoriasis", "RA",  "Vitiligo")

cols <- c("deepskyblue","burlywood4","black","darkorange","darkgreen","grey70","grey50", 
          "orange", "blue", "green", "pink", "yellow", "red", "darkorchid4", "white")


pct <- round(slices/sum(slices)*100)
pielabels<- paste(pct, "%", sep="")
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 


pdf(file="/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/figs/Navy_miRSNPs_AutoImmune.pdf", height=8, width=12)
pie(slices,labels = pielabels, main=" Distribution of miRSNPs associated to Autoimmune Diseases",col = cols,
    cex = 0.8, density = 200)
legend("topright", c("Alopecia Areata","Ankylosing Spondylitis","Crohn's disease", "IBD", "Ulcerative Colitis","Celiac Disease","Type 1 Diabetes","IgE", 
                     "Juvenile Arthritis", "Systemic Lupus Erythematosus ","Multiple Sclerosis", "Primary Biliary", "Psoriasis", "Rheumatoid Arthritis",  "Vitiligo"), cex=0.8, fill=cols)
dev.off()

#pie(slices,labels = lbls, col=rainbow(length(lbls)),
 #   main="miRSNPs in Immune Diseases", cex = 0.8, density = 100) 

#########Circusplot###########
###First put all AID miRSNPs together
aid_list <- paste(unique(alop_navy_haplo$rsID),unique(anky_navy_haplo$rsID),
                  unique(crohn_navy_haplo$rsID),unique(ibd_navy_haplo$rsID),
                  unique(uc_navy_haplo$rsID),unique(cel_navy_haplo$rsID),
                  unique(t1d_navy_haplo$rsID),unique(ige_navy_haplo$rsID),
                  unique(juv_navy_haplo$rsID),unique(lup_navy_haplo$rsID),
                  unique(ms_navy_haplo$rsID),unique(prim_navy_haplo$rsID),
                  unique(psori_navy_haplo$rsID),unique(ra_navy_haplo$rsID),
                    unique(vit_navy_haplo$rsID))
                  
aid_list
length(aid_list)
length(unique(aid_list))

aid_mirsnps <-unique(aid_list)
aid_mirsnps <- as.data.frame(t(aid_mirsnps))
head(aid_mirsnps)
aid_mirsnps <- t(aid_mirsnps)

colnames(aid_mirsnps)
length(unique(aid_mirsnps$aid_mirsnps))



