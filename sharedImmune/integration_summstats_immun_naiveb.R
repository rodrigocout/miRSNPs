####Analysis miRSNP Naive Bayes and Autoimmune########
setwd('/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/summary_statistics')

##Summary statistics from Immunobase###########
cro <- read.table("hg38_gwas_cro_franke_4_19_1.tab",  header = T)
class(cro)
dim(cro)
colnames(cro)
head(cro, n =10)
length(unique(cro$Marker))

atd <- read.table('hg38_gwas_ic_atd_cooper_4_19_1.tab', header = T)
colnames(atd)
length(unique(atd$Marker))

cel <- read.table('hg38_gwas_ic_cel_trynka_4_19_1.tab', header = T)
length(unique(cel$Marker))

jia <- read.table('hg38_gwas_ic_jia_hinks_UK_4_19_1.tab', header = T)
length(unique(jia$Marker))

nar <- read.table('hg38_gwas_ic_nar_faraco_4_19_1.tab', header = T)
length(unique(nar$Marker))

pbc <- read.table('hg38_gwas_ic_pbc_liu_4_19_1.tab', header = T)
length(unique(pbc$Marker))

pso <- read.table('hg38_gwas_ic_pso_tsoi_4_19_1.tab', header = T)
length(unique(pso$Marker))

ra <- read.table('hg38_gwas_ic_ra_eyre_4_19_1.tab', header = T)
length(unique(ra$Marker))

t1d <- read.table('hg38_gwas_ic_t1d_onengut_cc_4_19_1.tab', header = T)
length(unique(t1d$Marker))

t1d_meta <- read.table('hg38_gwas_ic_t1d_onengut_meta_4_19_1.tab', header = T)
length(unique(t1d_meta$Marker))

sle <- read.table('hg38_gwas_sle_bentham_4_20_0.tab', header = T)
length(unique(sle$Marker))

uc <- read.table('hg38_gwas_uc_anderson_4_19_2.tab', header = T)
length(unique(uc$Marker))

ic_cro <- read.table('gwas_ichip_meta_release.txt', header = T)
dim(ic_cro)
colnames(ic_cro)
head(ic_cro, n = 10)
length(unique(ic_cro$ICHIP_SNP))

cro_meta <- read.table('cd-meta.txt', header = T)
colnames(cro_meta)
length(unique(cro_meta$SNP))

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



#############Intersect ICHIP with miRSNP Navy Bayes ##############################
load(file = "/Users/lgmh/Documents/rodrigo/IdeaS/Navy_miRSNP_data/navyResultsData.RData")
colnames(navyD)
length(unique(navyD$SNP))
length(unique(navyD$miRNA))

unionSeveral <- function(...) { Reduce(union, list(...)) } 
ichip <- unionSeveral(cro$Marker, atd$Marker, cel$Marker, 
                 jia$Marker, pbc$Marker, nar$Marker,
                 pbc$Marker,pso$Marker,ra$Marker,
                 sle$Marker,t1d_meta$Marker, uc$Marker)
length(ichip) 

###Or #######################################################
ichip_union <- Reduce(union, list(cro$Marker, atd$Marker, cel$Marker, 
                                      jia$Marker, pbc$Marker, nar$Marker,
                                      pbc$Marker,pso$Marker,ra$Marker,
                                      sle$Marker,t1d_meta$Marker, uc$Marker))

length(unique(ichip_union))


mirsnps_ichip <- intersect(ichip_union, navyD$SNP)
length(unique(mirsnps_ichip))

###To intersect per each disease and extract P-values ##############
##CD with mirSNPs########################
inter_cro <- intersect(cro$Marker,navyD$SNP)
length(unique(inter_cro))

inter_cd <- as.data.frame(inter_cro)
colnames(inter_cd) = 'Marker'
#colnames(inter_cd) = 'SNP'

cd_mirsnps <- (merge(inter_cd, cro, by = 'Marker'))
head(cd_mirsnps)

cd_mirsnps <- cd_mirsnps[order(cd_mirsnps$PValue),]
head(cd_mirsnps_order, n = 20)

save(cd_mirsnps, file = "/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/tables/cd_mirsnps.RData")
write.table()

cd <- cd_mirsnps[cd_mirsnps$PValue < 1.0e-04,]
class(cd)
dim(cd)
head(cd)
colnames(cd)

cd <- cd[order(cd$PValue),]
head(cd)
dim(cd)
colnames(cd) <- c('SNP', "Chr", "Position",        
                  "PValue", "OR.MinAllele.", "LowerOR",        
                  "UpperOR","Alleles.Maj.Min")

crohn_navyscore <- (merge(cd, navyD, by = 'SNP'))

save(cd, file = "/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/tables/cron_mirsnps_all.RData")
save(crohn_navyscore, file = "/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/tables/cron_mirsnps_signf.RData")



##############################################################
##ATD and miRNAS
inter_atd <- intersect(atd$Marker, navyD$SNP)
length(unique(inter_atd))

inter_atd <- as.data.frame(inter_atd)
colnames(inter_atd) = 'Marker'

atd_mirsnps <- (merge(inter_atd, atd, by = 'Marker'))
head(atd_mirsnps)

##ordering by p-values
atd_mirsnps <- atd_mirsnps[order(atd_mirsnps$PValue),]
head(atd_mirsnps)

write.table(atd_mirsnps, "/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/tables/atd_mirsnps.txt", sep = '\t', quote = F, row.names = F)
save(atd_mirsnps, file = "/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/tables/atd_mirsnps.Rdata")

atd <- atd_mirsnps[atd_mirsnps$PValue < 1.0e-04,]
class(atd)
dim(atd)
head(atd)
colnames(atd)

atd <- atd[order(atd$PValue),]
head(atd)
dim(atd)
colnames(atd) <- c('SNP', "Chr", "Position",        
                  "PValue", "OR.MinAllele.", "LowerOR",        
                  "UpperOR","Alleles.Maj.Min")

##For the second analysis
colnames(inter_atd) = 'SNP'

atd_navyscore <- (merge(inter_atd, navyD, by = 'SNP'))

atd_navy <- merge(atd_navyscore, atd, by = 'SNP')
head(atd_navy)
save(atd_navy, file = "/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/tables/atd_navyscore.RData")

#################################################
##Celiac and miRSNPs
inter_cel <- intersect(cel$Marker, navyD$SNP)
length(unique(inter_cel))

inter_cel <- as.data.frame(inter_cel)
colnames(inter_cel) = 'Marker'

cel_mirsnps <- (merge(inter_cel, cel, by = 'Marker'))
head(cel_mirsnps)

##ordering by p-values
cel_mirsnps <- cel_mirsnps[order(cel_mirsnps$PValue),]
head(cel_mirsnps)

write.table(cel_mirsnps,"/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/tables/cel_mirsnps.txt", sep = '\t', quote = F, row.names = F)

##Subset SNPs with suggestive P-values
cel_mirs <- cel_mirsnps[cel_mirsnps$PValue < 1.0e-04,]
class(cel_mirs)
dim(cel_mirs)
head(cel_mirs)

cel_mirs <- cel_mirs[order(cel_mirs$PValue),]
head(cel_mirs)
dim(cel_mirs)
colnames(cel_mirs) <- c('SNP', "Chr", "Position",        
                   "PValue", "OR.MinAllele.", "LowerOR",        
                   "UpperOR","Alleles.Maj.Min")

##For the second analysis
colnames(inter_cel) = 'SNP'

cel_navyscore <- (merge(inter_cel, navyD, by = 'SNP'))

cel_navy <- merge(cel_navyscore, cel_mirs, by = 'SNP')
head(cel_navy)
length(unique(cel_navy$miRNA))

save(cel_navy, file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/tables/cel_mirsnps.RData")

#################################################
##JIA and miRSNPs
inter_jia <- intersect(jia$Marker, navyD$SNP)
length(unique(inter_jia))

inter_jia <- as.data.frame(inter_jia)
colnames(inter_jia) = 'Marker'

jia_mirsnps <- (merge(inter_jia, jia, by = 'Marker'))
head(jia_mirsnps)

##ordering by p-values
jia_mirsnps <- jia_mirsnps[order(jia_mirsnps$PValue),]
head(jia_mirsnps)

write.table(jia_mirsnps,"/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/tables/jia_mirsnps.txt", sep = '\t', quote = F, row.names = F)

##Subset SNPs with suggestive P-values
jia_mirs <- jia_mirsnps[jia_mirsnps$PValue < 1.0e-04,]
class(jia_mirs)
dim(jia_mirs)
head(jia_mirs)

jia_mirs <- jia_mirs[order(jia_mirs$PValue),]
head(jia_mirs)
dim(jia_mirs)
colnames(jia_mirs) <- c('SNP', "Chr", "Position",        
                        "PValue", "OR.MinAllele.", "LowerOR",        
                        "UpperOR","Alleles.Maj.Min")

##For the second analysis
colnames(inter_jia) = 'SNP'

jia_navyscore <- (merge(inter_jia, navyD, by = 'SNP'))

jia_navy <- merge(jia_navyscore, jia_mirs, by = 'SNP')
head(jia_navy)
length(unique(jia_navy$miRNA))

save(jia_navy, file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/tables/jia_mirsnps.RData")

#################################################
##Narcolepsy and miRSNPs
inter_nar <- intersect(nar$Marker, navyD$SNP)
length(unique(inter_nar))

inter_nar <- as.data.frame(inter_nar)
colnames(inter_nar) = 'Marker'

nar_mirsnps <- (merge(inter_nar, nar, by = 'Marker'))
head(nar_mirsnps)

##ordering by p-values
nar_mirsnps <- nar_mirsnps[order(nar_mirsnps$PValue),]
head(nar_mirsnps)

write.table(nar_mirsnps, "/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/tables/nar_mirsnps.txt", sep = '\t', quote = F, row.names = F)

##Subset SNPs with suggestive P-values
nar_mirs <- nar_mirsnps[nar_mirsnps$PValue < 1.0e-04,]
class(nar_mirs)
dim(nar_mirs)
head(nar_mirs)

nar_mirs <- nar_mirs[order(nar_mirs$PValue),]
head(nar_mirs)
dim(nar_mirs)
colnames(nar_mirs) <- c('SNP', "Chr", "Position",        
                        "PValue", "OR.MinAllele.", "LowerOR",        
                        "UpperOR","Alleles.Maj.Min")

##For the second analysis
colnames(inter_nar) = 'SNP'

nar_navyscore <- (merge(inter_nar, navyD, by = 'SNP'))

nar_navy <- merge(nar_navyscore, nar_mirs, by = 'SNP')
head(nar_navy)
length(unique(nar_navy$miRNA))

save(nar_navy, file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/tables/nar_mirsnps.RData")

#################################################
##PBC and miRSNPs
inter_pbc <- intersect(pbc$Marker, navyD$SNP)
length(unique(inter_pbc))

inter_pbc <- as.data.frame(inter_pbc)
colnames(inter_pbc) = 'Marker'

pbc_mirsnps <- (merge(inter_pbc, pbc, by = 'Marker'))
head(pbc_mirsnps)

##ordering by p-values
pbc_mirsnps <- pbc_mirsnps[order(pbc_mirsnps$PValue),]
head(pbc_mirsnps)

write.table(pbc_mirsnps, "/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/tables/pbc_mirsnps.txt", sep = '\t', quote = F, row.names = F)

##Subset SNPs with suggestive P-values
pbc_mirs <- pbc_mirsnps[pbc_mirsnps$PValue < 1.0e-04,]
class(pbc_mirs)
dim(pbc_mirs)
head(pbc_mirs)

pbc_mirs <- pbc_mirs[order(pbc_mirs$PValue),]
head(pbc_mirs)
dim(pbc_mirs)
colnames(pbc_mirs) <- c('SNP', "Chr", "Position",        
                        "PValue", "OR.MinAllele.", "LowerOR",        
                        "UpperOR","Alleles.Maj.Min")

##For the second analysis
colnames(inter_pbc) = 'SNP'

pbc_navyscore <- (merge(inter_pbc, navyD, by = 'SNP'))

pbc_navy <- merge(pbc_navyscore, pbc_mirs, by = 'SNP')
head(pbc_navy)
length(unique(pbc_navy$miRNA))

save(pbc_navy, file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/tables/pbc_mirsnps.RData")


#################################################
##PSO and miRSNPs
inter_pso <- intersect(pso$Marker, navyD$SNP)
length(unique(inter_pso))

inter_pso <- as.data.frame(inter_pso)
colnames(inter_pso) = 'Marker'

pso_mirsnps <- (merge(inter_pso, pso, by = 'Marker'))
head(pso_mirsnps)

##ordering by p-values
pso_mirsnps <- pso_mirsnps[order(pso_mirsnps$PValue),]
head(pso_mirsnps)

write.table(pso_mirsnps, "/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/tables/pso_mirsnps.txt", sep = '\t', quote = F, row.names = F)

##Subset SNPs with suggestive P-values
pso_mirs <- pso_mirsnps[pso_mirsnps$PValue < 1.0e-04,]
class(pso_mirs)

colnames(pso_mirs) <- c('SNP', "Chr", "Position",        
                       "PValue", "OR.MinAllele.", "LowerOR",        
                       "UpperOR","Alleles.Maj.Min")

##For the second analysis
colnames(inter_pso) = 'SNP'

pso_navyscore <- (merge(inter_pso, navyD, by = 'SNP'))

pso_navy <- merge(pso_navyscore, pbc_mirs, by = 'SNP')
head(pso_navy)
length(unique(pso_navy$miRNA))

save(pso_navy, file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/tables/pso_mirsnps.RData")


#################################################
##RA and miRSNPs
inter_ra <- intersect(ra$Marker, navyD$SNP)
length(unique(inter_ra))

inter_ra <- as.data.frame(inter_ra)
colnames(inter_ra) = 'Marker'

ra_mirsnps <- (merge(inter_ra, ra, by = 'Marker'))
head(ra_mirsnps)

##ordering by p-values
ra_mirsnps <- ra_mirsnps[order(ra_mirsnps$PValue),]
head(ra_mirsnps)

write.table(ra_mirsnps, "/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/tables/ra_mirsnps.txt", sep = '\t', quote = F, row.names = F)

##Subset SNPs with suggestive P-values
ra_mirs <- ra_mirsnps[ra_mirsnps$PValue < 1.0e-04,]
class(ra_mirs)
dim(ra_mirs)
head(ra_mirs)

ra_mirs <- ra_mirs[order(ra_mirs$PValue),]
head(ra_mirs)
dim(ra_mirs)
colnames(ra_mirs) <- c('SNP', "Chr", "Position",        
                       "PValue", "OR.MinAllele.", "LowerOR",        
                       "UpperOR","Alleles.Maj.Min")

##For the second analysis
colnames(inter_ra) = 'SNP'

ra_navyscore <- (merge(inter_ra, navyD, by = 'SNP'))

ra_navy <- merge(ra_navyscore, ra_mirs, by = 'SNP')
head(ra_navy)
length(unique(ra_navy$miRNA))

save(ra_navy, file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/tables/ra_mirsnps.RData")

#################################################
##T1D and miRSNPs
inter_t1d <- intersect(t1d$Marker, navyD$SNP)
length(unique(inter_t1d))

inter_t1d <- as.data.frame(inter_t1d)
colnames(inter_t1d) = 'Marker'

t1d_mirsnps <- (merge(inter_t1d, t1d, by = 'Marker'))
head(t1d_mirsnps)

##ordering by p-values
t1d_mirsnps <- t1d_mirsnps[order(t1d_mirsnps$PValue),]
head(t1d_mirsnps)

write.table(t1d_mirsnps, "/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/tables/t1d_mirsnps.txt", sep = '\t', quote = F, row.names = F)


##Subset SNPs with suggestive P-values
t1d_mirs <- t1d_mirsnps[t1d_mirsnps$PValue < 1.0e-04,]
class(t1d_mirs)
dim(t1d_mirs)
head(t1d_mirs)

t1d_mirs <- t1d_mirs[order(t1d_mirs$PValue),]
head(t1d_mirs)
dim(t1d_mirs)
colnames(t1d_mirs) <- c('SNP', "Chr", "Position",        
                        "PValue", "OR.MinAllele.", "LowerOR",        
                        "UpperOR","Alleles.Maj.Min")

##For the second analysis
colnames(inter_t1d) = 'SNP'

t1d_navyscore <- (merge(inter_t1d, navyD, by = 'SNP'))

t1d_navy <- merge(t1d_navyscore, t1d_mirs, by = 'SNP')
head(t1d_navy)
length(unique(t1d_navy$miRNA))

save(t1d_navy, file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/tables/t1d_mirsnps.RData")

#################################################
#################################################
##SLE and miRSNPs
inter_sle <- intersect(sle$Marker, navyD$SNP)
length(unique(inter_sle))

inter_sle <- as.data.frame(inter_sle)
colnames(inter_sle) = 'Marker'

sle_mirsnps <- (merge(inter_sle, sle, by = 'Marker'))
head(sle_mirsnps)

##ordering by p-values
sle_mirsnps <- sle_mirsnps[order(sle_mirsnps$PValue),]
head(sle_mirsnps)

write.table(sle_mirsnps, "/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/tables/sle_mirsnps.txt", sep = '\t', quote = F, row.names = F)

##Subset SNPs with suggestive P-values
sle_mirs <- sle_mirsnps[sle_mirsnps$PValue < 1.0e-04,]
class(sle_mirs)
dim(sle_mirs)
head(sle_mirs)

sle_mirs <- sle_mirs[order(sle_mirs$PValue),]
head(sle_mirs)
dim(sle_mirs)
colnames(sle_mirs) <- c('SNP', "Chr", "Position",        
                       "PValue", "OR.MinAllele.", "LowerOR",        
                       "UpperOR","Alleles.Maj.Min")

##For the second analysis
colnames(inter_sle) = 'SNP'

sle_navyscore <- (merge(inter_sle, navyD, by = 'SNP'))

sle_navy <- merge(sle_navyscore, sle_mirs, by = 'SNP')
head(sle_navy)
length(unique(sle_navy$miRNA))

save(sle_navy, file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/tables/sle_mirsnps.RData")

################################
##UC and miRSNPs
inter_uc <- intersect(uc$Marker, navyD$SNP)
length(unique(inter_uc))

inter_uc <- as.data.frame(inter_uc)
colnames(inter_uc) = 'Marker'

uc_mirsnps <- (merge(inter_uc, uc, by = 'Marker'))
head(uc_mirsnps)

##ordering by p-values
uc_mirsnps <- uc_mirsnps[order(uc_mirsnps$PValue),]
head(uc_mirsnps)

write.table(uc_mirsnps,"/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/tables/uc_mirsnps.txt", sep = '\t', quote = F, row.names = F)

##Subset SNPs with suggestive P-values
uc_mirs <- uc_mirsnps[uc_mirsnps$PValue < 1.0e-04,]
class(uc_mirs)
dim(uc_mirs)
head(uc_mirs)

uc_mirs <- uc_mirs[order(uc_mirs$PValue),]
head(uc_mirs)
dim(uc_mirs)
colnames(uc_mirs) <- c('SNP', "Chr", "Position",        
                       "PValue", "OR.MinAllele.", "LowerOR",        
                       "UpperOR","Alleles.Maj.Min")

##For the second analysis
colnames(inter_uc) = 'SNP'

uc_navyscore <- (merge(inter_uc, navyD, by = 'SNP'))

uc_navy <- merge(uc_navyscore, uc_mirs, by = 'SNP')
head(uc_navy)
length(unique(uc_navy$miRNA))

save(uc_navy, file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/tables/uc_mirsnps.RData")
