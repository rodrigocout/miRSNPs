load("/Users/lgmh/Dropbox/BED_files_miRSNPS/databases/navyResultsData.RData")
colnames(navyD)
length(unique(navyD$SNP))
length(unique(navyD$miRNA))

#####Load significant diseases tables############
top_atd <- read.table("/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/sign_mirsnps_atd_annotate.txt.annot", heade = T)
colnames(top_atd)
head(top_atd)

###Merge with significant mirsnps###############
mirsnps_navy_atd <- (merge(top_atd, navyD, by = 'SNP'))
mirsnps_navy_atd <- mirsnps_navy_atd[order(mirsnps_navy_atd$P),]

write.table(mirsnps_navy_atd, '/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/mirsnps_atd_navybayes.txt', sep = '\t', quote = F, row.names = F)

#####Load significant diseases tables############
nonhla_atd <- read.table("/Users/lgmh//Dropbox/shared_immune_mirsnps/tables/top_mirsnps/nonhla_sig_mirsnps_atd_annotate.txt.annot", heade = T)
colnames(nonhla_atd)
head(nonhla_atd)

###Merge with significant mirsnps###############
nonhla_mirsnps_navy_atd <- (merge(nonhla_atd, navyD, by = 'SNP'))
nonhla_mirsnps_navy_atd <- nonhla_mirsnps_navy_atd[order(nonhla_mirsnps_navy_atd$P),]

write.table(nonhla_mirsnps_navy_atd, '/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/nonhla_mirsnps_atd_navybayes.txt', sep = '\t', quote = F, row.names = F)

#####Load significant diseases tables############
nonhla_cd <- read.table("/Users/lgmh//Dropbox/shared_immune_mirsnps/tables/top_mirsnps/nonhla_sig_mirsnps_cd_annotate.txt.annot", header = T)
colnames(nonhla_cd)
head(nonhla_cd)

###Merge with significant mirsnps###############
nonhla_mirsnps_navy_cd <- (merge(nonhla_cd, navyD, by = 'SNP'))
nonhla_mirsnps_navy_cd <- nonhla_mirsnps_navy_cd[order(nonhla_mirsnps_navy_cd$P),]

write.table(nonhla_mirsnps_navy_cd, '/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/nonhla_mirsnps_cd_navybayes.txt', sep = '\t', quote = F, row.names = F)

#####Load significant diseases tables############
nonhla_cel <- read.table("/Users/lgmh//Dropbox/shared_immune_mirsnps/tables/top_mirsnps/nonhla_sig_mirsnps_cel_annotate.txt.annot", header = T)
colnames(nonhla_cel)
head(nonhla_cel)

###Merge with significant mirsnps###############
nonhla_mirsnps_navy_cel <- (merge(nonhla_cel, navyD, by = 'SNP'))
nonhla_mirsnps_navy_cel <- nonhla_mirsnps_navy_cel[order(nonhla_mirsnps_navy_cel$P),]

write.table(nonhla_mirsnps_navy_cel, '/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/nonhla_mirsnps_cel_navybayes.txt', sep = '\t', quote = F, row.names = F)

#####Load significant diseases tables############
nonhla_jia <- read.table("/Users/lgmh//Dropbox/shared_immune_mirsnps/tables/top_mirsnps/nonhla_sig_mirsnps_jia_annotate.txt.annot", header = T)
colnames(nonhla_jia)
head(nonhla_jia)

###Merge with significant mirsnps###############
nonhla_mirsnps_navy_jia <- (merge(nonhla_jia, navyD, by = 'SNP'))
nonhla_mirsnps_navy_jia <- nonhla_mirsnps_navy_jia[order(nonhla_mirsnps_navy_jia$P),]

write.table(nonhla_mirsnps_navy_jia, '/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/nonhla_mirsnps_jia_navybayes.txt', sep = '\t', quote = F, row.names = F)

#####Load significant diseases tables############
nonhla_pbc <- read.table("/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/nonhla_sig_mirsnps_pbc_annotate.txt.annot", header = T)

###Merge with significant mirsnps###############
nonhla_mirsnps_navy_pbc <- (merge(nonhla_pbc, navyD, by = 'SNP'))
nonhla_mirsnps_navy_pbc <- nonhla_mirsnps_navy_pbc[order(nonhla_mirsnps_navy_pbc$P),]

write.table(nonhla_mirsnps_navy_pbc, '/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/nonhla_mirsnps_pbc_navybayes.txt', sep = '\t', quote = F, row.names = F)

#####Load significant diseases tables############
nonhla_pso <- read.table("/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/nonhla_sig_mirsnps_pso_annotate.txt.annot", header = T)

###Merge with significant mirsnps###############
nonhla_mirsnps_navy_pso <- (merge(nonhla_pso, navyD, by = 'SNP'))
nonhla_mirsnps_navy_pso <- nonhla_mirsnps_navy_pso[order(nonhla_mirsnps_navy_pso$P),]

write.table(nonhla_mirsnps_navy_pso, '/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/nonhla_mirsnps_pso_navybayes.txt', sep = '\t', quote = F, row.names = F)


#####Load significant diseases tables############
nonhla_sle <- read.table("/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/nonhla_sig_mirsnps_sle_annotate.txt.annot", header = T)

###Merge with significant mirsnps###############
nonhla_mirsnps_navy_sle <- (merge(nonhla_sle, navyD, by = 'SNP'))
nonhla_mirsnps_navy_sle <- nonhla_mirsnps_navy_sle[order(nonhla_mirsnps_navy_sle$P),]

write.table(nonhla_mirsnps_navy_sle, '/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/nonhla_mirsnps_sle_navybayes.txt', sep = '\t', quote = F, row.names = F)

#####Load significant diseases tables############
nonhla_t1d <- read.table("/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/nonhla_sig_mirsnps_t1d_annotate.txt.annot", header = T)

###Merge with significant mirsnps###############
nonhla_mirsnps_navy_t1d <- (merge(nonhla_t1d, navyD, by = 'SNP'))
nonhla_mirsnps_navy_t1d <- nonhla_mirsnps_navy_t1d[order(nonhla_mirsnps_navy_t1d$P),]

write.table(nonhla_mirsnps_navy_t1d, '/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/nonhla_mirsnps_t1d_navybayes.txt', sep = '\t', quote = F, row.names = F)


#####Load significant diseases tables############
nonhla_uc <- read.table("/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/nonhla_sig_mirsnps_uc_annotate.txt.annot", header = T)

###Merge with significant mirsnps###############
nonhla_mirsnps_navy_uc <- (merge(nonhla_uc, navyD, by = 'SNP'))
nonhla_mirsnps_navy_uc <- nonhla_mirsnps_navy_uc[order(nonhla_mirsnps_navy_uc$P),]

write.table(nonhla_mirsnps_navy_uc, '/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/nonhla_mirsnps_uc_navybayes.txt', sep = '\t', quote = F, row.names = F)

#####Load significant diseases tables############
nonhla_ms <- read.table("/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/nonhla_sig_mirsnps_ms_annotate.txt.annot", header = T)

###Merge with significant mirsnps###############
nonhla_mirsnps_navy_ms <- (merge(nonhla_ms, navyD, by = 'SNP'))
nonhla_mirsnps_navy_ms <- nonhla_mirsnps_navy_ms[order(nonhla_mirsnps_navy_ms$P),]

write.table(nonhla_mirsnps_navy_ms, '/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/nonhla_mirsnps_ms_navybayes.txt', sep = '\t', quote = F, row.names = F)



