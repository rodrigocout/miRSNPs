##Subset Navy Results#####
load(file = "~/Dropbox/BED_files_miRSNPS/databases/navyResultsData.RData")

##Combined score >= 0.8
subset_0.8 <- navyD[navyD$combined_score >= 0.8, ]
subset_0.8 <- subset_0.8[order(subset_0.8$combined_score),]
subset_0.8_SNPs <- unique(subset_0.8$SNP)
length(subset_0.8_SNPs)
write.table(subset_0.8_SNPs, 'navy_mirsnps08.txt', sep = '\t', quote = F, row.names = F)

##Combined score >= 0.9
subset_0.9 <- navyD[navyD$combined_score >= 0.9, ]
subset_0.9 <- subset_0.9[order(subset_0.9$combined_score),]
subset_0.9_SNPs <- unique(subset_0.9$SNP)
length(subset_0.9_SNPs)
write.table(subset_0.9_SNPs, 'navy_mirsnps09.txt', sep = '\t', quote = F, row.names = F)


##Combined score >=0.95
subset_0.95 <- navyD[navyD$combined_score >= 0.95, ]
subset_0.95 <- subset_0.95[order(subset_0.9$combined_score),]
subset_0.95_SNPs <- unique(subset_0.95$SNP)
length(subset_0.95_SNPs)
write.table(subset_0.95_SNPs, 'navy_mirsnps095.txt', sep = '\t', quote = F, row.names = F)

##Combined score >= 0.99
subset_0.99 <- navyD[navyD$combined_score >= 0.99, ]
subset_0.99 <- subset_0.99[order(subset_0.99$combined_score),]
subset_0.99_SNPs <- unique(subset_0.99$SNP)
length(subset_0.99_SNPs)
write.table(subset_0.99_SNPs, 'navy_mirsnps099.txt', sep = '\t', quote = F, row.names = F)

###########Gradient with GWAS results #######################

gwas <- read.delim('~/Dropbox/BED_files_miRSNPS/databases/Table_SNPs_GWAS_BreastCancer_NavyScore_eQTL.txt', header = T)
class(gwas)
colnames(gwas)
hist(gwas$combined_Navy_score)

##Naive Combined score >= 0.8
subset_0.8 <- gwas[gwas$combined_Navy_score >= 0.8, ]
subset_0.8 <- subset_0.8[order(subset_0.8$combined_Navy_score),]

subset_0.8_SNPs <- unique(subset_0.8$Proxy_predicted_miRSNP)
length(subset_0.8_SNPs)

##Combined score >= 0.9
subset_0.9 <- gwas[gwas$combined_Navy_score >= 0.9, ]
subset_0.9 <- subset_0.9[order(subset_0.9$combined_Navy_score),]
length(subset_0.9)
dim(subset_0.9)

subset_0.9_SNPs <- unique(subset_0.9$Proxy_predicted_miRSNP)
length(subset_0.9_SNPs)

##Combined score >=0.95
subset_0.95 <- gwas[gwas$combined_Navy_score >= 0.95, ]
subset_0.95 <- subset_0.95[order(subset_0.9$combined_Navy_score),]
dim(subset_0.95)

subset_0.95_SNPs <- unique(subset_0.95$Proxy_predicted_miRSNP)
length(subset_0.95_SNPs)

##Combined score >= 0.99
subset_0.99 <- gwas[gwas$combined_Navy_score >= 0.99, ]
subset_0.99 <- subset_0.99[order(subset_0.99$combined_Navy_score),]
dim(subset_0.99)

subset_0.99_SNPs <- unique(subset_0.99$Proxy_predicted_miRSNP)
length(subset_0.99_SNPs)
