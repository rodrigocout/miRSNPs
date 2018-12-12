###Integration with GWAS results#####
###GWAS Results from GWAS Catalog####
gwas <- read.delim('~/Dropbox/BED_files_miRSNPS/databases/GWAS_breast_cancer_Proxies.txt', header = T)
head(gwas)
colnames(gwas)
#I took LD info (r2 >= 0.8 and D' = 1 in CEU 1KGenomas) from SNAP for breast cancer associated SNPs 
length(gwas$Proxies)
gwas_bcancer <- as.data.frame(unique(gwas$Proxies))
colnames(gwas_bcancer) = 'Proxies'

load(file = "~/Dropbox/BED_files_miRSNPS/databases/navyResultsData.RData")
colnames(navyD)
length(navyD$SNP)

#Intersect
inter_gwas <- intersect(gwas_bcancer$Proxies,navyD$SNP)
length(inter_gwas)
inter_gwas <- as.data.frame(inter_gwas)
colnames(inter_gwas) = 'SNP'

save(inter_gwas, file = "~/Dropbox/BED_files_miRSNPS/databases/inter_gwas.RData")

gwas_navyscore <- (merge(inter_gwas, navyD, by = 'SNP'))

gwas_navy_unique <- unique(gwas_navyscore$SNP)

save(gwas_nayvscore, file = "~/Dropbox/BED_files_miRSNPS/databases/gwas_navyscore.RData")
