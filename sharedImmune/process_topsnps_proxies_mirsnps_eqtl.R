####Process Associated SNPs and SNPs in LD###
setwd("/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/proceed_data/list_assoc_snps/proxies_functional_annotation/")


##Load Navy Bayes results#################################
load(file = "/Users/lgmh/Documents/rodrigo/IdeaS/Navy_miRSNP_data/navyResultsData.RData")
colnames(navyD)
length(unique(navyD$SNP))
length(unique(navyD$miRNA))

##Load proxies top SNPs per disease###
#############ATD#############################################
atd_ld <- read.delim("ATD_proxies_haploreg.txt", header = T)
#re-order columns for better visualize
atd_ld <- atd_ld[,c(1:4,33,6,5,7:32)]
colnames(atd_ld)

##Intersect###
inter_prox_atd <- intersect(atd_ld$rsID, navyD$SNP)
length(unique(inter_prox_atd))
###No intersection####

##########CRO##########################
cro_ld <- read.delim("cro_proxies_haploreg.txt", header = T)
#re-order columns for better visualize
cro_ld <- cro_ld[,c(1:4,33,6,5,7:32)]
colnames(cro_ld)

#Intersect with proxies
inter_prox_cro <- intersect(cro_ld$rsID, navyD$SNP)
length(unique(inter_prox_cro))

##Merge intersection and mirnas scores
#subset only mirSNPs intersected
cro_sub <- subset(cro_ld[which(cro_ld$rsID %in% inter_prox_cro),])
length(unique(cro_sub$rsID))
colnames(cro_sub)
#Check 
all(cro_sub$rsID %in% inter_prox_cro)

###process ld dataframe
cro_sub <- subset(cro_sub, select = c("chr", "pos_hg38", "query_snp_rsid","rsID" , "r2",  "D.", "ref", "alt", "EUR", "GENCODE_name", "eQTL"))

colnames(cro_sub) <- c("chr", "pos_hg38", "query_snp_rsid","SNP" , "r2",  "D.","ref", "alt", "MAF","GENCODE_name", "eQTL")
head(cro_sub)

cro_df_merg <- merge(cro_sub, navyD, by = "SNP")
cro_df_merg <- subset(cro_df_merg, select = c("chr", "pos_hg38", "query_snp_rsid","SNP" , "r2",  "D.","ref", "alt", "MAF","GENCODE_name","miRNA", "combined_score",  "eQTL"))

colnames(cro_df_merg) <- c("chr", "pos_hg38", "Top_SNP","miRSNP" , "r2",  "D.","ref", "alt", "MAF","GENCODE_name","miRNA", "combined_score",  "eQTL")

###Wrinting table###
write.table (cro_df_merg, "cro_mirsnps_proxies_scores_eqtl.txt", sep = '\t', quote = F, row.names = F)

####Check which disease overlapps with miRSNPs
assoc_cro <- read.delim("/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/proceed_data/top_associated_snps_perdisease/CRO_assoc_variants.txt", header = T, stringsAsFactors = F)
colnames(assoc_cro)
head(assoc_cro)
###Subset
assoc_cro <- subset(assoc_cro, select = c("Rs.Id", "Disease", "P.Value", "Odds.Ratio"))
#Rename to merge
colnames(assoc_cro) <- c("Top_SNP", "Other Disease", "P-value", "OR")

####Merge
cro_res <- merge(assoc_cro, cro_df_merg, by = "Top_SNP") 
colnames(cro_res)
cro_res <- subset(cro_res, select = c("miRSNP","chr", "pos_hg38", "Top_SNP","r2" ,"D.", "Other Disease" , "P-value", "OR",
                                      "ref", "alt", "MAF", "GENCODE_name", "miRNA", "combined_score", "eQTL" )) 
  
#Writing table
write.table(cro_res, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/results/cron_results.txt", sep = '\t', quote = F, row.names = F)
##################################################

####CEL############################
cel_ld <- read.delim("cel_proxies_haploreg.txt", header = T)
#re-order columns for better visualize
cel_ld <- cel_ld[,c(1:4,33,6,5,7:32)]
colnames(cel_ld)
##Col query_snp_rsid is the associated SNP
##Colr rsID is the LD SNP

#Intersect Navy Bayes with proxies
inter_prox_cel <- intersect(cel_ld$rsID, navyD$SNP)
length(unique(inter_prox_cel))

##Merge intersection and mirnas scores
#subset only mirSNPs intersected
cel_sub <- subset(cel_ld[which(cel_ld$rsID %in% inter_prox_cel),])
length(unique(cel_sub$rsID))
colnames(cel_sub)
#Check 
all(cel_sub$rsID %in% inter_prox_cel)

###process ld dataframe
cel_df <- subset(cel_sub, select = c("chr", "pos_hg38", "query_snp_rsid","rsID" , "r2",  "D.", "ref", "alt", "EUR", "GENCODE_name", "eQTL"))

colnames(cel_df) <- c("chr", "pos_hg38", "query_snp_rsid","SNP" , "r2",  "D.","ref", "alt", "MAF","GENCODE_name", "eQTL")
head(cel_df)

cel_df_merg <- merge(cel_df, navyD, by = "SNP")
colnames(cel_df_merg)

cel_df_merg <- subset(cel_df_merg, select = c("chr", "pos_hg38", "query_snp_rsid","SNP" , "r2",  "D.","ref", "alt", "MAF","GENCODE_name","miRNA", "combined_score",  "eQTL"))

colnames(cel_df_merg) <- c("chr", "pos_hg38", "Top_SNP","miRSNP" , "r2",  "D.","ref", "alt", "MAF","GENCODE_name","miRNA", "combined_score",  "eQTL")


####Check which disease overlapps with miRSNPs
assoc_cel <- read.delim("/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/proceed_data/top_associated_snps_perdisease/CEL_assoc_variants.txt", header = T, stringsAsFactors = F)
colnames(assoc_cel)
head(assoc_cel)
###Subset
assoc_cel <- subset(assoc_cel, select = c("Rs.Id", "Disease", "P.Value", "Odds.Ratio"))
#Rename to merge
colnames(assoc_cel) <- c("Top_SNP", "Other Disease", "P-value", "OR")

####Merge
cel_res <- merge(assoc_cel, cel_df_merg, by = "Top_SNP") 


cel_res <- subset(cel_res, select = c("miRSNP","chr", "pos_hg38", "Top_SNP","r2" ,"D.", "Other Disease" , "P-value", "OR",
                                      "ref", "alt", "MAF", "GENCODE_name", "miRNA", "combined_score", "eQTL" )) 
##Ordering by p-value
cel_res <- cel_res[order(cel_res$`P-value`),]
#Writing table
write.table(cel_res, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/results/cel_results.txt", sep = '\t', quote = F, row.names = F)

###################################################################
####jia###########################

jia_ld <- read.delim("jia_proxies_haploreg.txt", header = T)
#re-order columns for better visualize
jia_ld <- jia_ld[,c(1:4,33,6,5,7:32)]
colnames(jia_ld)
##Col query_snp_rsid is the associated SNP
##Colr rsID is the LD SNP

#With proxies
inter_prox_jia <- intersect(jia_ld$rsID, navyD$SNP)
length(unique(inter_prox_jia))

##Merge intersections and mirnas scores
#subset only mirSNPs intersected
jia_sub <- subset(jia_ld[which(jia_ld$rsID %in% inter_prox_jia),])
length(unique(jia_sub$rsID))
colnames(jia_sub)
#Check 
all(jia_sub$rsID %in% inter_prox_jia)

###process ld dataframe
jia_df <- subset(jia_sub, select = c("chr", "pos_hg38", "query_snp_rsid","rsID" , "r2",  "D.", "ref", "alt", "EUR", "GENCODE_name", "eQTL"))

colnames(jia_df) <- c("chr", "pos_hg38", "query_snp_rsid","SNP" , "r2",  "D.","ref", "alt", "MAF","GENCODE_name", "eQTL")
head(jia_df)

jia_df_merg <- merge(jia_df, navyD, by = "SNP")
colnames(jia_df_merg)

jia_df_merg <- subset(jia_df_merg, select = c("chr", "pos_hg38", "query_snp_rsid","SNP" , "r2",  "D.","ref", "alt", "MAF","GENCODE_name","miRNA", "combined_score",  "eQTL"))

colnames(jia_df_merg) <- c("chr", "pos_hg38", "Top_SNP","miRSNP" , "r2",  "D.","ref", "alt", "MAF","GENCODE_name","miRNA", "combined_score",  "eQTL")

####Check which disease overlapps with miRSNPs
assoc_jia <- read.delim("/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/proceed_data/top_associated_snps_perdisease/JIA_assoc_variants.txt", header = T, stringsAsFactors = F)
colnames(assoc_jia)
head(assoc_jia)
###Subset
assoc_jia <- subset(assoc_jia, select = c("Rs.Id", "Disease", "P.Value", "Odds.Ratio"))
#Rename to merge
colnames(assoc_jia) <- c("Top_SNP", "Other Disease", "P-Value", "OR")

####Merge
jia_res <- merge(assoc_jia, jia_df_merg, by = "Top_SNP") 

colnames(jia_res)

jia_res <- subset(jia_res, select = c("miRSNP","chr", "pos_hg38", "Top_SNP","r2" ,"D.", "Other Disease" , "P-Value", "OR",
                                      "ref", "alt", "MAF", "GENCODE_name", "miRNA", "combined_score", "eQTL" )) 

##Ordering by p-value
jia_res <- jia_res[order(jia_res$`P-Value`),]

#Writing table
write.table(jia_res, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/results/jia_results.txt", sep = '\t', quote = F, row.names = F)

#########################################################################
#######MS##################

ms_ld <- read.delim("ms_proxies_haploreg.txt", header = T)
#re-order columns for better visualize
ms_ld <- ms_ld[,c(1:4,33,6,5,7:32)]
colnames(ms_ld)
##Col query_snp_rsid is the associated SNP
##Colr rsID is the LD SNP

##Intersect with miRSNPs with proxies

inter_prox_ms <- intersect(ms_ld$rsID, navyD$SNP)
length(unique(inter_prox_ms))

##Merge intersection with p-values and mirnas scores
#subset only mirSNPs intersected
ms_sub <- subset(ms_ld[which(ms_ld$rsID %in% inter_prox_ms),])
length(unique(ms_sub$rsID))
colnames(ms_sub)
#Check 
all(ms_sub$rsID %in% inter_prox_ms)

###process ld dataframe
ms_df <- subset(ms_sub, select = c("chr", "pos_hg38", "query_snp_rsid","rsID" , "r2",  "D.", "ref", "alt", "EUR", "GENCODE_name", "eQTL"))

colnames(ms_df) <- c("chr", "pos_hg38", "query_snp_rsid","SNP" , "r2",  "D.","ref", "alt", "MAF","GENCODE_name", "eQTL")
head(ms_df)

ms_df_merg <- merge(ms_df, navyD, by = "SNP")
colnames(ms_df_merg)

####Change column name 
colnames(ms_df_merg) <- c("miRSNP" ,"chr", "pos_hg38", "Top_SNP", "r2",  "D.","ref", "alt", "MAF","GENCODE_name","miRNA", "combined_score",  "eQTL")

####Check which disease overlapps with miRSNPs
assoc_ms <- read.delim("/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/proceed_data/top_associated_snps_perdisease/MS_assoc_variants.txt", header = T, stringsAsFactors = F)
colnames(assoc_ms)
head(assoc_ms)
###Subset
assoc_ms <- subset(assoc_ms, select = c("Rs.Id", "Disease", "P.Value", "Odds.Ratio"))
#Rename to merge
colnames(assoc_ms) <- c("Top_SNP", "Other Disease", "P-value", "OR")

####Merge
ms_res <- merge(assoc_ms, ms_df_merg, by = "Top_SNP") 

ms_res <- subset(ms_res, select = c("miRSNP","chr", "pos_hg38", "Top_SNP","r2" ,"D.", "Other Disease" , "P-value", "OR",
                                      "ref", "alt", "MAF", "GENCODE_name", "miRNA", "combined_score", "eQTL" )) 

##Ordering by p-value
ms_res <- ms_res[order(ms_res$`P-value`),]

#Writing table
write.table(ms_res, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/results/ms_results.txt", sep = '\t', quote = F, row.names = F)

###################################################################
####Narcolepsy##############################

nar_ld <- read.delim("nar_proxies_haploreg.txt", header = T)
#re-order columns for better visualize
nar_ld <- nar_ld[,c(1:4,33,6,5,7:32)]
colnames(nar_ld)
##Col query_snp_rsid is the associated SNP
##Colr rsID is the LD SNP

#With proxies
inter_prox_nar <- intersect(nar_ld$rsID, navyD$SNP)
length(unique(inter_prox_nar))
#No intersection with Narcolepsy

###################################################################
######PBC########################################################

pbc_ld <- read.delim("pbc_proxies_haploreg.txt", header = T)
#re-order columns for better visualize
pbc_ld <- pbc_ld[,c(1:4,33,6,5,7:32)]
colnames(pbc_ld)
##Col query_snp_rsid is the associated SNP
##Colr rsID is the LD SNP

#With proxies
inter_prox_pbc <- intersect(pbc_ld$rsID, navyD$SNP)
length(unique(inter_prox_pbc))

##Merge intersection with p-values and mirnas scores
#subset only mirSNPs intersected
pbc_sub <- subset(pbc_ld[which(pbc_ld$rsID %in% inter_prox_pbc),])
length(unique(pbc_sub$rsID))
colnames(pbc_sub)
#Check 
all(pbc_sub$rsID %in% inter_prox_pbc)

###process ld dataframe
pbc_df <- subset(pbc_sub, select = c("chr", "pos_hg38", "query_snp_rsid","rsID" , "r2",  "D.", "ref", "alt", "EUR", "GENCODE_name", "eQTL"))

colnames(pbc_df) <- c("chr", "pos_hg38", "query_snp_rsid","SNP" , "r2",  "D.","ref", "alt", "MAF","GENCODE_name", "eQTL")
head(pbc_df)

pbc_df_merg <- merge(pbc_df, navyD, by = "SNP")
colnames(pbc_df_merg) <- c("miRSNP" ,"chr", "pos_hg38", "Top_SNP", "r2",  "D.","ref", "alt", "MAF","GENCODE_name","miRNA", "combined_score",  "eQTL")

####Check which disease overlapps with miRSNPs
assoc_pbc <- read.delim("/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/proceed_data/top_associated_snps_perdisease/PBC_assoc_variants.txt", header = T, stringsAsFactors = F)
colnames(assoc_pbc)
head(assoc_pbc)
###Subset
assoc_pbc <- subset(assoc_pbc, select = c("Rs.Id", "Disease", "P.Value", "Odds.Ratio"))
#Rename to merge
colnames(assoc_pbc) <- c("Top_SNP", "Other Disease", "P-value", "OR")

####Merge
pbc_res <- merge(assoc_pbc, pbc_df_merg, by = "Top_SNP") 


pbc_res <- subset(pbc_res, select = c("miRSNP","chr", "pos_hg38", "Top_SNP","r2" ,"D.", "Other Disease" , "P-value", "OR",
                                    "ref", "alt", "MAF", "GENCODE_name", "miRNA", "combined_score", "eQTL" )) 

##Ordering by p-value
pbc_res <- pbc_res[order(pbc_res$`P-value`),]
View(pbc_res)

#Writing table
write.table(pbc_res, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/results/pbc_results.txt", sep = '\t', quote = F, row.names = F)

########################################################################
######AS##############################################

as_ld <- read.delim("AS_proxies_haploreg.txt", header = T)
#re-order columns for better visualize
as_ld <- as_ld[,c(1:4,33,6,5,7:32)]
colnames(as_ld)

#With proxies
inter_prox_as <- intersect(as_ld$rsID, navyD$SNP)
length(unique(inter_prox_as))

##Merge intersection and mirnas scores
#subset only mirSNPs intersected
as_sub <- subset(as_ld[which(as_ld$rsID %in% inter_prox_as),])
length(unique(as_sub$rsID))
colnames(as_sub)
#Check 
all(as_sub$rsID %in% inter_prox_as)

###process ld dataframe
as_df <- subset(as_sub, select = c("chr", "pos_hg38", "query_snp_rsid","rsID" , "r2",  "D.", "ref", "alt", "EUR", "GENCODE_name", "eQTL"))

colnames(as_df) <- c("chr", "pos_hg38", "query_snp_rsid","SNP" , "r2",  "D.","ref", "alt", "MAF","GENCODE_name", "eQTL")
head(as_df)

as_df_merg <- merge(as_df, navyD, by = "SNP")
colnames(as_df_merg)

as_df_merg <- as_df_merg[,c(1:10,12,13,11)]
####Change column name 
colnames(as_df_merg) <- c("miRSNP", "chr", "pos_hg38", "Top_SNP","r2","D_prime","ref", "alt",
                           "MAF", "GENCODE_name", "miRNA", "combined_score", "eQTL")


####Check which disease overlapps with miRSNPs
assoc_as <- read.delim("/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/proceed_data/top_associated_snps_perdisease/AS_assoc_variants.txt", header = T, stringsAsFactors = F)
colnames(assoc_as)
head(assoc_as)
###Subset
assoc_as <- subset(assoc_as, select = c("Rs.Id", "Disease", "P.Value", "Odds.Ratio"))
#Rename to merge
colnames(assoc_as) <- c("Top_SNP", "Other Disease", "P-value", "OR")

####Merge
as_res <- merge(assoc_as, as_df_merg, by = "Top_SNP") 
colnames(as_res)

as_res <- as_res[,c(5:7, 1, 8, 9, 2:4,10:16)]


##Ordering by p-value
as_res <- as_res[order(as_res$`P-value`),]

#Writing table
write.table(as_res, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/results/as_results.txt", sep = '\t', quote = F, row.names = F)

##################################################
##################PSO############################

pso_ld <- read.delim("PSO_proxies_haploreg.txt", header = T)
#re-order columns for better visualize
pso_ld <- pso_ld[,c(1:4,33,6,5,7:32)]
colnames(pso_ld)
##Col query_snp_rsid is the psosociated SNP
##Colr rsID is the LD SNP

##Intersect 
#With proxies
inter_prox_pso <- intersect(pso_ld$rsID, navyD$SNP)
length(unique(inter_prox_pso))

##Merge intersection with p-values and mirnas scores
#subset only mirSNPs intersected
pso_sub <- subset(pso_ld[which(pso_ld$rsID %in% inter_prox_pso),])
length(unique(pso_sub$rsID))
colnames(pso_sub)
#Check 
all(pso_sub$rsID %in% inter_prox_pso)

###process ld dataframe
pso_df <- subset(pso_sub, select = c("chr", "pos_hg38", "query_snp_rsid","rsID" , "r2",  "D.", "ref", "alt", "EUR", "GENCODE_name", "eQTL"))

colnames(pso_df) <- c("chr", "pos_hg38", "query_snp_rsid","SNP" , "r2",  "D.","ref", "alt", "MAF","GENCODE_name", "eQTL")
head(pso_df)

pso_df_merg <- merge(pso_df, navyD, by = "SNP")
colnames(pso_df_merg)
###Re-ordering columns e subseting
pso_df_merg <- subset(pso_df_merg, select = c( "SNP","chr", "pos_hg38", "query_snp_rsid", "r2","D.","ref", "alt", "MAF",
                                               "GENCODE_name", "miRNA", "combined_score",  "eQTL"))

####Change column name 
colnames(pso_df_merg) <- c( "miRSNP","chr", "pos_hg38", "Top_SNP", "r2","D.","ref", "alt", "MAF",
                            "GENCODE_name", "miRNA", "combined_score",  "eQTL")


####Check which disease overlapps with miRSNPs
assoc_pso <- read.delim("/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/proceed_data/top_associated_snps_perdisease/PSO_assoc_variants.txt", header = T, stringsAsFactors = F)
colnames(assoc_pso)
head(assoc_pso)
###Subset
assoc_pso <- subset(assoc_pso, select = c("Rs.Id", "Disease", "P.Value", "Odds.Ratio"))
#Rename to merge
colnames(assoc_pso) <- c("Top_SNP", "Other Disease", "P.Value", "OR")

####Merge
pso_res <- merge(assoc_pso, pso_df_merg, by = "Top_SNP") 


##Ordering by p-value
pso_res <- pso_res[order(pso_res$P.Value),]

##Changing order
pso_res <- pso_res[,c(5:7, 1, 8, 9, 2:4,10:16)]


#Writing table
write.table(pso_res, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/results/pso_results.txt", sep = '\t', quote = F, row.names = F)

#########################################################
#################RA#####################################
ra_ld <- read.delim("RA_proxies_haploreg.txt", header = T)
#re-order columns for better visualize
ra_ld <- ra_ld[,c(1:4,33,6,5,7:32)]
colnames(ra_ld)
##Col query_snp_rsid is the associated SNP
##Colr rsID is the LD SNP

##Intersect 
#With proxies
inter_prox_ra <- intersect(ra_ld$rsID, navyD$SNP)
length(unique(inter_prox_ra))

##Merge intersection with p-values and mirnas scores
#subset only mirSNPs intersected
ra_sub <- subset(ra_ld[which(ra_ld$rsID %in% inter_prox_ra),])
length(unique(ra_sub$rsID))
colnames(ra_sub)
#Check 
all(ra_sub$rsID %in% inter_prox_ra)

###process ld dataframe
ra_df <- subset(ra_sub, select = c("chr", "pos_hg38", "query_snp_rsid","rsID" , "r2",  "D.", "ref", "alt", "EUR", "GENCODE_name", "eQTL"))

colnames(ra_df) <- c("chr", "pos_hg38", "query_snp_rsid","SNP" , "r2",  "D.","ref", "alt", "MAF","GENCODE_name", "eQTL")
head(ra_df)

ra_df_merg <- merge(ra_df, navyD, by = "SNP")
colnames(ra_df_merg)
###Re-ordering columns e subseting
ra_df_merg <- subset(ra_df_merg, select = c( "SNP","chr", "pos_hg38", "query_snp_rsid", "r2","D.","ref", "alt",
                                             "MAF", "GENCODE_name", "miRNA", "combined_score", "eQTL"))

####Change column name 
colnames(ra_df_merg) <- c( "miRSNP","chr", "pos_hg38", "Top_SNP", "r2","D_prime","ref", "alt",
                           "MAF","GENCODE_name","miRNA", "combined_score", "eQTL")


####Check which disease overlapps with miRSNPs
assoc_ra <- read.delim("/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/proceed_data/top_associated_snps_perdisease/RA_assoc_variants.txt", header = T, stringsAsFactors = F)
colnames(assoc_ra)
head(assoc_ra)
###Subset
assoc_ra <- subset(assoc_ra, select = c("Rs.Id", "Disease", "P.Value", "Odds.Ratio"))
#Rename to merge
colnames(assoc_ra) <- c("Top_SNP", "Other Disease", "P.Value", "OR")

####Merge
ra_res <- merge(assoc_ra, ra_df_merg, by = "Top_SNP") 

##Ordering by p-value
ra_res <- ra_res[order(ra_res$P.Value),]

##Changing order
ra_res <- ra_res[,c(5:7, 1, 8, 9, 2:4,10:16)]

#Writing table
write.table(ra_res, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/results/ra_results.txt", sep = '\t', quote = F, row.names = F)

#####################################################
##########SLE#######################################
sle_ld <- read.delim("SLE_proxies_haploreg.txt", header = T)
#re-order columns for better visualize
sle_ld <- sle_ld[,c(1:4,33,6,5,7:32)]
colnames(sle_ld)
##Col query_snp_rsid is the associated SNP
##Colr rsID is the LD SNP

##Intersect 
#With proxies
inter_prox_sle <- intersect(sle_ld$rsID, navyD$SNP)
length(unique(inter_prox_sle))

##Merge intersection with p-values and mirnas scores
#subset only mirSNPs intersected
sle_sub <- subset(sle_ld[which(sle_ld$rsID %in% inter_prox_sle),])
length(unique(sle_sub$rsID))
colnames(sle_sub)
#Check 
all(sle_sub$rsID %in% inter_prox_sle)

###process ld datafsleme
sle_df <- subset(sle_sub, select = c("chr", "pos_hg38", "query_snp_rsid","rsID" , "r2",  "D.", "ref", "alt", "EUR", "GENCODE_name", "eQTL"))

colnames(sle_df) <- c("chr", "pos_hg38", "query_snp_rsid","SNP" , "r2",  "D.","ref", "alt", "MAF","GENCODE_name", "eQTL")
head(sle_df)

sle_df_merg <- merge(sle_df, navyD, by = "SNP")
colnames(sle_df_merg)
###Re-ordering columns e subseting
sle_df_merg <- subset(sle_df_merg, select = c( "SNP","chr", "pos_hg38", "query_snp_rsid", "r2","D.","ref", "alt",
                                               "MAF", "GENCODE_name","miRNA", "combined_score", "eQTL"))

####Change column name 
colnames(sle_df_merg) <- c( "miRSNP","chr", "pos_hg38", "Top_SNP", "r2","D_prime","ref", "alt",
                            "MAF", "GENCODE_name","miRNA", "combined_score",  "eQTL")


####Check which disease overlapps with miRSNPs
assoc_sle <- read.delim("/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/proceed_data/top_associated_snps_perdisease/SLE_assoc_variants.txt", header = T, stringsAsFactors = F)
colnames(assoc_sle)
head(assoc_sle)
###Subset
assoc_sle <- subset(assoc_sle, select = c("Rs.Id", "Disease", "P.Value", "Odds.Ratio"))
#Rename to merge
colnames(assoc_sle) <- c("Top_SNP", "Disease", "P.Value", "OR")

####Merge
sle_res <- merge(assoc_sle, sle_df_merg, by = "Top_SNP") 

##Changing order
sle_res <- sle_res[,c(5:7, 1, 8, 9, 2:4,10:16)]

##Ordering by p-value
sle_res <- sle_res[order(sle_res$P.Value),]

#Writing table
write.table(sle_res, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/results/sle_results.txt", sep = '\t', quote = F, row.names = F)

##################################################
##T1D#######

t1d_ld <- read.delim("T1D_proxies_haploreg.txt", header = T)
#re-order columns for better visualize
t1d_ld <- t1d_ld[,c(1:4,33,6,5,7:32)]
colnames(t1d_ld)
##Col query_snp_rsid is the associated SNP
##Colr rsID is the LD SNP

##Intersect 
#With proxies
inter_prox_t1d <- intersect(t1d_ld$rsID, navyD$SNP)
length(unique(inter_prox_t1d))

##Merge intersection with p-values and mirnas scores
#subset only mirSNPs intersected
t1d_sub <- subset(t1d_ld[which(t1d_ld$rsID %in% inter_prox_t1d),])
length(unique(t1d_sub$rsID))
colnames(t1d_sub)

#Check 
all(t1d_sub$rsID %in% inter_prox_t1d)

###process ld dataft1dme
t1d_df <- subset(t1d_sub, select = c("chr", "pos_hg38", "query_snp_rsid","rsID" , "r2",  "D.", "ref", "alt", "EUR", "GENCODE_name", "eQTL"))

colnames(t1d_df) <- c("chr", "pos_hg38", "query_snp_rsid","SNP" , "r2",  "D.","ref", "alt", "MAF","GENCODE_name", "eQTL")
head(t1d_df)

t1d_df_merg <- merge(t1d_df, navyD, by = "SNP")
colnames(t1d_df_merg)
###Re-ordering columns e subseting
t1d_df_merg <- subset(t1d_df_merg, select = c("SNP", "chr", "pos_hg38", "query_snp_rsid", "r2","D.","ref", "alt",
                                               "MAF","GENCODE_name","miRNA", "combined_score",  "eQTL"))
##Ordering by P-values
t1d_df_merg <- t1d_df_merg[order(t1d_df_merg$PValue),]

####Change column name 
colnames(t1d_df_merg) <- c( "miRSNP","chr", "pos_hg38", "Top_SNP", "r2","D_prime","ref", "alt",
                            "MAF","GENCODE_name", "miRNA", "combined_score",  "eQTL")

####Check which disease overlapps with miRSNPs
assoc_t1d <- read.delim("/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/proceed_data/top_associated_snps_perdisease/T1D_assoc_variants.txt", header = T, stringsAsFactors = F)
colnames(assoc_t1d)
head(assoc_t1d)
###Subset
assoc_t1d <- subset(assoc_t1d, select = c("Rs.Id", "Disease", "P.Value",  "Odds.Ratio"))
#Rename to merge
colnames(assoc_t1d) <- c("Top_SNP", "Disease", "P.Value", "OR")

####Merge
t1d_res <- merge(assoc_t1d, t1d_df_merg, by = "Top_SNP") 


##Changing order
t1d_res <- t1d_res[,c(5:7, 1, 8, 9, 2:4,10:16)]

##Ordering by p-value
t1d_res <- t1d_res[order(t1d_res$P.Value),]

#Writing table
write.table(t1d_res, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/results/t1d_results.txt", sep = '\t', quote = F, row.names = F)

##################################################
#############UC##################################

uc_ld <- read.delim("UC_proxies_haploreg.txt", header = T)
#re-order columns for better visualize
uc_ld <- uc_ld[,c(1:4,33,6,5,7:32)]
colnames(uc_ld)
##Col query_snp_rsid is the associated SNP
##Colr rsID is the LD SNP

##Intersect 
#With proxies
inter_prox_uc <- intersect(uc_ld$rsID, navyD$SNP)
length(unique(inter_prox_uc))

##Merge intersection with p-values and mirnas scores
#subset only mirSNPs intersected
uc_sub <- subset(uc_ld[which(uc_ld$rsID %in% inter_prox_uc),])
length(unique(uc_sub$rsID))
colnames(uc_sub)
#Check 
all(uc_sub$rsID %in% inter_prox_uc)

###process ld datafucme
uc_df <- subset(uc_sub, select = c("chr", "pos_hg38", "query_snp_rsid","rsID" , "r2",  "D.", "ref", "alt", "EUR", "GENCODE_name", "eQTL"))

colnames(uc_df) <- c("chr", "pos_hg38", "query_snp_rsid","SNP" , "r2",  "D.","ref", "alt", "MAF","GENCODE_name", "eQTL")
head(uc_df)

uc_df_merg <- merge(uc_df, navyD, by = "SNP")
colnames(uc_df_merg)

###Re-ordering columns e subseting
uc_df_merg <- subset(uc_df_merg, select = c("SNP", "chr", "pos_hg38", "query_snp_rsid", "r2","D.","ref", "alt",
                                             "MAF", "GENCODE_name", 
                                             "miRNA", "combined_score", "eQTL"))
##Ordering by P-values
uc_df_merg <- uc_df_merg[order(uc_df_merg$PValue),]

####Change column name 
colnames(uc_df_merg) <- c( "miRSNP","chr", "pos_hg38", "Top_SNP", "r2","D_prime","ref", "alt",
                           "MAF", "GENCODE_name", "miRNA", "combined_score",  "eQTL")


####Check which disease overlapps with miRSNPs
assoc_uc <- read.delim("/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/proceed_data/top_associated_snps_perdisease/UC_assoc_variants.txt", header = T, stringsAsFactors = F)
colnames(assoc_uc)
head(assoc_uc)
###Subset
assoc_uc <- subset(assoc_uc, select = c("Rs.Id", "Disease", "P.Value", "Odds.Ratio"))
#Rename to merge
colnames(assoc_uc) <- c("Top_SNP", "Disease", "P.Value", "OR")

####Merge
uc_res <- merge(assoc_uc, uc_df_merg, by = "Top_SNP") 

##Changing order
uc_res <- uc_res[,c(5:7, 1, 8, 9, 2:4,10:16)]

##Ordering by p-value
uc_res <- uc_res[order(uc_res$P.Value),]

#Writing table
write.table(uc_res, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/results/uc_results.txt", sep = '\t', quote = F, row.names = F)

##################################################
