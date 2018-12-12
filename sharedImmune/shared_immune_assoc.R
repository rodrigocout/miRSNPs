###Shared association###########################
######Association Immunochip######################
cd_assoc <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/raw_data/CRO_assoc_variants.txt', header = T, stringsAsFactors = F)
head(cd_assoc)
colnames(cd_assoc)

cel_assoc <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/raw_data/CEL_assoc_variants.txt', header = T, stringsAsFactors = F)

atd_assoc <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/raw_data/ATD_assoc_variants.txt', header = T, stringsAsFactors = F)

jia_assoc <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/raw_data/JIA_assoc_variants.txt', header = T, stringsAsFactors = F)

ms_assoc <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/raw_data/MS_assoc_variants.txt', header = T, stringsAsFactors = F)

nar_assoc <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/raw_data/NAR_assoc_variants.txt', header = T, stringsAsFactors = F)

pbc_assoc <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/raw_data/PBC_assoc_variants.txt', header = T, stringsAsFactors = F)

pso_assoc <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/raw_data/PSO_assoc_variants.txt', header = T, stringsAsFactors = F)

ra_assoc <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/raw_data/RA_assoc_variants.txt', header = T, stringsAsFactors = F)

sle_assoc <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/raw_data/SLE_assoc_variants.txt', header = T, stringsAsFactors = F)

t1d_assoc <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/raw_data/T1D_assoc_variants.txt', header = T, stringsAsFactors = F)

uc_assoc <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/raw_data/UC_assoc_variants.txt', header = T, stringsAsFactors = F)

###########################################################################
##### Subset dataframes to use############################################
colnames(cd_assoc)
cro_assoc_df <- subset(cd_assoc, select = c("Region", "Rs.Id", "Position", "Alleles", "Disease",    
                                            "P.Value", "Odds.Ratio" ))
atd_assoc_df <- subset(atd_assoc, select = c("Region", "Rs.Id", "Position", "Alleles", "Disease",    
                                            "P.Value", "Odds.Ratio" ))

cel_assoc_df <- subset(cel_assoc, select = c("Region", "Rs.Id", "Position", "Alleles", "Disease",    
                                            "P.Value", "Odds.Ratio" ))

jia_assoc_df <- subset(jia_assoc, select = c("Region", "Rs.Id", "Position", "Alleles", "Disease",    
                                            "P.Value", "Odds.Ratio" ))

ms_assoc_df <- subset(ms_assoc, select = c("Region", "Rs.Id", "Position", "Alleles", "Disease",    
                                           "P.Value", "Odds.Ratio" ))

nar_assoc_df <- subset(nar_assoc, select = c("Region", "Rs.Id", "Position", "Alleles", "Disease",    
                                            "P.Value", "Odds.Ratio" ))

pbc_assoc_df <- subset(pbc_assoc, select = c("Region", "Rs.Id", "Position", "Alleles", "Disease",    
                                            "P.Value", "Odds.Ratio" ))

pso_assoc_df <- subset(pso_assoc, select = c("Region", "Rs.Id", "Position", "Alleles", "Disease",    
                                            "P.Value", "Odds.Ratio" ))

ra_assoc_df <- subset(ra_assoc, select = c("Region", "Rs.Id", "Position", "Alleles", "Disease",    
                                           "P.Value", "Odds.Ratio" ))

sle_assoc_df <- subset(sle_assoc, select = c("Region", "Rs.Id", "Position", "Alleles", "Disease",    
                                            "P.Value", "Odds.Ratio" ))

t1d_assoc_df <- subset(t1d_assoc, select = c("Region", "Rs.Id", "Position", "Alleles", "Disease",    
                                            "P.Value", "Odds.Ratio" ))

uc_assoc_df <- subset(uc_assoc, select = c("Region", "Rs.Id", "Position", "Alleles", "Disease",    
                                           "P.Value", "Odds.Ratio" ))

####################
###MirSNPs per disease##############
load(file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/mirsnps_ichip/atd_navyscore.RData")
load(file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/mirsnps_ichip/cd_mirsnps.RData")
load(file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/mirsnps_ichip/cel_mirsnps.RData")
load(file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/mirsnps_ichip/jia_mirsnps.RData")
load(file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/mirsnps_ichip/ms_mirsnps.RData")
load(file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/mirsnps_ichip/nar_mirsnps.RData")
load(file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/mirsnps_ichip/pbc_mirsnps.RData")
load(file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/mirsnps_ichip/pso_mirsnps.RData")
load(file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/mirsnps_ichip/ra_mirsnps.RData")
load(file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/mirsnps_ichip/sle_mirsnps.RData")
load(file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/mirsnps_ichip/t1d_mirsnps.RData")
load(file = "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/mirsnps_ichip/uc_mirsnps.RData")

#############################################################
####Intersect data ############################
###Replace colnames to merge################
colnames(cro_assoc_df) <- c("Region", "SNP", "Position",  "Alleles", "Disease", 
                            "P.Value" , "Odds.Ratio")

colnames(atd_assoc_df) <- c("Region", "SNP", "Position",  "Alleles", "Disease", 
                            "P.Value" , "Odds.Ratio")

colnames(cel_assoc_df) <- c("Region", "SNP", "Position",  "Alleles", "Disease", 
                            "P.Value" , "Odds.Ratio")

colnames(jia_assoc_df) <- c("Region", "SNP", "Position",  "Alleles", "Disease", 
                            "P.Value" , "Odds.Ratio")
colnames(ms_assoc_df) <- c("Region", "SNP", "Position",  "Alleles", "Disease", 
                           "P.Value" , "Odds.Ratio")

colnames(nar_assoc_df) <- c("Region", "SNP", "Position",  "Alleles", "Disease", 
                            "P.Value" , "Odds.Ratio")

colnames(pbc_assoc_df) <- c("Region", "SNP", "Position",  "Alleles", "Disease", 
                            "P.Value" , "Odds.Ratio")

colnames(pso_assoc_df) <- c("Region", "SNP", "Position",  "Alleles", "Disease", 
                            "P.Value" , "Odds.Ratio")

colnames(ra_assoc_df) <- c("Region", "SNP", "Position",  "Alleles", "Disease", 
                           "P.Value" , "Odds.Ratio")

colnames(sle_assoc_df) <- c("Region", "SNP", "Position",  "Alleles", "Disease", 
                            "P.Value" , "Odds.Ratio")

colnames(t1d_assoc_df) <- c("Region", "SNP", "Position",  "Alleles", "Disease", 
                            "P.Value" , "Odds.Ratio")

colnames(uc_assoc_df) <- c("Region", "SNP", "Position",  "Alleles", "Disease", 
                           "P.Value" , "Odds.Ratio")

####Concatnate with mirSNP data #################

cro_navy_assoc <- merge(cro_assoc_df, cd_navy, by = 'SNP', all = TRUE)
write.table(cro_navy_assoc, "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/mirsnps_assoc/CRO_mirsnps_assoc.txt", sep = '\t', quote = FALSE, row.names = FALSE)

atd_navy_assoc <- merge(atd_assoc_df, atd_navy, by = 'SNP', all = TRUE)
write.table(atd_navy_assoc, "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/mirsnps_assoc/ATD_mirsnps_assoc.txt", sep = '\t', quote = FALSE, row.names = FALSE)

cel_navy_assoc <- merge(cel_assoc_df, cel_navy, by = 'SNP', all = TRUE)
write.table(atd_navy_assoc, "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/mirsnps_assoc/CEL_mirsnps_assoc.txt", sep = '\t', quote = FALSE, row.names = FALSE)

jia_navy_assoc <- merge(jia_assoc_df, jia_navy, by = 'SNP', all = TRUE)
write.table(jia_navy_assoc, "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/mirsnps_assoc/JIA_mirsnps_assoc.txt", sep = '\t', quote = FALSE, row.names = FALSE)

ms_navy_assoc <- merge(ms_assoc_df, ms_navy, by = 'SNP', all = TRUE)
write.table(ms_navy_assoc, "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/mirsnps_assoc/MS_mirsnps_assoc.txt", sep = '\t', quote = FALSE, row.names = FALSE)

nar_navy_assoc <- merge(nar_assoc_df, nar_navy, by = 'SNP', all = TRUE)
write.table(nar_navy_assoc, "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/mirsnps_assoc/NAR_mirsnps_assoc.txt", sep = '\t', quote = FALSE, row.names = FALSE)

pbc_navy_assoc <- merge(pbc_assoc_df, pbc_navy, by = 'SNP', all = TRUE)
write.table(pbc_navy_assoc, "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/mirsnps_assoc/PBC_mirsnps_assoc.txt", sep = '\t', quote = FALSE, row.names = FALSE)

pso_navy_assoc <- merge(pso_assoc_df, pso_navy, by = 'SNP', all = TRUE)
write.table(pso_navy_assoc, "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/mirsnps_assoc/PSO_mirsnps_assoc.txt", sep = '\t', quote = FALSE, row.names = FALSE)

ra_navy_assoc <- merge(ra_assoc_df, ra_navy, by = 'SNP', all = TRUE)
write.table(ra_navy_assoc, "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/mirsnps_assoc/RA_mirsnps_assoc.txt", sep = '\t', quote = FALSE, row.names = FALSE)

sle_navy_assoc <- merge(sle_assoc_df, sle_navy, by = 'SNP', all = TRUE)
write.table(sle_navy_assoc, "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/mirsnps_assoc/SLE_mirsnps_assoc.txt", sep = '\t', quote = FALSE, row.names = FALSE)

t1d_navy_assoc <- merge(t1d_assoc_df, t1d_navy, by = 'SNP', all = TRUE)
write.table(t1d_navy_assoc, "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/mirsnps_assoc/T1D_mirsnps_assoc.txt", sep = '\t', quote = FALSE, row.names = FALSE)

uc_navy_assoc <- merge(uc_assoc_df, uc_navy, by = 'SNP', all = TRUE)
write.table(uc_navy_assoc, "/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/mirsnps_assoc/UC_mirsnps_assoc.txt", sep = '\t', quote = FALSE, row.names = FALSE)

###Calculate LD with top SNP#################
library(proxysnps)
atd_ld <- get_proxies(query = "rs2843403", pop = "CEU")

atd_
plot(atd$POS, atd$R.squared, main="rs2843403", xlab="Position", ylab=bquote("R"^2))

#or
#d <- get_proxies(chrom = "12", pos = 583090, window_size = 1e5, pop = "AFR")

