####eQTL data integration with miRSNPs per disease###############

blood_cis_eqtl <- read.table("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/eqtl_data/blood_eqtl/Blood_eQTL_CisAssociationsProbeLevelFDR0.5.txt", header = T)

blood_trans_eqtld <- read.table("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/eqtl_data/blood_eqtl/Blood_eQTL_TransEQTLsFDR0.5.txt", header = T)

tcd4_cis_eqtl_eu <- read.table("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/eqtl_data/eqtls_tcd4/tableS4_eu_cd4T_cis_fdr05.tsv", header =T, stringsAsFactors = F)

tcd4_cis_eqtl_ea <- read.table("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/eqtl_data/eqtls_tcd4/tableS4_ea_cd4T_cis_fdr05.tsv", header =T, stringsAsFactors = F)

tcd4_cis_eqtl_aa <- read.table("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/eqtl_data/eqtls_tcd4/tableS4_aa_cd4T_cis_fdr05.tsv", header =T, stringsAsFactors = F)

gtex_ebv <- read.delim("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/eqtl_data/gtex/GTEx_Analysis_V6_eQTLs/Cells_EBV-transformed_lymphocytes_Analysis.snpgenes", header =T, stringsAsFactors = F)
head(gtex_ebv)
class(gtex_ebv)
colnames(gtex_ebv)

gtex_blood <- read.delim("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/eqtl_data/gtex/GTEx_Analysis_V6_eQTLs/Whole_Blood_Analysis.snpgenes", header =T, stringsAsFactors = F)
colnames(gtex_blood)
dim(gtex_blood)
head(gtex_blood)

gtex_intestine <- read.delim("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/eqtl_data/gtex/GTEx_Analysis_V6_eQTLs/Small_Intestine_Terminal_Ileum_Analysis.snpgenes", header =T, stringsAsFactors = F)
colnames(gtex_intestine)

##################################################################################
###MirSNPs from Ichip data in each disease#####################

