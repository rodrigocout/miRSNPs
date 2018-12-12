#GEUVADIS miRNA mRNA correlations RNA-seq 312 individuos PBMCs
setwd ('/Users/Rodrigo/Documents/Pos_doc/Projeto_miRNA_SNPs/GEUVADIS_RNA_seq_public')
#Normalized miRNA data
mir <- read.delim('geuvadis_mature_sense.TMM.txt', head = T, row.names = 1, as.is=T)
mir <- as.matrix(mir)
str(mir)
dim(mir)

#Normalized mRNA data
rna <- read.delim('mRNA_expression_table.genelevel.v71.htseq.batch9.TMM.geuvadis.txt', head = T, row.names = 1, as.is=T) 
rna <- as.matrix(rna)
str(rna)
dim(rna)

id <- read.delim('miRNA_to_mRNA_sample_conversion.txt',head = F, as.is = T)

mir_df <- mir[,which(colnames(mir) %in% id[,1])]

mir_df <- mir_df[,order(colnames(mir_df))]

id <- id[order(id[,1]),]

all(id[,1] == colnames(mir_df))

colnames(mir_df) <- id[,2]

rna_df <- rna[,order(colnames(rna))]
mir_df <- mir_df[,order(colnames(mir_df))]

mir_df <- mir_df[,which(colnames(mir_df) %in% colnames(rna_df))]
rna_df <- rna_df[,which(colnames(rna_df) %in% colnames(mir_df))]
all(colnames(mir_df) == colnames(rna_df))

dim(mir_df)
dim(rna_df)

#UBE2L3
mir_4741 <- mir_df[c('hsa-miR-4741'),]
mir_4741 <- as.matrix(mir_4741)
length(mir_4741)

ube2l3 <- rna_df[c('ENSG00000185651'),]
ube2l3 <- as.matrix(ube2l3)
ube2l3 <- log2(ube2l3)
length(ube2l3)
head(ube2l3)

mir_1207 <- mir_df[c('hsa-miR-1207-5p'),]
mir_1207 <- as.matrix(mir_1207)
length(mir_1207)

mir_4763 <- mir_df[c('hsa-miR-4763-3p'),]
mir_4763 <- as.matrix(mir_4763)

mir_3918 <- mir_df[c('hsa-miR-3918'),]   
mir_3918 <- as.matrix(mir_3918)
mir_3918 <- log2(mir_3918)
head(mir_3918)

plot(ube2l3, mir_4741)
plot(mir_4741, ube2l3, xlim = c(5,20))
plot(mir_1207, ube2l3)
plot(mir_4763, ube2l3)
plot(mir_3918,ube2l3, main = "Norm Log2 expression")


cor(mir_3181, masp1, method = "spearman", use = "pairwise")
cor(mir_3181, masp1, method = "pearson", use = "pairwise")
cor.test(mir_3181, masp1, method = "pearson")
cor.test(mir_3181, masp1, method = "spearman", use = "pairwise")

pdf("./correlationplotMASP1.pdf")

plot(mir_3181, masp1, main = 'miRNA/mRNA expression', ylab = c('MASP1'), xlab = c('miR-3181'), 
     pch = 21, lwd=1, bty='n', col = 'dodgerblue')
dev.off()

#correlation matrix
cor <- read.delim('geuvadis_miRNA_mRNA_correlation.txt', head = T)
dim(cor)
class(cor)

mirnas_ai <- read.delim('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/shared_aid/mirnas.txt', header = F)
dim(mirnas_ai)
head(mirnas_ai)
class(mirnas_ai)
mirnas_ai <- as.matrix(mirnas_ai)
mirnas_ai <- as.vector(mirnas_ai[1:251,])


mirs_exp <- mir_df[mirnas_ai,]
head(mirs_exp)
head(mir_df)
