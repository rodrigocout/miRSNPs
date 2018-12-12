#GEUVADIS miRNA mRNA correlations
setwd ('/home/lgmh/rodrigo/resources/RNA_seq_public_data/GEUVADIS_RNA_seq_public')

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

##Make all IDs equal to the mRNA data#####
id <- read.delim('miRNA_to_mRNA_sample_conversion.txt',head = F, as.is = T)

mir_df <- mir[,which(colnames(mir) %in% id[,1])]

mir_df <- mir_df[,order(colnames(mir_df))]

id <- id[order(id[,1]),]

all(id[,1] == colnames(mir_df))

colnames(mir_df) <- id[,2]

###Now order mRNA and microRNA data#####
rna_df <- rna[,order(colnames(rna))]
mir_df <- mir_df[,order(colnames(mir_df))]

mir_df <- mir_df[,which(colnames(mir_df) %in% colnames(rna_df))]
rna_df <- rna_df[,which(colnames(rna_df) %in% colnames(mir_df))]
all(colnames(mir_df) == colnames(rna_df))

###################
##Write tables####
mir_df <- as.data.frame(mir_df)
write.table(mir_df, "geuvadis_miRNAS_matched_rna.txt", sep = '\t', quote = F, row.names = T)

rna_df <- as.data.frame(rna_df)
write.table(rna_df, "geuvadis_RNAS_matched_mirna.txt", sep = '\t', quote = F, row.names = T)

###############################################################
cor(t(mir_df), t(rna_df))

cor <- read.delim('geuvadis_miRNA_mRNA_correlation.txt', head = T)

cor_ube2l3 <- cor[cor$Matrix2RowName == 'ENSG00000185651' & cor$P.value < 0.05,]
class(cor_ube2l3)

write.table(cor_ube2l3, "Correlation_table_UBE2L3.txt", sep = '\t', quote = )
#################################################
##Read data with gene expression and microRNAS expression##
################################################
gxrna <- read.delim("~/Dropbox/shared_immune_mirsnps/data/Gene_names_mirsnps_geuvadis.txt")
#gxrna <- gxrna[rowMeans(gxrna) > 1,]

gxmir <- read.delim("~/Dropbox/shared_immune_mirsnps/data/geuvadis_miRNAS_matched_rna.txt")
#gxmir <- gxmir[rowMeans(gxmir) > 1,]

###Now take the 1000 Genomes IDs######
gv <- read.delim("~/Dropbox/shared_immune_mirsnps/data/E-GEUV-1.sdrf.txt", header = T)
class(gv)
colnames(gv)
gv <- gv[,c(1,28,33)]
colnames(gv) <- c('1KG', 'ID', 'Pop')
gv <- unique(gv)

mir_1kg <- gxmir[,which(colnames(gxmir) %in% gv[,2])]

mir_1kg <- mir_1kg[,order(colnames(mir_1kg))]

gv <- gv[which(gv[,2] %in% colnames(gxmir)),]

gv <- gv[order(gv[,2]),]

all(gv[,2] == colnames(mir_1kg))
all(gv$ID %in% colnames(mir_1kg))

colnames(mir_1kg) <- gv[,1]


rna_1kg <- gxrna[,which(colnames(gxrna) %in% gv[,2])]

rna_1kg <- rna_1kg[,order(colnames(rna_1kg))]

all(gv[,2] == colnames(rna_1kg))
all(gv$ID %in% colnames(rna_1kg))

colnames(rna_1kg) <- gv[,1]

write.table(mir_1kg, "~/Dropbox/shared_immune_mirsnps/data/microRNAs_Geuvadis_1KG_expression.txt", sep = '\t', quote = F, row.names = T)

write.table(rna_1kg, "~/Dropbox/shared_immune_mirsnps/data/mRNAs_Geuvadis_1KG_expression.txt", sep = '\t', quote = F, row.names = T)

#geuv_1kg <- read.delim("list_ind_mirsnps_geuvadis_1kg.txt", header = F)

#mirs_1000G <- read.delim("microRNAs_Geuvadis_1KG_expression.txt", header = T)

#rnas_1000G <- read.delim("mRNAs_Geuvadis_1KG_expression.txt", header = T)

###Read genotypes ############################
map <- read.table('~/Dropbox/shared_immune_mirsnps/data/mirnas_geuvadis_compound.map', stringsAsFactors = F)
map <- map[,c(1,2,4)]
colnames(map) <- c('chr', 'SNP', 'BP')
head(map)

snps <- read.table('~/Dropbox/shared_immune_mirsnps/data/mirnas_geuvadis_compound.ped', stringsAsFactors = F, row.names = 1)  
head(snps, n = 2)
dim(snps)
##1 = A##2 = C##3 = G##4 = T
snps <- snps[,c(6:34)]
rownames(snps)
colnames(snps)

head(map)
length(map$SNP)

length(colnames(snps))

colnames(snps) <- map$SNP

all(map$SNP == colnames(snps) )

##Subset miRNAs samples
mirsnps <- mir_1kg[,which(colnames(mir_1kg) %in% rownames(snps))]

dim(mirsnps)
colnames(mirsnps)

##Subset mRNA samples
rna_snps <- rna_1kg[, which(colnames(rna_1kg) %in% rownames(snps))]


all(rownames(snps) %in% colnames(mirsnps))
all(rownames(snps) %in% colnames(rna_snps))

write.table(mirsnps, "~/Dropbox/shared_immune_mirsnps/data/micronas_snps_geuvadis.txt", sep = '\t', quote = F)

write.table(rna_snps, "~/Dropbox/shared_immune_mirsnps/data/mRNAs_snps_geuvadis.txt", sep = '\t', quote = F)

##################################

mirsnps <- read.delim("~/Dropbox/shared_immune_mirsnps/data/micronas_snps_geuvadis.txt")
rna_snps <- read.delim("~/Dropbox/shared_immune_mirsnps/data/mRNAs_snps_geuvadis.txt")

all(colnames(rna_snps) %in% colnames(mirsnps))
all(colnames(rna_snps) == colnames(mirsnps))

lista <- read.table("~/Dropbox/shared_immune_mirsnps/tables/results/lista_microRNAs.txt")

miRInter <- intersect(rownames(mirsnps), lista$V1)
class(miRInter)

##----gexp intersections
ids <- colnames (mirsnps)
#rna_snps <- rna_snps[,ids]
#mirsnps_mrna <- rbind (rna_snps[gensInter,], mirsnps[miRInter, ])

mirsnps <- mirsnps[miRInter,]
#mirsnps <- log2(mirsnps)
head(mirsnps)
rownames(mirsnps)

#df <-  mirsnps[rowMeans(mirsnps) > 1,]

#make all the same number of samples
snps <- snps[which(rownames(snps) %in% colnames(mirsnps)),]
head(snps)

##########################################################
snps$rs1054037[snps$rs1054037 == 'CC'] <- '0'
snps$rs1054037[snps$rs1054037 == 'CT'] <- '1'
snps$rs1054037[snps$rs1054037 == 'TT'] <- '2'

e1 = as.numeric(mirsnps[c("hsa-miR-660-5p"),])
s1 = as.numeric(snps$rs1054037)

lm1 = lm(e1 ~ s1)

##Plot the data
pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_mir660_5p.pdf")

boxplot(e1 ~ s1, lwd = 2, ylab = 'Normalized Expression', xlab = "rs1054037", main = "hsa-miR-660-5p",
        axes = F)
stripchart(e1 ~ s1, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("CC","CT","TT"))
axis(side=2)
lines(lm1$fitted ~ s1,type="b",pch=15,col="darkgrey")

dev.off()
####################
e2 = as.numeric(mirsnps["hsa-miR-7151-5p",])

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_mir7151_5p.pdf")
boxplot(e2 ~ s1, lwd = 2, ylab = 'Normalized Expression', xlab = "rs1054037", main = "hsa-miR-7151-5p",
        axes = F)
stripchart(e2 ~ s1, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("CC","CT","TT"))
axis(side=2)

dev.off()
#####################################################
e3 = as.numeric(mirsnps["hsa-miR-3686",])

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_mir3686.pdf")
boxplot(e3 ~ s1, lwd = 2, ylab = 'Normalized Expression', xlab = "rs1054037", main = "hsa-miR-3686",
        axes = F)
stripchart(e3 ~ s1, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("CC","CT","TT"))
axis(side=2)
dev.off()
########################################################
##miRSNP rs4013 ################################
snps$rs4013[snps$rs4013 == 'CC'] <- '0'
snps$rs4013[snps$rs4013 == 'CT'] <- '1'
snps$rs4013[snps$rs4013 == 'TT'] <- '2'

e4 = as.numeric(mirsnps[c("hsa-miR-4742-3p"),])
s2 = as.numeric(snps$rs4013)

##Plot the data
pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_mir4742_3p.pdf")
boxplot(e4 ~ s2, lwd = 2, ylab = 'Normalized Expression', xlab = "rs4013", main = "hsa-miR-4742-3p",
        axes = F)
stripchart(e4 ~ s2, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("CC","CT","TT"))
axis(side=2)
dev.off()
#------------
e5 = as.numeric(mirsnps[c("hsa-miR-630"),])

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_mir630.pdf")
boxplot(e5 ~ s2, lwd = 2, ylab = 'Normalized Expression', xlab = "rs4013", main = "hsa-miR-630",
        axes = F)
stripchart(e5 ~ s2, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("CC","CT","TT"))
axis(side=2)
dev.off()

#######################################
##miRSNP rs1054029
snps$rs1054029[snps$rs1054029 == 'GG'] <- '0'
snps$rs1054029[snps$rs1054029 == 'GA'] <- '1'
snps$rs1054029[snps$rs1054029 == 'AA'] <- '2'

e6 = as.numeric(mirsnps[c("hsa-miR-124-5p"),])
s3 = as.numeric(snps$rs1054029)

##Plot the data
pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_mir124_5p.pdf")
boxplot(e6 ~ s3, lwd = 2, ylab = 'Normalized Expression', xlab = "rs1054029", main = "hsa-miR-124-5p",
        axes = F)
stripchart(e6 ~ s3, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("GG","GA","AA"))
axis(side=2)
dev.off()
#----------------------
e7 = as.numeric(mirsnps[c("hsa-miR-4766-5p"),])

##Plot the data
pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_mir4766_5p.pdf")
boxplot(e7 ~ s3, lwd = 2, ylab = 'Normalized Expression', xlab = "rs1054029", main = "hsa-miR-4766-5p",
        axes = F)
stripchart(e7 ~ s3, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("GG","GA","AA"))
axis(side=2)
dev.off()
#---------------------
e8 = as.numeric(mirsnps[c("hsa-miR-124"),])

##Plot the data
pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_mir124.pdf")
boxplot(e8 ~ s3, lwd = 2, ylab = 'Normalized Expression', xlab = "rs1054029", main = "hsa-miR-124",
        axes = F)
stripchart(e8 ~ s3, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("GG","GA","AA"))
axis(side=2)
dev.off()
#########################################
#miRSNP rs39602

snps$rs39602[snps$rs39602 == 'CC'] <- '0'
snps$rs39602[snps$rs39602 == 'CG'] <- '1'
snps$rs39602[snps$rs39602 == 'GG'] <- '2'

e9 = as.numeric(mirsnps[c("hsa-miR-6800-5p"),])
s4 = as.numeric(snps$rs39602)

##Plot the data
pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_mir6800_5p.pdf")
boxplot(e9 ~ s4, lwd = 2, ylab = 'Normalized Expression', xlab = "rs39602", main = "hsa-miR-6800-5p",
        axes = F)
stripchart(e9 ~ s4, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("CC","CG","GG"))
axis(side=2)
dev.off()
######################################################
##miRSNPs rs2070197

snps$rs2070197[snps$rs2070197 == 'TT'] <- '0'
snps$rs2070197[snps$rs2070197 == 'CT'] <- '1'
snps$rs2070197[snps$rs2070197 == 'CC'] <- '2'

e10 = as.numeric(mirsnps[c("hsa-miR-3136-3p"),])
s5 = as.numeric(snps$rs2070197)

##Plot the data
pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs2070197_mir3136.pdf")
boxplot(e10 ~ s5, lwd = 2, ylab = 'Normalized Expression', xlab = "rs2070197", main = "hsa-miR-3136-3p",
        axes = F)
stripchart(e10 ~ s5, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("TT","CT","CC"))
axis(side=2)
dev.off()
#------------------------------
e11 = as.numeric(mirsnps[c("hsa-miR-7155-3p"),])

##Plot the data
pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs2070197_mir7155.pdf")
boxplot(e11 ~ s5, lwd = 2, ylab = 'Normalized Expression', xlab = "rs2070197", main = "hsa-miR-7155-3p",
        axes = F)
stripchart(e10 ~ s5, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("TT","CT","CC"))
axis(side=2)
dev.off()
##############################################################
##MiRSNP rs10114470 #########################################

snps$rs10114470[snps$rs10114470 == 'TT'] <- '0'
snps$rs10114470[snps$rs10114470 == 'TC'] <- '1'
snps$rs10114470[snps$rs10114470 == 'CC'] <- '2'

e11 = as.numeric(mirsnps[c("hsa-miR-376a-3p"),])
s6 = as.numeric(snps$rs10114470)

##Plot the data
pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs10114470_hsa-miR-376a-3p.pdf")

boxplot(e11 ~ s6, lwd = 2, ylab = 'Normalized Expression', xlab = "rs10114470", main = "hsa-miR-376a-3p",
        axes = F)
stripchart(e11 ~ s6, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("TT","TC","CC"))
axis(side=2)

dev.off()
#-----------------
e12 = as.numeric(mirsnps["hsa-miR-4753-5p",])

##Plot the data
pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs10114470_hsa-miR-4753-5p.pdf")

boxplot(e12 ~ s6, lwd = 2, ylab = 'Normalized Expression', xlab = "rs10114470", main = "hsa-miR-4753-5p",
        axes = F)
stripchart(e12 ~ s6, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("TT","TC","CC"))
axis(side=2)

dev.off()
###################################################
####miRSNP rs3088081 #############################

snps$rs3088081[snps$rs3088081 == 'AA'] <- '0'
snps$rs3088081[snps$rs3088081 == 'GA'] <- '1'
snps$rs3088081[snps$rs3088081 == 'GG'] <- '2'

e13 = as.numeric(mirsnps[c("hsa-miR-3661"),])
s7 = as.numeric(snps$rs3088081)

##Plot the data
pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs3088081_hsa-miR-3661.pdf")

boxplot(e13 ~ s7, lwd = 2, ylab = 'Normalized Expression', xlab = "rs3088081", main = "hsa-miR-3661",
        axes = F)
stripchart(e13 ~ s7, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("AA","GA","GG"))
axis(side=2)
dev.off()
#######################################
##miRSNP rs9943 ######################

snps$rs9943[snps$rs9943 == 'AA'] <- '0'
snps$rs9943[snps$rs9943 == 'GA'] <- '1'
snps$rs9943[snps$rs9943 == 'GG'] <- '2'

e14 = as.numeric(mirsnps[c("hsa-miR-628-5p"),])
s8 = as.numeric(snps$rs9943)
lm8 = lm(e14 ~ s8)

summary(lm8)
names(summary(lm8))

##Plot the data
pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs9943_hsa-miR-628-5p.pdf")

boxplot(e14 ~ s8, lwd = 2, ylab = 'Normalized Expression', xlab = "rs9943", main = "hsa-miR-628-5p",
        axes = F)
stripchart(e14 ~ s8, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("AA","GA","GG"))
axis(side=2)
lines(lm8$fitted ~ s8,type="b", pch=15,col="darkgrey")

dev.off()
######################################################
###miRSNP rs907091###################################
snps$rs907091[snps$rs907091 == 'CC'] <- '0'
snps$rs907091[snps$rs907091 == 'TC'] <- '1'
snps$rs907091[snps$rs907091 == 'TT'] <- '2'

e15 = as.numeric(mirsnps[c("hsa-miR-4497"),])
s9 = as.numeric(snps$rs907091)

##Plot the data
pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs907091_hsa-miR-4497.pdf")

boxplot(e15 ~ s9, lwd = 2, ylab = 'Normalized Expression', xlab = "rs907091", main = "hsa-miR-4497",
        axes = F, ylim = c(0, 12))
stripchart(e15 ~ s9, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("CC","TC","TT"), )
axis(side=2)
dev.off()
#--------------------
e16 = as.numeric(mirsnps[c("hsa-miR-4518"),])

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs907091_hsa-miR-4518.pdf")

boxplot(e16 ~ s9, lwd = 2, ylab = 'Normalized Expression', xlab = "rs907091", main = "hsa-miR-4518",
        axes = F)
stripchart(e16 ~ s9, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("CC","TC","TT"), )
axis(side=2)
dev.off()
#----------------------
e17 = as.numeric(mirsnps["hsa-miR-330-5p",])

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs907091_hsa-miR-330-5p.pdf")

boxplot(e17 ~ s9, lwd = 2, ylab = 'Normalized Expression', xlab = "rs907091", main = "hsa-miR-330-5p",
        axes = F)
stripchart(e17 ~ s9, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("CC","TC","TT"), )
axis(side=2)
dev.off()
#---------------------
e18 = as.numeric(mirsnps["hsa-miR-326",])

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs907091_hsa-miR-326.pdf")

boxplot(e18 ~ s9, lwd = 2, ylab = 'Normalized Expression', xlab = "rs907091", main = "hsa-miR-326",
        axes = F)
stripchart(e18 ~ s9, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("CC","TC","TT"), )
axis(side=2)
dev.off()
#########################################
#miRSNP rs16940681 #####################
snps$rs16940681[snps$rs16940681 == 'GG'] <- '0'
snps$rs16940681[snps$rs16940681 == 'CG'] <- '1'
snps$rs16940681[snps$rs16940681 == 'CC'] <- '2'

e19 = as.numeric(mirsnps[c("hsa-miR-6740-5p"),])
s10 = as.numeric(snps$rs16940681)

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs16940681_hsa-miR-6740-5p.pdf")

boxplot(e19 ~ s10, lwd = 2, ylab = 'Normalized Expression', xlab = "rs16940681", main = "hsa-miR-6740-5p",
        axes = F)
stripchart(e19 ~ s10, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("GG","CG","CC"), )
axis(side=2)
dev.off()

##########################################
#miRSNP rs2316765 ########################

snps$rs2316765[snps$rs2316765 == 'TT'] <- '0'
snps$rs2316765[snps$rs2316765 == 'CT'] <- '1'
snps$rs2316765[snps$rs2316765 == 'CC'] <- '2'

e20 = as.numeric(mirsnps[c("hsa-miR-30c-2-3p"),])
s11 = as.numeric(snps$rs2316765)

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs2316765_hsa-miR-30c-2-3p.pdf")

boxplot(e20 ~ s11, lwd = 2, ylab = 'Normalized Expression', xlab = "rs2316765", main = "hsa-miR-30c-2-3p",
        axes = F)
stripchart(e20 ~ s11, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("TT","CT","CC"), )
axis(side=2)
dev.off()
#-----------------------------
e21 = as.numeric(mirsnps["hsa-miR-30c-1-3p",])

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs2316765_hsa-miR-30c-1-3p.pdf")

boxplot(e21 ~ s11, lwd = 2, ylab = 'Normalized Expression', xlab = "rs2316765", main = "hsa-miR-30c-1-3p",
        axes = F)
stripchart(e21 ~ s11, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("TT","CT","CC"), )
axis(side=2)
dev.off()
#####################################
###miRSNP rs878886 #################

snps$rs878886[snps$rs878886 == 'CC'] <- '0'
snps$rs878886[snps$rs878886 == 'GC'] <- '1'
snps$rs878886[snps$rs878886 == 'GG'] <- '2'

e22 = as.numeric(mirsnps[c("hsa-miR-4685-5p"),])
s12 = as.numeric(snps$rs878886)


pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs878886_hsa-miR-4685-5p.pdf")
boxplot(e22 ~ s12, lwd = 2, ylab = 'Normalized Expression', xlab = "rs878886", main = "hsa-miR-4685-5p",
        axes = F)
stripchart(e22 ~ s12, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("CC","GC","GG"), )
axis(side=2)
dev.off()
#--------------------
e23 = as.numeric(mirsnps["hsa-miR-1915-3p",])

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs878886_hsa-miR-1915.pdf")
boxplot(e23 ~ s12, lwd = 2, ylab = 'Normalized Expression', xlab = "rs878886", main = "hsa-miR-1915",
        axes = F)
stripchart(e23 ~ s12, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("CC","GC","GG"), )
axis(side=2)
dev.off()
#--------------------
e24 = as.numeric(mirsnps["hsa-miR-3918",])

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs878886_hsa-miR-3918.pdf")
boxplot(e24 ~ s12, lwd = 2, ylab = 'Normalized Expression', xlab = "rs878886", main = "hsa-miR-3918",
        axes = F)
stripchart(e24 ~ s12, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("CC","GC","GG"), )
axis(side=2)
dev.off()
###################################
##miRSNP rs878887 ################

snps$rs878887[snps$rs878887 == 'CC'] <- '0'
snps$rs878887[snps$rs878887 == 'TC'] <- '1'
snps$rs878887[snps$rs878887 == 'TT'] <- '2'

e25 = as.numeric(mirsnps[c("hsa-miR-3186-5p"),])
s13 = as.numeric(snps$rs878887)

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs878887_hsa-miR-3186-5p.pdf")
boxplot(e25 ~ s13, lwd = 2, ylab = 'Normalized Expression', xlab = "rs878887", main = "hsa-miR-3186-5p",
        axes = F)
stripchart(e25 ~ s13, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("CC","TC","TT"), )
axis(side=2)
dev.off()
#-------------------------------------
e26 = as.numeric(mirsnps["hsa-miR-136",])

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs878887_hsa-miR-136.pdf")
boxplot(e26 ~ s13, lwd = 2, ylab = 'Normalized Expression', xlab = "rs878887", main = "hsa-miR-136",
        axes = F)
stripchart(e26 ~ s13, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("CC","TC","TT"), )
axis(side=2)
dev.off()
##############################################
#miRSNP rs878888 #############################

snps$rs878888[snps$rs878888 == 'AA'] <- '0'
snps$rs878888[snps$rs878888 == 'GA'] <- '1'
snps$rs878888[snps$rs878888 == 'GG'] <- '2'

e27 = as.numeric(mirsnps["hsa-miR-5708",])
s14 = as.numeric(snps$rs878888)

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs878888_hsa-miR-5708.pdf")
boxplot(e27 ~ s14, lwd = 2, ylab = 'Normalized Expression', xlab = "rs878888", main = "hsa-miR-5708",
        axes = F)
stripchart(e27 ~ s14, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("AA","GA","GG"), )
axis(side=2)
dev.off()
#-------------------------------------------
e28 = as.numeric(mirsnps["hsa-miR-1226-5p",])

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs878888_hsa-miR-1226-5p.pdf")
boxplot(e28 ~ s14, lwd = 2, ylab = 'Normalized Expression', xlab = "rs878888", main = "hsa-miR-1226-5p",
        axes = F)
stripchart(e28 ~ s14, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("AA","GA","GG"), )
axis(side=2)
dev.off()
###############################################
#miRSNP rs4640231 #############################

snps$rs4640231[snps$rs4640231 == 'GG'] <- '0'
snps$rs4640231[snps$rs4640231 == 'CG'] <- '1'
snps$rs4640231[snps$rs4640231 == 'CC'] <- '2'

e29 = as.numeric(mirsnps["hsa-miR-6755-5p",])
s15 = as.numeric(snps$rs4640231)

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs4640231_hsa-miR-6755-5p.pdf")
boxplot(e29 ~ s15, lwd = 2, ylab = 'Normalized Expression', xlab = "rs4640231", main = "hsa-miR-6755-5p",
        axes = F)
stripchart(e29 ~ s15, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("GG","CG","CC"))
axis(side=2)
dev.off()
##############################################
#mirsnp rs4482334 ############################

snps$rs4482334[snps$rs4482334 == 'TT'] <- '0'
snps$rs4482334[snps$rs4482334 == 'CT'] <- '1'
snps$rs4482334[snps$rs4482334 == 'CC'] <- '2'

e30 = as.numeric(mirsnps["hsa-miR-6890-5p",])
s16 = as.numeric(snps$rs4482334)

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs4482334_hsa-miR-6890-5p.pdf")
boxplot(e30 ~ s16, lwd = 2, ylab = 'Normalized Expression', xlab = "rs4482334", main = "hsa-miR-6890-5p",
        axes = F)
stripchart(e30 ~ s16, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("TT","CT","CC"))
axis(side=2)
dev.off()
#-----------------------------
e31 = as.numeric(mirsnps["hsa-miR-6742-5p",])

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs4482334_hsa-miR-6742-5p.pdf")
boxplot(e31 ~ s16, lwd = 2, ylab = 'Normalized Expression', xlab = "rs4482334", main = "hsa-miR-6742-5p",
        axes = F)
stripchart(e31 ~ s16, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("TT","CT","CC"))
axis(side=2)
dev.off()
#----------------
e32 = as.numeric(mirsnps["hsa-miR-4722-5p",])

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs4482334_hsa-miR-4722-5p.pdf")
boxplot(e32 ~ s16, lwd = 2, ylab = 'Normalized Expression', xlab = "rs4482334", main = "hsa-miR-4722-5p",
        axes = F)
stripchart(e32 ~ s16, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("TT","CT","CC"))
axis(side=2)
dev.off()
#---------------------------------
e33 = as.numeric(mirsnps["hsa-miR-6796-5p",])

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs4482334_miR-6796-5p.pdf")
boxplot(e33 ~ s16, lwd = 2, ylab = 'Normalized Expression', xlab = "rs4482334", main = "miR-6796-5p",
        axes = F)
stripchart(e33 ~ s16, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("TT","CT","CC"))
axis(side=2)
dev.off()
#-----------------------------------
e34 = as.numeric(mirsnps["hsa-miR-4459",])

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs4482334_miR-4459.pdf")
boxplot(e34 ~ s16, lwd = 2, ylab = 'Normalized Expression', xlab = "rs4482334", main = "miR-4459",
        axes = F)
stripchart(e34 ~ s16, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("TT","CT","CC"))
axis(side=2)
dev.off()
##################################################
#miRSNP rs12373168 ##############################

snps$rs12373168[snps$rs12373168 == 'AA'] <- '0'
snps$rs12373168[snps$rs12373168 == 'CA'] <- '1'
snps$rs12373168[snps$rs12373168 == 'CC'] <- '2'

e35 = as.numeric(mirsnps["hsa-miR-33b-3p",])
s17 = as.numeric(snps$rs12373168)

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs12373168_miR-33b-3p.pdf")
boxplot(e35 ~ s17, lwd = 2, ylab = 'Normalized Expression', xlab = "rs12373168", main = "miR-33b-3p",
        axes = F)
stripchart(e35 ~ s17, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("AA","CA","CC"))
axis(side=2)
dev.off()
#####################################################
#miRSNP rs45450798 #################################

snps$rs45450798[snps$rs45450798 == 'CC'] <- '0'
snps$rs45450798[snps$rs45450798 == 'GC'] <- '1'
snps$rs45450798[snps$rs45450798 == 'GG'] <- '2'

e36 = as.numeric(mirsnps["hsa-miR-4531",])
s18 = as.numeric(snps$rs45450798)

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs45450798_miR-4531.pdf")
boxplot(e36 ~ s18, lwd = 2, ylab = 'Normalized Expression', xlab = "rs45450798", main = "miR-4531",
        axes = F)
stripchart(e36 ~ s18, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("CC","CC","GG"))
axis(side=2)
dev.off()
############################################
##miRSNP rs9950174 ####################

snps$rs9950174[snps$rs9950174 == "TT"] <- "0"
snps$rs9950174[snps$rs9950174 == "CT"] <- "1"
snps$rs9950174[snps$rs9950174 == "CC"] <- "2"

e37 = as.numeric(mirsnps["hsa-miR-5189-3p",])
s19 = as.numeric(snps$rs9950174)

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs9950174_hsa-miR-5189-3p.pdf")
boxplot(e37 ~ s19, lwd = 2, ylab = 'Normalized Expression', xlab = "rs9950174", main = "miR-5189-3p",
        axes = F)
stripchart(e37 ~ s19, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("TT","CT","CC"))
axis(side=2)
dev.off()
#####################################
##mirSNP rs727088

snps$rs727088[snps$rs727088 == "GG"] <- "0"
snps$rs727088[snps$rs727088 == "AG"] <- "1"
snps$rs727088[snps$rs727088 == "AA"] <- "2"

e38 = as.numeric(mirsnps["hsa-miR-513a-3p",])
s20 = as.numeric(snps$rs727088)

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs727088_hsa-miR-513a-3p.pdf")
boxplot(e38 ~ s20, lwd = 2, ylab = 'Normalized Expression', xlab = "rs727088", main = "miR-513a-3p",
        axes = F)
stripchart(e38 ~ s20, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("TT","CT","CC"))
axis(side=2)
dev.off()
#############################
###miRSNP rs571689 ##########

snps$rs571689[snps$rs571689 == "CC"] <- "0"
snps$rs571689[snps$rs571689 == "TC"] <- "1"
snps$rs571689[snps$rs571689 == "TT"] <- "2"

e39 = as.numeric(mirsnps["hsa-miR-552-3p",])
s21 = as.numeric(snps$rs571689)

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs571689_hsa-miR-552-3p.pdf")
boxplot(e39 ~ s21, lwd = 2, ylab = 'Normalized Expression', xlab = "rs571689", main = "miR-552-3p",
        axes = F)
stripchart(e39 ~ s21, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("CC","TC","TT"))
axis(side=2)
dev.off()
##############################################
#miRSNP rs570794 ###############

snps$rs570794[snps$rs570794 == "TT"] <- "0"
snps$rs570794[snps$rs570794 == "CT"] <- "1"
snps$rs570794[snps$rs570794 == "CC"] <- "2"

e40 = as.numeric(mirsnps["hsa-miR-4430",])
s22 = as.numeric(snps$rs570794)

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs570794_hsa-miR-4430.pdf")
boxplot(e40 ~ s22, lwd = 2, ylab = 'Normalized Expression', xlab = "rs570794", main = "miR-4430",
        axes = F)
stripchart(e40 ~ s22, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("TT","CT","CC"))
axis(side=2)
dev.off()
#-----------
e41 = as.numeric(mirsnps["hsa-miR-1295b-5p",])

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs570794_hsa-miR-1295b_5p.pdf")
boxplot(e41 ~ s22, lwd = 2, ylab = 'Normalized Expression', xlab = "rs570794", main = "miR-1295b-5p",
        axes = F)
stripchart(e41 ~ s22, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("TT","CT","CC"))
axis(side=2)
dev.off()
#------------------------------------------
e42 = as.numeric(mirsnps["hsa-miR-4463",])

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs570794_hsa-miR-4463.pdf")
boxplot(e42 ~ s22, lwd = 2, ylab = 'Normalized Expression', xlab = "rs570794", main = "miR-4463",
        axes = F)
stripchart(e42 ~ s22, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("TT","CT","CC"))
axis(side=2)
dev.off()
#######################################
########miRSNP rs507766 ##############

snps$rs507766[snps$rs507766 == "TT"] <- "0"
snps$rs507766[snps$rs507766 == "CT"] <- "1"
snps$rs507766[snps$rs507766 == "CC"] <- "2"

e43 = as.numeric(mirsnps["hsa-miR-136-5p",])
s23 = as.numeric(snps$rs507766)

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs507766_hsa-miR-136-5p.pdf")
boxplot(e43 ~ s23, lwd = 2, ylab = 'Normalized Expression', xlab = "rs507766", main = "miR-136-5p",
        axes = F)
stripchart(e43 ~ s23, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("TT","CT","CC"))
axis(side=2)
dev.off()
#############################
#miRSNP rs506897 ###########

snps$rs506897[snps$rs506897 == "GG"] <- "0"
snps$rs506897[snps$rs506897 == "CG"] <- "1"
snps$rs506897[snps$rs506897 == "CC"] <- "2"

e44 = as.numeric(mirsnps["hsa-miR-4530",])
s24 = as.numeric(snps$rs506897)

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs506897_hsa-miR-4530.pdf")
boxplot(e44 ~ s24, lwd = 2, ylab = 'Normalized Expression', xlab = "rs506897", main = "miR-4530",
        axes = F)
stripchart(e44 ~ s24, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("GG","CG","CC"))
axis(side=2)
dev.off()
###########################################
##miRSNP rs503279 #######################

snps$rs503279[snps$rs503279 == "TT"] <- "0"
snps$rs503279[snps$rs503279 == "CT"] <- "1"
snps$rs503279[snps$rs503279 == "CC"] <- "2"

e45 = as.numeric(mirsnps["hsa-miR-675-3p",])
s25 = as.numeric(snps$rs503279)

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs503279_hsa-miR-675-3p.pdf")
boxplot(e45 ~ s25, lwd = 2, ylab = 'Normalized Expression', xlab = "rs503279", main = "miR-675-3p",
        axes = F)
stripchart(e45 ~ s25, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("TT","CT","CC"))
axis(side=2)
dev.off()
#############################
#miRSNP rs1056441 ######

snps$rs1056441[snps$rs1056441 == "TT"] <- "0"
snps$rs1056441[snps$rs1056441 == "TC"] <- "1"
snps$rs1056441[snps$rs1056441 == "CC"] <- "2"


e46 = as.numeric(mirsnps["hsa-miR-4745-3p",])
s26 = as.numeric(snps$rs1056441)

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs1056441_hsa-miR-4745-3p.pdf")
boxplot(e46 ~ s26, lwd = 2, ylab = 'Normalized Expression', xlab = "rs1056441", main = "miR-4745-3p",
        axes = F)
stripchart(e46 ~ s26, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("TT","TC","CC"))
axis(side=2)
dev.off()
#---------------
e47 =as.numeric(mirsnps["hsa-miR-1538",])

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs1056441_hsa-miR-1538.pdf")
boxplot(e47 ~ s26, lwd = 2, ylab = 'Normalized Expression', xlab = "rs1056441", main = "miR-1538",
        axes = F)
stripchart(e47 ~ s26, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("TT","TC","CC"))
axis(side=2)
dev.off()
#-------------------
e48 = as.numeric(mirsnps["hsa-miR-4467",])

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs1056441_hsa-miR-4467.pdf")
boxplot(e48 ~ s26, lwd = 2, ylab = 'Normalized Expression', xlab = "rs1056441", main = "miR-4467",
        axes = F)
stripchart(e48 ~ s26, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("TT","TC","CC"))
axis(side=2)
dev.off()
#--------------------
e49 = as.numeric(mirsnps["hsa-miR-6770-3p",])

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs1056441_hsa-miR-6770-3p.pdf")
boxplot(e49 ~ s26, lwd = 2, ylab = 'Normalized Expression', xlab = "rs1056441", main = "miR-6770-3p",
        axes = F)
stripchart(e49 ~ s26, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("TT","TC","CC"))
axis(side=2)
dev.off()
#---------------------
e50 = as.numeric(mirsnps["hsa-miR-3940-3p",])

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/differential/mirsnp_rs1056441_hsa-miR-3940-3p.pdf")
boxplot(e50 ~ s26, lwd = 2, ylab = 'Normalized Expression', xlab = "rs1056441", main = "miR-3940-3p",
        axes = F, ylim = c(0, 80))
stripchart(e50 ~ s26, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("TT","TC","CC"))
axis(side=2)
dev.off()
#----------------------
e51 = as.numeric(mirsnps["hsa-miR-762",])

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs1056441_hsa-miR-762.pdf")
boxplot(e51 ~ s26, lwd = 2, ylab = 'Normalized Expression', xlab = "rs1056441", main = "miR-762",
        axes = F)
stripchart(e51 ~ s26, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("TT","TC","CC"))
axis(side=2)
dev.off()
###############################################
#miRSNP rs7444 ####################
snps$rs7444[snps$rs7444 == "TT"] <- "0"
snps$rs7444[snps$rs7444 == "CT"] <- "1"
snps$rs7444[snps$rs7444 == "CC"] <- "2"

e52 = as.numeric(mirsnps[c("hsa-miR-4741"),])
s27 = as.numeric(snps$rs7444)

lm27 = lm(e52 ~ s27)
summary(lm27)

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs7444_hsa-miR-4741.pdf")
boxplot(e52 ~ s27, lwd = 2, ylab = 'Normalized Expression', xlab = "rs7444", main = "miR-4741",
        axes = F, ylim = c(0,20))
stripchart(e52 ~ s27, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21, ylim = c(0,20))
axis(1,at=c(1:3),labels=c("TT","CT","CC"))
axis(side=2)
dev.off()
#-------------
e53 = as.numeric(mirsnps[c("hsa-miR-4763-3p"),])

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs7444_hsa-miR-4763_3p.pdf")
boxplot(e53 ~ s27, lwd = 2, ylab = 'Normalized Expression', xlab = "rs7444", main = "miR-4763-3p",
        axes = F)
stripchart(e53 ~ s27, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("TT","CT","CC"))
axis(side=2)
dev.off()
#-------------
e54 =  as.numeric(mirsnps[c("hsa-miR-3918"),])

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs7444_hsa-miR-3918.pdf")
boxplot(e54 ~ s27, lwd = 2, ylab = 'Normalized Expression', xlab = "rs7444", main = "miR-3918",
        axes = F)
stripchart(e54 ~ s27, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("TT","CT","CC"))
axis(side=2)
dev.off()
#-----------
e55 =  as.numeric(mirsnps[c("hsa-miR-1207-5p"),])

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs7444_hsa-miR-1207_5p.pdf")
boxplot(e55 ~ s27, lwd = 2, ylab = 'Normalized Expression', xlab = "rs7444", main = "miR-1207-5p",
        axes = F)
stripchart(e55 ~ s27, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("TT","CT","CC"))
axis(side=2)
dev.off()
#####################
#miRSNP rs7445 #######
snps$rs7445[snps$rs7445 == "CC"] <- "0"
snps$rs7445[snps$rs7445 == "TC"] <- "1"
snps$rs7445[snps$rs7445 == "TT"] <- "2"

e56 =  as.numeric(mirsnps[c("hsa-miR-3064-5p"),])
s28 = as.numeric(snps$rs7445)

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/mirsnp_rs7445_hsa-miR-3064_5p.pdf")
boxplot(e56 ~ s28, lwd = 2, ylab = 'Normalized Expression', xlab = "rs7445", main = "miR-3064-5p",
        axes = F)
stripchart(e56 ~ s28, vertical = TRUE,  
           method = "jitter", add = TRUE, col = adjustcolor("blue", alpha.f = 0.3), 
           bg=adjustcolor("cyan", alpha.f = 0.3), pch=21)
axis(1,at=c(1:3),labels=c("CC","TC","TT"))
axis(side=2)
dev.off()

####################################################################################
##mRNA data##########
e5 = as.numeric(rna_snp1[c("ENSG00000185651"),])

e5 = log2(e5)

##Plot the data
plot(e1~ jitter(s1),
     col=(s1+1), xaxt="n",xlab="rs7444",ylab="Expression", main = "hsa-miR-4741")
axis(1,at=c(0:2),labels=c("TT","CT","CC"))

#lines(lm1$fitted ~ s1,type="b",pch=15,col="darkgrey")

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripchart_rs744_miRNAs.pdf")
stripchart(e1 ~ s1, vertical = TRUE, method = "jitter", pch = 1, add = F , col = "lightgray", 
           axes = F, lwd=2, ylab = "Expression", xlab = "rs7444", main = "hsa-miR-4741", y = c(0, 25))
axis(1,at=c(1:3),labels=c("TT","CT","CC"))
axis(side=2)

stripchart(e2 ~ s1, vertical = TRUE, method = "jitter", pch = 1, add = F , col = "lightgray", 
           axes = F, lwd=2, ylab = "Expression", xlab = "rs7444", main = "hsa-miR-4763-3p")
axis(1,at=c(1:3),labels=c("TT","CT","CC"))
axis(side=2)

stripchart(e3 ~ s1, vertical = TRUE, method = "jitter", pch = 1, add = F , col = "lightblue", 
           axes = F, lwd=2, ylab = "Expression", xlab = "rs7444", main = "hsa-miR-3918")
axis(1,at=c(1:3),labels=c("TT","CT","CC"))
axis(side=2)

stripchart(e4 ~ s1, vertical = TRUE, method = "jitter", pch = 1, add = F , col = "lightblue", 
           axes = F, lwd=2, ylab = "Expression", xlab = "rs7444", main = "hsa-miR-1207-5p")
axis(1,at=c(1:3),labels=c("TT","CT","CC"))
axis(side=2)

dev.off()
###
tiff("~/Dropbox/shared_immune_mirsnps/figs/stripchart_rs7444_miR4741.tiff")
stripchart(e1 ~ s1, vertical = TRUE, method = "jitter", pch = 1, add = F , col = "lightgray", 
           axes = F, lwd=2, ylab = "Expression", xlab = "rs7444", main = "hsa-miR-4741", y = c(0, 25))

axis(1,at=c(1:3),labels=c("TT","CT","CC"))
axis(side=2)
dev.off()

#boxplot(e1 ~ s1,axes = F, outline = F, add = F,  ylab = "Expression", xlab = "rs7444", main = "hsa-miR-4741", y = c(0, 25))
#boxplot(e1 ~ s1,axes = F, outline = F, add = T,col = "lightgray")

##mRNA
tiff("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/Geuvadis_UBE2L3_rs7444.tiff")
stripchart(e5 ~ s1, vertical = TRUE, method = "jitter", pch = 1, add = F , col = "lightblue", 
           axes = F, lwd=2, ylab = "Expression", xlab = "rs7444", main = "UBE2L3")
axis(1,at=c(1:3),labels=c("TT","CT","CC"))
axis(side=2)
dev.off()

pdf("~/Dropbox/shared_immune_mirsnps/figs/stripcharts/Geuvadis_UBE2L3_rs7444.pdf")
stripchart(e5 ~ s1, vertical = TRUE, method = "jitter", pch = 1, add = F , col = "lightblue", 
           axes = F, lwd=2, ylab = "Expression", xlab = "rs7444", main = "UBE2L3")
axis(1,at=c(1:3),labels=c("TT","CT","CC"))
axis(side=2)
dev.off()




###Correlation plots####
snp1 <- read.table ('~/rodrigo/resources/geuvadis_1000G/rs7444.ped', stringsAsFactors = F, row.names = 1)
colnames(snp1)
snp1$SNP1 <- paste(snp1$V7,snp1$V8)

snp1 <- snp1[,c(7,8)]
head(snp1)
snp1$V8 <- NULL


miR4741 <- mir_snp1[c("hsa-miR-4741"),]
miR4741 <- as.matrix(miR4741)

ubel3 <-rna_snp1[c('ENSG00000185651'), ]
ubel3 <- as.matrix(ubel3)

plot(miR4741,ubel3, xlim = c(0, 25))

cor.test(miR4741,ubel3)

snp1_hm <- snp1[snp1$SNP1 == 'C C', ]
snp1_hm <- as.data.frame(snp1_hm)
snp1_hm <- snp1_hm[which(rownames(snp1) %in% rownames(miR4741))]

snp1[which(rownames(snp1) %in% colnames(mir_snp1)),]

miR4741_hm <- miR4741[,]


##Subesting SNP1################
mirsnp1 <- as.matrix(mir_df[rownames(mir_df %in% snp1_hm) == 'hsa-miR-4741'])
ube2l3_snp1 <- as.data.frame(rna_df[rownames(rna_df %in% snp1_hm) == 'ENSG00000185651'])

plot (mirsnp1, ube2l3_snp1)

plot(mir_df[rownames(mir_df) == 'hsa-miR-4741'], rna_df[rownames(rna_df) == 'ENSG00000185651'])
