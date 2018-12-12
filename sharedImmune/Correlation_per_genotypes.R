mirsnps <- read.delim("~/Dropbox/shared_immune_mirsnps/data/micronas_snps_geuvadis.txt")
rna_snps <- read.delim("~/Dropbox/shared_immune_mirsnps/data/mRNAs_snps_geuvadis.txt")

all(colnames(rna_snps) %in% colnames(mirsnps))
all(colnames(rna_snps) == colnames(mirsnps))

lista <- read.table("~/Dropbox/shared_immune_mirsnps/tables/results/exploratory/lista_microRNAs.txt")
head(lista)

miRInter <- intersect(rownames(mirsnps), lista$V1)

mirsnps <- mirsnps[miRInter,]
#mirsnps <- log2(mirsnps)
head(mirsnps)
rownames(mirsnps)

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

#make all the same number of samples
snps <- snps[which(rownames(snps) %in% colnames(mirsnps)),]
head(snps)

###############################################
#miRSNP rs7444 ####################
snps$rs7444[snps$rs7444 == "TT"] <- "0"
snps$rs7444[snps$rs7444 == "CT"] <- "1"
snps$rs7444[snps$rs7444 == "CC"] <- "2"

e52 = as.numeric(mirsnps[c("hsa-miR-4741"),])
s27 = as.numeric(snps$rs7444)

lm27 = lm(e52 ~ s27)
summary(lm27)

rs7444_cc <- snps[snps$rs7444 == "CC",]
indiv_cc <- row.names(rs7444_cc)
indiv_cc

indiv_tt <- snps[snps$rs7444 == "TT",]
indiv_tt <- row.names(indiv_tt)


mir4741 <- mirsnps[c("hsa-miR-4741"),]
mir4741 <- as.matrix(mir4741)
mir4741 <- log2(mir4741)
mir4741[which(mir4741== "-Inf")] = NA

ube2l3 <- rna_snps[c("UBE2L3"),]
ube2l3 <- as.matrix(ube2l3)
ube2l3 <- log2(ube2l3)

#mir4741[which(is.nan(mir4741))] = NA
ube2l3[which(ube2l3==Inf)] = NA

pdf("~/Dropbox/shared_immune_mirsnps/figs/correlation/NewScatterplot_correlation_ube2l3_mir4741.pdf")

par(mar = c(10, 4, 4, 2) + 0.1, mfrow = c(1, 2))

plot( mir4741[,indiv_cc] , ube2l3[,indiv_cc],
      pch = 19, bty = "n", ylim = c(10, 13), xlim = c(-0.5, 6),
       xlab = c("miR-4741 Log2 expression"), 
     ylab = c("UBE2L3 Log2 expression"), main = "rs7444 C/C",
     abline(lsfit( mir4741[,indiv_cc], ube2l3[,indiv_cc])))
legend(c("topright"), c("r2= 0.52", "P = 0.01687" ), bty = 'n', 
       cex=0.8)


cor.test(mir4741[,indiv_cc] , ube2l3[,indiv_cc], method = "pearson")

summary(lm( mir4741[,indiv_cc] ~ ube2l3[,indiv_cc]))

    #abline(lsfit( mir4741[,indiv_cc], ube2l3[,indiv_cc]))

plot( mir4741[,indiv_tt],ube2l3[,indiv_tt],
     pch = 19, bty = "n", 
     xlim = c(-0.5, 8), ylim = c(10, 12.5), 
     xlab = c("miR-4741 Log2 expression"), 
     ylab = c("UBE2L3 Log2 expression"), main = "rs7444 T/T", 
     abline(lsfit( mir4741[,indiv_tt], ube2l3[,indiv_tt])))
legend(c("topright"), c("P-value = 0.006", "r2 = 0.21"), bty = 'n', 
       cex=0.8)

dev.off()

summary(lm(ube2l3[,indiv_tt] ~ mir4741[,indiv_tt]))
cor.test(ube2l3[,indiv_tt], mir4741[,indiv_tt], method = "pearson")
###################################################################
###rs9943_hsa-miR-628-5p#####

rs9943_gg <- snps[snps$rs9943 == "GG",]
indiv_gg <- row.names(rs9943_gg)
indiv_gg

indiv_aa <- snps[snps$rs9943 == "AA",]
indiv_aa <- row.names(indiv_aa)

snps$rs9943[snps$rs9943 == 'AA'] <- '0'
snps$rs9943[snps$rs9943 == 'GA'] <- '1'
snps$rs9943[snps$rs9943 == 'GG'] <- '2'


miR628_5p <- mirsnps[c("hsa-miR-628-5p"),]
miR628_5p <- as.matrix(miR628_5p)
miR628_5p <- log2(miR628_5p)
miR628_5p[which(miR628_5p== "-Inf")] = NA


COG6 <- rna_snps[c("COG6"),]
COG6 <- as.matrix(COG6)
COG6 <- log2(COG6)
COG6 [which(COG6 == "-Inf")] = NA


pdf("~/Dropbox/shared_immune_mirsnps/figs/correlation/NewScaterrplot_Pearson_correlation_cog6_mir628.pdf")

tiff("~/Dropbox/shared_immune_mirsnps/figs/correlation/NewScaterrplot_Pearson_correlation_cog6_mir628.tiff", res = 300)
par(mar = c(10, 4, 4, 2) + 0.1, mfrow = c(1, 2))


plot( miR628_5p[,indiv_gg],COG6[,indiv_gg],
      pch = 19, bty = "n", xlim = c(0, 5), ylim = c(8.5,10),
      xlab = c("miR-628-5p Log2 expression"), 
      ylab = c("COG6 Log2 expression"), 
      main = "rs9943 G/G",
      abline(lsfit(miR628_5p[,indiv_gg], COG6[,indiv_gg])))
legend(c("topleft"), c("r2 = 0.43","P = 0.0169"), bty = 'n', 
       cex=0.8)


lm1 <- lm(miR628_5p[,indiv_gg] ~ COG6[,indiv_gg])
summary(lm1)

cor.test(miR628_5p[,indiv_gg], COG6[,indiv_gg], method = "pearson")

plot( miR628_5p[,indiv_aa],COG6[,indiv_aa],
      pch = 19, bty = "n", xlim = c(0, 5), ylim = c(8.5, 10),
      xlab = c("miR-628-5p Log2 expression"), 
      ylab = c("COG6 Log2 expression"), 
      main = "rs9943 A/A",
      abline(lsfit(miR628_5p[,indiv_aa],COG6[,indiv_aa])))
legend(c("topleft"), c("r2 = 0.11","P-value = 0.22"), bty = 'n', 
       cex=0.8)

cor.test(miR628_5p[,indiv_aa], COG6[,indiv_aa], method = "pearson")

dev.off()

#################################

rs907091_tt <- snps[snps$rs907091 == "TT",]
rs907091_tt <- row.names(rs907091_tt)

rs907091_cc <- snps[snps$rs907091 == "CC", ]
rs907091_cc <- row.names(rs907091_cc)

snps$rs907091[snps$rs907091 == 'CC'] <- '0'
snps$rs907091[snps$rs907091 == 'TC'] <- '1'
snps$rs907091[snps$rs907091 == 'TT'] <- '2'


miR326 <- mirsnps["hsa-miR-326",]
miR326 <- log2(miR326)
miR326 <- as.matrix(miR326)
miR326[which(miR326== "-Inf")] = NA

IKZF3 <- rna_snps["IKZF3",]
IKZF3 <- log2(IKZF3)
IKZF3 <- as.matrix(IKZF3)

#IKZF3[which(is.nan(IKZF3)),] = NA
#IKZF3[which(IKZF3==Inf)] = NA

pdf("~/Dropbox/shared_immune_mirsnps/figs/correlation/New_Pearson_correlation_ikzf3_mir326.pdf")

par(mar = c(10, 4, 4, 2) + 0.1, mfrow = c(1, 2))

plot(miR326[,rs907091_tt], IKZF3[,rs907091_tt], 
     pch = 19, bty = "n", ylim = c(5, 12), xlim = c( -1, 6) ,
     ylab = c("IKZF3 Log2 expression"), 
     xlab = c("miR-326 Log2 expression"),
     main = "rs907091 T/T",
     abline(lsfit(miR326[,rs907091_tt], IKZF3[,rs907091_tt] )))
legend(c("topright"), c("r2 = -0.35", "P = 0.006479"), bty = 'n', 
       cex=0.8)
#p-value = 0.006479


cor.test(miR326[,rs907091_tt], IKZF3[,rs907091_tt], method = "pearson", use =  "complete.obs")

summary(lm(miR326[,rs907091_tt] ~ IKZF3[,rs907091_tt] ))


plot(miR326[,rs907091_cc], IKZF3[,rs907091_cc], 
     pch = 19, bty = "n", ylim = c(5, 12), xlim = c(-1, 6), 
     ylab = c("IKZF3 Log2 expression"), 
     xlab = c("miR-326 Log2 expression"),
     main = "rs907091 C/C",
     abline(lsfit(miR326[,rs907091_cc], IKZF3[,rs907091_cc] )))
legend(c("topright"), c("r2 = -0.24","P = 0.05" ), bty = 'n', 
       cex=0.8)
dev.off()

#p-value = 0.006479


cor.test(miR326[,rs907091_cc], IKZF3[,rs907091_cc], method = "pearson", use =  "complete.obs")

summary(lm(miR326[,rs907091_cc] ~ IKZF3[,rs907091_cc] ))
#---------------------------------
#miR-4497 negative correlated with IKZF3
miR4497 <- mirsnps["hsa-miR-4497",]
miR4497 <- log2(miR4497)
miR4497 <- as.matrix(miR4497)
miR4497[which(miR4497== "-Inf")] = NA

pdf("~/Dropbox/shared_immune_mirsnps/figs/correlation/Pearson_correlation_ikzf3_mir4497.pdf")

par(mar = c(10, 4, 4, 2) + 0.1, mfrow = c(1, 2))

plot(miR4497[,rs907091_tt], IKZF3[,rs907091_tt], 
     pch = 19, bty = "n", ylim = c(8, 12), xlim = c(0, 6), 
     ylab = c("IKZF3 Log2 expression"), 
     xlab = c("miR-4497 Log2 expression"),
     main = "Individuals rs907091 TT",
     abline(lsfit(miR4497[,rs907091_tt], IKZF3[,rs907091_tt] )))
legend(c("topright"), c("P-value = 0.37", "r2 = -0.16"), bty = 'n', 
       cex=0.8)

cor.test(miR4497[,rs907091_tt], IKZF3[,rs907091_tt], method = "pearson", use =  "complete.obs")

summary(lm(miR4497[,rs907091_tt] ~ IKZF3[,rs907091_tt] ))


plot(miR4497[,rs907091_cc], IKZF3[,rs907091_cc], 
     pch = 19, bty = "n", ylim = c(8, 12), xlim = c(0, 6), 
     ylab = c("IKZF3 Log2 expression"), 
     xlab = c("miR-4497 Log2 expression"),
     main = "Individuals rs907091 CC",
     abline(lsfit(miR4497[,rs907091_cc], IKZF3[,rs907091_cc] )))
legend(c("topright"), c("P-value = 0.97", "r2 = 0.004"), bty = 'n', 
       cex=0.8)
dev.off()

cor.test(miR4497[,rs907091_cc], IKZF3[,rs907091_cc], method = "pearson", use =  "complete.obs")

summary(lm(miR4497[,rs907091_cc] ~ IKZF3[,rs907091_cc] ))

#-------------------------------------
#miR???4518 postive correlated with IKZF3
miR4518 <- mirsnps["hsa-miR-4518",]
miR4518 <- log2(miR4518)
miR4518 <- as.matrix(miR4518)
miR4518[which(miR4518== "-Inf")] = NA

pdf("~/Dropbox/shared_immune_mirsnps/figs/correlation/Pearson_correlation_ikzf3_mir4518.pdf")

par(mar = c(10, 4, 4, 2) + 0.1, mfrow = c(1, 2))

plot(miR4518[,rs907091_tt], IKZF3[,rs907091_tt], 
     pch = 19, bty = "n", ylim = c(8, 12), xlim = c(0, 6), 
     ylab = c("IKZF3 Log2 expression"), 
     xlab = c("miR-4518 Log2 expression"),
     main = "Individuals rs907091 TT",
     abline(lsfit(miR4518[,rs907091_tt], IKZF3[,rs907091_tt] )))
legend(c("topright"), c("P-value = 0.02", "r2 = 0.2"), bty = 'n', 
       cex=0.8)

cor.test(miR4518[,rs907091_tt], IKZF3[,rs907091_tt], method = "pearson", use =  "complete.obs")

summary(lm(miR4518[,rs907091_tt] ~ IKZF3[,rs907091_tt] ))


plot(miR4518[,rs907091_cc], IKZF3[,rs907091_cc], 
     pch = 19, bty = "n", ylim = c(8, 12), xlim = c(0, 6), 
     ylab = c("IKZF3 Log2 expression"), 
     xlab = c("miR-4518 Log2 expression"),
     main = "Individuals rs907091 CC",
     abline(lsfit(miR4518[,rs907091_cc], IKZF3[,rs907091_cc] )))
legend(c("topright"), c("P-value = 0.69", "r2 = 0.04"), bty = 'n', 
       cex=0.8)
dev.off()

cor.test(miR4518[,rs907091_cc], IKZF3[,rs907091_cc], method = "pearson", use =  "complete.obs")

summary(lm(miR4518[,rs907091_cc] ~ IKZF3[,rs907091_cc] ))

#----------------------------------
#hsa-miR-3661 negative correlated with SNAPC4##
snps$rs3088081[snps$rs3088081 == 'AA'] <- '0'
snps$rs3088081[snps$rs3088081 == 'GA'] <- '1'
snps$rs3088081[snps$rs3088081 == 'GG'] <- '2'

rs3088081_aa <- snps[snps$rs3088081 == 'AA',]
rs3088081_aa <- rownames(rs3088081_aa)

rs3088081_gg <- snps[snps$rs3088081 == 'GG',]
rs3088081_gg <- rownames(rs3088081_gg)

miR3661 <- mirsnps[c("hsa-miR-3661"),]
miR3661 <- log2(miR3661)
miR3661 <- as.matrix(miR3661)
miR3661[which(miR3661== "-Inf")] = NA

SNAPC4 <- rna_snps["SNAPC4",]
SNAPC4 <- log2(SNAPC4)
SNAPC4 <- as.matrix(SNAPC4)


pdf("~/Dropbox/shared_immune_mirsnps/figs/correlation/Pearson_correlation_scatterplot_snapc4_mir3661.pdf")

par(mar = c(10, 4, 4, 2) + 0.1, mfrow = c(1, 2))

plot(miR3661[,rs3088081_aa], SNAPC4[,rs3088081_aa], 
     pch = 19, bty = "n",
     ylab = c("SNAPC4 Log2 expression"), 
     xlab = c("miR-3661 Log2 expression"),
     main = "Individuals rs30881 AA",
     abline(lsfit(miR3661[,rs3088081_aa], SNAPC4[,rs3088081_aa] )))
legend(c("topright"), c("P-value = 0.21", "r2 = -0.11"), bty = 'n', 
       cex=0.8)

cor.test(miR3661[,rs3088081_aa], SNAPC4[,rs3088081_aa], method = "pearson", use =  "complete.obs")

summary(lm(miR3661[,rs3088081_aa] ~ SNAPC4[,rs3088081_aa] ))


plot(miR3661[,rs3088081_gg], SNAPC4[,rs3088081_gg], 
     pch = 19, bty = "n", 
     ylab = c("SNPAPC4 Log2 expression"), 
     xlab = c("miR-3661 Log2 expression"),
     main = "Individuals rs3088081 GG",
     abline(lsfit(miR3661[,rs3088081_gg], SNAPC4[,rs3088081_gg] )))
legend(c("topright"), c("P-value = 0.7", "r2 = -0.05"), bty = 'n', 
       cex=0.8)
dev.off()

cor.test(miR3661[,rs3088081_gg], SNAPC4[,rs3088081_gg], method = "pearson", use =  "complete.obs")

summary(lm(miR3661[,rs3088081_gg] ~ SNAPC4[,rs3088081_gg]))

#------------------------------------
##miR-4463 positive correlated with FUT2##

rs570794_tt <- row.names(snps[snps$rs570794 == "TT",])
rs570794_ct <- row.names(snps[snps$rs570794 == "CT",])
rs570794_cc <- row.names(snps[snps$rs570794 == "CC",])


miR4463 <- mirsnps[c("hsa-miR-4463"),]
miR4463 <- log2(miR4463)
miR4463 <- as.matrix(miR4463)
miR4463[which(miR4463== "-Inf")] = NA

FUT2 <- rna_snps["FUT2",]
FUT2 <- log2(FUT2)
FUT2 <- as.matrix(FUT2)
FUT2[which(FUT2 == "-Inf")] = NA

pdf("~/Dropbox/shared_immune_mirsnps/figs/correlation/Pearson_correlation_scatterplot_fut2_mir4463.pdf")

par(mar = c(10, 4, 4, 2) + 0.1, mfrow = c(1, 2))

plot(miR4463[,rs570794_tt], FUT2[,rs570794_tt], 
     pch = 19, bty = "n",
     ylab = c("FUT2 Log2 expression"), 
     xlab = c("miR-4463 Log2 expression"),
     main = "Individuals rs570797 TT",
     abline(lsfit(miR4463[,rs570794_tt], FUT2[,rs570794_tt] )))
legend(c("topleft"), c("P-value = 0.26", "r2 = -0.14"), bty = 'n', 
       cex=0.8)

cor.test(miR4463[,rs570794_tt], FUT2[,rs570794_tt], method = "pearson", use =  "complete.obs")

summary(lm(miR4463[,rs570794_tt] ~ FUT2[,rs570794_tt] ))


plot(miR4463[,rs570794_cc], FUT2[,rs570794_cc], 
     pch = 19, bty = "n", 
     ylab = c("FUT2 Log2 expression"), 
     xlab = c("miR-4463 Log2 expression"),
     main = "Individuals rs5707941 CC",
     abline(lsfit(miR4463[,rs570794_cc], FUT2[,rs570794_cc] )))
legend(c(1,0), c("P-value = 0.19", "r2 = -0.2"), bty = 'n', 
       cex=0.8)
dev.off()

cor.test(miR4463[,rs570794_cc], FUT2[,rs570794_cc], method = "pearson", use =  "complete.obs")

summary(lm(miR4463[,rs570794_cc] ~ FUT2[,rs570794_cc]))

##--------------------------------------
#miR???455???3p negative correlated with IL18RAP
#SNP not present in GEUVADIS#####

#------------------------------------------
#miR-660-5p positive correlated with MANBA 

rs1054037_cc <- rownames(snps[snps$rs1054037 == 'CC',])
rs1054037_tt <- rownames(snps[snps$rs1054037 == 'TT',])

miR660 <- mirsnps[c("hsa-miR-660-5p"),]
miR660 <- log2(miR660)
miR660 <- as.matrix(miR660)
miR660[which(miR660== "-Inf")] = NA

MANBA <- rna_snps["MANBA",]
MANBA <- log2(MANBA)
MANBA <- as.matrix(MANBA)
MANBA[which(MANBA == "-Inf")] = NA

pdf("~/Dropbox/shared_immune_mirsnps/figs/correlation/Pearson_correlation_scatterplot_manba_mir660_5p.pdf")

par(mar = c(10, 4, 4, 2) + 0.1, mfrow = c(1, 2))

plot(miR660[,rs1054037_cc], MANBA[,rs1054037_cc], 
     pch = 19, bty = "n", ylim = c(8.5, 11),
     ylab = c("MANBA Log2 expression"), 
     xlab = c("miR-660-5p Log2 expression"),
     main = "Individuals rs1054037 CC",
     abline(lsfit(miR660[,rs1054037_cc], MANBA[,rs1054037_cc])))
legend(c("topleft"), c("P-value = 0.004", "r2 = 0.32"), bty = 'n', 
       cex=0.8)

cor.test(miR660[,rs1054037_cc], MANBA[,rs1054037_cc], method = "pearson", use =  "complete.obs")

summary(lm(miR660[,rs1054037_cc] ~ MANBA[,rs1054037_cc] ))


plot(miR660[,rs1054037_tt], MANBA[,rs1054037_tt], 
     pch = 19, bty = "n",ylim = c(8.5, 11), 
     ylab = c("MANBA Log2 expression"), 
     xlab = c("miR-660-5p Log2 expression"),
     main = "Individuals rs1054037 TT",
     abline(lsfit(miR660[,rs1054037_tt], MANBA[,rs1054037_tt] )))
legend(c("topleft"), c("P-value = 0.02", "r2 = 0.25"), bty = 'n', 
       cex=0.8)
dev.off()

cor.test(miR660[,rs1054037_tt], MANBA[,rs1054037_tt], method = "pearson", use =  "complete.obs")

summary(lm(miR660[,rs1054037_cc] ~ MANBA[,rs1054037_cc]))


