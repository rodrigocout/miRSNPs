#########################
### Load
#########################
load ("~/Dropbox/BED_files_miRSNPS/databases/bigdata.RData")

#########################
### pre-Process
#########################
bdata <- bdata[,-3]

#---bdata2
bdata$miRSNP <- paste (bdata$SNP, bdata$miRNA, sep = "~")
bdata2 <- data.frame (miRSNP = unique (bdata$miRSNP), mirsnpscpre = NA,
                      mirsnptarget = NA, polymirts = NA, stringsAsFactors = FALSE)
rownames (bdata2) <- bdata2$miRSNP

#########################
### Process
#########################

#---mirsnpscore
tmp <- bdata [bdata$BD == "mirsnpscpre",]
bdata2 [tmp$miRSNP,"mirsnpscpre"] <- tmp$socre

#---mirsnptarget
tmp <- bdata [bdata$BD == "mirsnptarget",]
bdata2 [tmp$miRSNP,"mirsnptarget"] <- tmp$socre

#---polymirts
tmp <- bdata [bdata$BD == "polymirts",]
bdata2 [tmp$miRSNP,"polymirts"] <- tmp$socre

#---SUB miRNA and SNP
SNP <- sub ("(\\~)(.)*","", bdata2$miRSNP, perl=T)
miRNA <- sub ("(.)*(\\~)", "", bdata2$miRSNP, perl=T)

#---
bdata2 <- bdata2[,-1]
datanavy <- data.frame (SNP = SNP, miRNA = miRNA, mirsnpscore = bdata2$mirsnpscpre, mirsnptarget =
                          bdata2$mirsnptarget, polymirts = bdata2$polymirts, stringsAsFactors = FALSE)
rownames (datanavy) <- rownames (bdata2)

save (datanavy, file = "~/Dropbox/BED_files_miRSNPS/databases/navyData.RData")
