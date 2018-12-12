#########################
### Load
#########################
load ("~/Dropbox/BED_files_miRSNPS/databases/bigdata.RData")
load ("~/Dropbox/BED_files_miRSNPS/databases/navyResultsData.RData")

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
### Consisten test
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

#---
bdata2 <- bdata2[,-1]
bdata2 [bdata2 != 0] <- 1

#--links 2d3
bdata2$Sum <- rowSums (bdata2, na.rm = TRUE)
idx <- which (bdata2$Sum > 1)
links2d3 <- rownames (bdata2[idx,])

#---links navy
navyLinks <- rownames (navyD)

#---test
(length (intersect (links2d3, navyLinks)) / length (links2d3)) * 100
test <- which (links2d3 %in% navyLinks)
