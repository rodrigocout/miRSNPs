###############################################
########PreprocessII concatenate datas
###############################################

###########################
########Load
###########################
load(file="~/Dropbox/BED_files_miRSNPS/databases/mirsnptarget.RData")
load(file="~/Dropbox/BED_files_miRSNPS/databases/newmirsnpscoreREV.RData")
load(file="~/Dropbox/BED_files_miRSNPS/databases/polymirts_merge.RData")

#---mirsnpscore
nm<-rep("mirsnpscpre", nrow(dtmir))
dtmir$miRNA <- paste ("hsa-", dtmir$miRNA, sep = "")
dtmir$miRNA <- sub ("(\\*)(.)*","", dtmir$miRNA, perl=T)
bdata<-cbind(dtmir[,c(10,6,7)],nm)
colnames(bdata)<-c("SNP","miRNA","score","BD")

#---mirsnptarget
nm<-rep("mirsnptarget", nrow(mirsnptarget))
tmp<-cbind(mirsnptarget[,c(1,3,5)],nm)
colnames(tmp) <- c("SNP","miRNA","score","BD")
bdata<-rbind(bdata,tmp)

#---polymirts
nm<-rep("polymirts", nrow(polymirts_merge))
tmp <- cbind(polymirts_merge, nm)
colnames(tmp) <- c("SNP","miRNA","score", "BD")
bdata <- rbind(bdata,tmp)

#---as.numeric
bdata$score [bdata$score == "\\N"] <- NA
bdata[,"socre"] <- as.numeric(bdata[,"score"])

#---remove NAs ?

save(bdata,file="~/Dropbox/BED_files_miRSNPS/databases/bigdata.RData")
