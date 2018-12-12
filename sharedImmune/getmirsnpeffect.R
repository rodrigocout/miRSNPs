####Get miRSNP effect####
setwd("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data")
mirsnptarget<-read.table(file = "~/Dropbox/BED_files_miRSNPS/databases/MirSNPInTarget.txt",sep="",header =F,stringsAsFactors = FALSE)
tmp<-paste("rs",mirsnptarget[,1],sep="")
mirsnptarget[,1]<-tmp
head(mirsnptarget)
mirsnptarget<-mirsnptarget[,-c(5,6,12,13,15:20)]
colnames(mirsnptarget)<-c("SNP","Gene","miRNA","transcript","alleles","effect", "allele","score", "energy","conservation")
write.table(mirsnptarget, "mirsnpintarget.txt", sep = "\t", quote = F)

mirsnp_ait <- read.delim("mirsnps.txt", header = T)
colnames(mirsnp_ait)

df <- merge(mirsnp_ait, mirsnptarget, by = 'SNP')
head(df)

write.table(df, 'mirsnp_effects.txt', sep = "\t", quote = F)

polymirts<-read.table(file = "~/Dropbox/BED_files_miRSNPS/databases/Polymirts.txt", sep="\t", quote = '',header=TRUE,stringsAsFactors = FALSE)
head(polymirts, n = 5)

mirsnpscore<-read.table(file = "~/Dropbox/BED_files_miRSNPS/databases/mirsnpscore.txt",sep="#",header=F,stringsAsFactors=FALSE)
head(mirsnpscore)
