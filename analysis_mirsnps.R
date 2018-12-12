###Analysis miRSNPs in PF #####################
setwd("/Users/lgmh/Documents/rodrigo/Projeto_GWAS_PF/gwas_pf/mirsnps_assoc/")

load("/Users/lgmh/Dropbox/BED_files_miRSNPS/databases/navyResultsData.RData")
head(navyD)
length(unique(navyD$SNP))
navySNPs <- as.data.frame(unique(navyD$SNP))

write.table(navySNPs, "/Users/lgmh/Documents/rodrigo/Projeto_GWAS_PF/gwas_pf/mirsnps_assoc/NavyBayesSNPS.txt",  quote = F , row.names = F, col.names = F)

mirsnps <- read.table("", header = T)
mirsnps <- mirsnps[mirsnps$TEST == "ADD",]
mirsnps <- mirsnps[order(mirsnps$P),]
mirsnps <- subset(mirsnps, select = c( 'CHR' , 'SNP', 'BP', 'P' ))
mirsnps <- mirsnps[complete.cases(mirsnps),]


##MH plot#########

pdf("mh_mirsnps_pf.pdf")
manhattan(mirsnps, main = "miRSNPs in Pemphigus",  cex = 0.9, 
          cex.axis = 0.9,col = c("blue4", "lightskyblue"), suggestiveline = -log10(1e-03), genomewideline = -log10(1e-05))
dev.off()

#Creating QQ-Plot
library(qqman)
qq(mirsnps$P)

qq(mirsnps$P, main = "Q-Q plot of GWAS p-values", xlim = c(0, 7), ylim = c(0, 12), 
   pch = 18, col = "blue4", cex = 1.5, las = 1)
