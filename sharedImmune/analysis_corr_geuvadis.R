###mirComb analysis in GEUVADIS RNA-seq data######################
setwd('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/geuvadis_rnaseq/')
library(miRComb)

mirna <- read.delim('geuvadis_rnaseq_mirna_expression.txt', header = T)
head(mirna)
class(mirna)
mirna <- as.matrix(mirna)
mirna <- log2(mirna)
dim(mirna)


genes <- read.delim('geuvadis_rnaseq_gene_expression.txt', header = T)
head(genes)
class(genes)
genes <- as.matrix(genes)
dim(genes)

data.obj<-new("corObject",dat.miRNA=as.matrix(mirna),dat.mRNA=as.matrix(genes))

###Subset miRNAs
data.obj<-addSig(data.obj,"miRNA",manual=c("hsa-miR-4763-3p",	"hsa-miR-3918",	"hsa-miR-1207-5p",	
                                           "hsa-miR-4741",	"hsa-miR-3064-5p",	"hsa-miR-4301",	"hsa-miR-136",	"hsa-miR-3156-3p",	
                                           "hsa-miR-376a-3p",	"hsa-miR-4753-5p",	"hsa-miR-3149",	"hsa-miR-3713",	"hsa-miR-887",	"hsa-miR-191-5p",	"hsa-miR-887-3p",	
                                           "hsa-miR-6780a-3p",	"hsa-miR-4502",	"hsa-miR-6881-3p",	"hsa-miR-7111-3p",	"hsa-miR-4314",	"hsa-miR-4497",	"hsa-miR-326",	"hsa-miR-330-5p",	"hsa-miR-3649",	
                                           "hsa-miR-4518",	"hsa-miR-3192",	"hsa-miR-518c",	"hsa-miR-5001-3p",	"hsa-miR-4531",	"hsa-miR-20a",	"hsa-miR-217",	"hsa-miR-1200",	"hsa-miR-1253",	"hsa-miR-4324",	"hsa-miR-449b",	
                                           "hsa-miR-449b-3p",	"hsa-miR-4691-3p",	"hsa-miR-4802-5p",	"hsa-miR-769-5p",	"hsa-miR-513c-3p",	"hsa-miR-513a-3p",	"hsa-miR-181c",	"hsa-miR-5189-3p",	"hsa-miR-4259",	"hsa-miR-4635",	
                                           "hsa-miR-4793-5p",	"hsa-miR-4433b-5p",	"hsa-miR-4449",	"hsa-miR-2114",	"hsa-miR-2114-5p",	"hsa-miR-1251",	"hsa-miR-554",	"hsa-miR-517",	"hsa-miR-6738-3p",	"hsa-miR-544b",	"hsa-miR-4682",	"hsa-miR-4680-3p",	
                                           "hsa-miR-3140-5p",	"hsa-miR-7155-3p",	"hsa-miR-3136-3p",	"hsa-miR-588",	"hsa-miR-1245a",	"hsa-miR-8079",	"hsa-miR-5587-5p",	"hsa-miR-649",	"hsa-miR-10b-5p",	"hsa-miR-10a-5p",	"hsa-miR-3620-3p",	"hsa-miR-6865-3p",	
                                           "hsa-miR-6879-3p",	"hsa-miR-4748",	"hsa-miR-4464",	"hsa-miR-6839-3p"))

data.obj <- addSig(data.obj, "mRNA", manual = c("ENSG00000185651",	"ENSG00000115607",	"ENSG00000181634",	
                                                "ENSG00000168496",	"ENSG00000176920",	"ENSG00000161405",	"ENSG00000111321" ,	"ENSG00000175354" ,
                                                "ENSG00000185650",	"ENSG00000184293" ,	"ENSG00000176142",	"ENSG00000121594",	"ENSG00000150637",	"ENSG00000123427",	
                                                "ENSG00000123297",	"ENSG00000197991",	"ENSG00000145626" ,	"ENSG00000128604",	"ENSG00000135407" ,	"ENSG00000121797",	"ENSG00000110848",	
                                                "ENSG00000160013"))

###Correlation#########
data.obj<-addCorrelation(data.obj,alternative="less")
data.obj@cor[1:3,1:3]
data.obj@pval[1:3,1:3]
class(data.obj@cor)

plotCorrelation(data.obj,miRNA="hsa-miR-191-5p",mRNA="ENSG00000111321",type="cor",
                col.color="group",sample.names=TRUE)

plotCorrelation(data.obj,miRNA="hsa-miR-4763-3p",mRNA="ENSG00000185651",type="cor",
                col.color="group",sample.names=FALSE)

##Organize the pairs in rows
data.obj<-addNet(data.obj)
head(data.obj@net)

#################
data(microCosm_v5_18)
data(targetScan_v6.2_18)

##To use this you have to convert genes ID
data.obj<-addDatabase(data.obj,database=c("microCosm_v5_18","targetScan_v6.2_18"))
head(data.obj@net)

#####################
##Save results
writeCsv(data.obj,"geuvadis_mirna_mrna_mircomb.txt")
write.table(data.obj@net,"autoimmune_geuvadis_mirna_mrna_pearson_cor_mircomb.txt", sep = '\t', quote = F, row.names = F)

##Network
plotNetwork(data.obj,dat.sum=1)


##########################################################
#Plot
mir_191 <- mirna[c('hsa-miR-191-5p'),]
mir_191 <- as.matrix(mir_191)

fut2 <- genes[c('ENSG00000176920'),]
fut2 <- as.matrix(fut2)

plot(log2(fut2), log2(mir_191))
plot(mir_191, fut2)

plot(mir_191,fut2, main = 'FUT2 and miR-191-5p', ylab = c('FUT2 Norm Expression'), xlab = c('miR-191-5p Norm Expression'), 
     pch = 19, lwd=1, bty='n', abline(lm(mir_191 ~ fut2)))

plot(fut2,mir_191, main = 'FUT2 and miR-191-5p', ylab = c('FUT2 Norm Expression'), xlab = c('miR-191-5p Norm Expression'), 
     pch = 19, lwd=1, bty='n', abline(lm(fut2 ~ mir_191)))


# xlim = c(4.54, 4.99), ylim = c(4.6, 5),

#abline(lm(I(y-q) ~ I(x-p) + 0, data=test), col="red")

cor(mir_191, fut2, method = "spearman", use = "pairwise")

cor(mir_191, fut2, method = "pearson", use = "pairwise")
cor.test(mir_191, fut2, method = "pearson")
cor.test(mir_191, fut2, method = "spearman", use = "pairwise")

miR_449b <- mirna[c("hsa-miR-449b-3p"),]
miR_449b <- as.matrix(miR_449b)

tmem <- genes[c("ENSG00000176142"),]
tmem <- as.matrix(tmem)

plot(miR_449b, tmem)
cor.test(miR_449b, tmem, method = "pearson")

mir10 <- mirna[c('hsa-miR-10a-5p'),]
mir10 <- as.matrix(mir10)

ptpn2 <- genes[c('ENSG00000175354'),]
ptpn2 <- as.matrix(ptpn2)

plot(log2(ptpn2), log2(mir10))
plot(ptpn2, mir10)

plot(log2(mir10),log2(ptpn2), main = 'PTPN2 and miR-10a-5p', ylab = c('PTPN2 Norm Expression'), xlab = c('miR-10-5p Norm Expression'), 
     pch = 19, lwd=1, bty='n', xlim = c(10, 19), abline(lm( ptpn2 ~ mir10)))

plot(mir10,ptpn2, main = 'PTPN2 and miR-10a-5p', ylab = c('PTPN2 Norm Expression'), xlab = c('miR-10a-5p Norm Expression'), 
     pch = 19, lwd=1, bty='n', xlim = c(10, 19), abline(lm( ptpn2 ~ mir10)))

cor(mir10, ptpn2, method = "pearson", use = "pairwise")
cor.test(mir10, ptpn2, method = "pearson")




#####################################################################
#Subset
sub_mirnas <- mirna[c("hsa-miR-4763-3p",	"hsa-miR-3918",	"hsa-miR-1207-5p",	
                      "hsa-miR-4741","hsa-miR-3064-5p","hsa-miR-4301","hsa-miR-376a-3p",
                      "hsa-miR-4753-5p",	"hsa-miR-3149",	"hsa-miR-3713",	"hsa-miR-887"),]
                      "hsa-miR-191-5p"),]	
                      "hsa-miR-887-3p"),]
                      "hsa-miR-4502"),]
"hsa-miR-6881-3p",	"hsa-miR-7111-3p",	"hsa-miR-4314",	"hsa-miR-4497",	"hsa-miR-326",	"hsa-miR-330-5p",	"hsa-miR-3649"
"hsa-miR-4518",	"hsa-miR-3192",	"hsa-miR-518c",	"hsa-miR-5001-3p",	"hsa-miR-4531",	"hsa-miR-20a",	"hsa-miR-217",	"hsa-miR-1200",	"hsa-miR-1253",	"hsa-miR-4324",	"hsa-miR-449b",	
"hsa-miR-449b-3p",	"hsa-miR-4691-3p",	"hsa-miR-4802-5p",	"hsa-miR-769-5p",	"hsa-miR-513c-3p",	"hsa-miR-513a-3p",	"hsa-miR-181c",	"hsa-miR-5189-3p",	"hsa-miR-4259",	"hsa-miR-4635",	
"hsa-miR-4793-5p",	"hsa-miR-4433b-5p",	"hsa-miR-4449",	"hsa-miR-2114",	"hsa-miR-2114-5p",	"hsa-miR-1251",	"hsa-miR-554",	"hsa-miR-517",	"hsa-miR-6738-3p",	"hsa-miR-544b",	"hsa-miR-4682",	"hsa-miR-4680-3p",	
"hsa-miR-3140-5p",	"hsa-miR-7155-3p",	"hsa-miR-3136-3p",	"hsa-miR-588",	"hsa-miR-1245a",	"hsa-miR-8079",	"hsa-miR-5587-5p",	"hsa-miR-649",	"hsa-miR-10b-5p",	"hsa-miR-10a-5p",	"hsa-miR-3620-3p",	"hsa-miR-6865-3p",	
"hsa-miR-6879-3p",	"hsa-miR-4748",	"hsa-miR-4464",	"hsa-miR-6839-3p"),]





"hsa-miR-6780a-3p")
"hsa-miR-3156-3p"),]	
"hsa-miR-136"),]

                      