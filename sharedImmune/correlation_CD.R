setwd('/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/raw_data/cd_mirna_mrna_rnaseq/')
library(miRComb)

################################################################################
matrix1 <- read.table ("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/raw_data/cd_mirna_mrna_rnaseq/reajuda/matrix_test_mirna.txt", header = T, stringsAsFactors = F, sep = '\t')
rownames(matrix1) <- matrix1[,1]
matrix1 <- matrix1[,-1]
matrix1 <- log2 (matrix1)
head(matrix1)
dim(matrix1)

###############################
####matrix2
###############################
matrix2 <- read.csv ("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/raw_data/cd_mirna_mrna_rnaseq/reajuda/GSE66208_Supplemental_Table_8_smRNAseq_ExpressionData.csv", header = T, sep = ',')
rownames(matrix2) <- matrix2[,1]
matrix2 <- matrix2[,-1]
matrix2 <- log2 (matrix2)

###############################
####matrix3
###############################
matrix3 <- read.csv ("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/raw_data/cd_mirna_mrna_rnaseq/reajuda/GSE66208_Supplemental_Table_8_smRNAseq_ExpressionData.csv", header = F, skip = 1, sep = ',')
rownames(matrix3) <- matrix3[,1]
matrix3 <- matrix3[,-1]
matrix3 <- log2 (matrix3)
head(matrix3)


###########################################################################################
matrix1 <- read.table ("matrix_mirna.txt", header = T, stringsAsFactors = F, sep = '\t')
rownames(matrix1) <- matrix1[,1]
matrix1 <- matrix1[,-1]
matrix1 <- log2 (matrix1)
head(matrix1)
dim(matrix1)


mirna <- read.table('matrix_mirna.txt', header = T, stringsAsFactors = F, sep = '\t')
rownames(mirna) <- mirna[,1]
mirna <- mirna[,-1]
mirna <- log2 (mirna)

head(mirna)
class(mirna)
colnames(mirna)
row.names(mirna)
dim(mirna)
str(mirna)

genes <- read.csv ("Normalized_AllCD_AllControl_GEO.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)
rownames(genes) <- genes[,1]
genes <- genes[,-1]
genes <- log2 (genes)
head(genes)
dim(genes)

colnames(genes)

all(colnames(genes) %in% colnames(mirna))

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

##Plot
mir_191 <- mirna[c('hsa-miR-191-5p'),]
mir_191 <- as.matrix(mir_191)

fut2 <- genes[c('ENSG00000176920'),]
fut2 <- as.matrix(fut2)



plot(mirna,genes, main = '', ylab = c('FUT2 Norm Expression'), xlab = c('miR-191-5p Norm Expression'), 
     pch = 19, lwd=1, bty='n', abline(lm(mir_191 ~ fut2)))

