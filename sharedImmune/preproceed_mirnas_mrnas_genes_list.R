#####Proceed data analysis#############

targets <- read.table("~/Dropbox/shared_immune_mirsnps/tables/results/Lista_target_genes_mirsnps.txt")
#---get annotation
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- targets$V1

gene_id <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
head(gene_id)

lista_mirs <- read.table("~/Dropbox/shared_immune_mirsnps/tables/results/lista_microRNAs.txt")
class(lista_mirs)

#########microRNAS#######################################
m1 <- read.delim("~/Dropbox/shared_immune_mirsnps/data/geuvadis_miRNAS_matched_rna.txt")

lista2 <- intersect(rownames(m1), lista_mirs$V1)
class(lista2)

m1 <- m1[lista2,]

all(rownames(m1) %in% lista2)
all(rownames(m1) == lista2)

###60 miRNAS expressos in GEUVADIS##########################
write.table(m1, "~/Dropbox/shared_immune_mirsnps/data/microRNAS_in_Geuvadis.txt", sep = '\t', quote = F)


###mRNA##################################
r1 <- read.delim("~/Dropbox/shared_immune_mirsnps/data/geuvadis_RNAS_matched_mirna.txt")

lista3 <- intersect(rownames(r1), gene_id$ensembl_gene_id)

r1 <- r1[lista3,]


r1_df <- r1[which(rownames(r1) %in% gene_id$ensembl_gene_id),]

r1_df <- r1_df[order(rownames(r1_df)),]

id <- gene_id[order(gene_id[,1]),]
id <- id[which(id[,1] %in% rownames(r1_df)),]

all(id[,1] == rownames(r1_df))
all(id[,1] %in% rownames(r1_df))

rownames(r1_df) <- id[,2]

write.table(r1_df, "~/Dropbox/shared_immune_mirsnps/data/Gene_names_mirsnps_geuvadis.txt", sep = '\t', quote = F)

r1 <- read.delim("~/Dropbox/shared_immune_mirsnps/data/Gene_names_mirsnps_geuvadis.txt")
m1 <- read.delim("~/Dropbox/shared_immune_mirsnps/data/micronas_snps_geuvadis.txt")

manba <- r1[c("MANBA"),]
manba <- log2(manba)

miR660_5p <- m1[c("hsa-miR-660-5p"),]
miR660_5p <- log2(miR660_5p)

pdf("~/Dropbox/shared_immune_mirsnps/figs/Correlation_plot_mir660_manba.pdf")
plot( miR660_5p, manba)
dev.off()

par("mar")
par(mar=c(1,1,1,1))

par(mar = rep(2, 4))
par(mfrow=c(4,2))
