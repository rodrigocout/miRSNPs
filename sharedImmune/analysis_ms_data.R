### Analysis miRNA mRNA in MS T cells #############################

# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.36.0, limma 3.26.8
# R scripts generated  Tue May 24 12:08:18 EDT 2016

#Server: www.ncbi.nlm.nih.gov
#Query: acc=GSE43590&platform=GPL14613&type=txt&groups=&colors=&selection=XXXXXXXXXXXXXXXXXXXX&padj=fdr&logtransform=auto&columns=ID&columns=adj.P.Val&columns=P.Value&columns=F&columns=miRNA_ID_LIST&columns=SPOT_ID&num=250&annot=submitter

# Unable to generate script analyzing differential expression.
#      Invalid input: at least two groups of samples should be selected.

################################################################
#   Boxplot for selected GEO samples
library(Biobase)

source("https://bioconductor.org/biocLite.R")
biocLite("GEOquery")
library(GEOquery)

# load series and platform data from GEO

#miRNA dataset
##MicroRNA regulate immune pathways in T-cells in multiple sclerosis (MS) miRNA
gset <- getGEO("GSE43590", GSEMatrix =TRUE)
if (length(gset) > 1) idx <- grep("GPL14613", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
show(gset)
show(pData(phenoData(gset)))
head(exprs(gset))

# names of all the GSM objects contained in the GSE
names(GSMList(gset))

#mRNA dataset
##MicroRNA regulate immune pathways in T-cells in multiple sclerosis (MS) mRNA
gset1 <- getGEO("GSE43591", GSEMatrix = TRUE)

#plataform GPL570
if (length(gset1) > 1) idx <- grep("GPL570", attr(gset1, "names")) else idx <- 1
gset1 <- gset1[[idx]]
head(exprs(gset1))






# set parameters and draw the plot

dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE43590", '/', annotation(gset), " selected samples", sep ='')
boxplot(exprs(gset), boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)
