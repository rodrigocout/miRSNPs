#################
######FuncSNP###
##############
#setwd("C:/Users/max/Documents/dados_R/")

options(width=80)
library(FunciSNP)
package.version("FunciSNP")
biocLite("FunciSNP.data")


library("FunciSNP.data")

library("Rsamtools")
library("rtracklayer")
library("GGtools")
library("methods")
library("ChIPpeakAnno")
library("GenomicRanges")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("VariantAnnotation")


##Criar o dataframe para ler os SNPs
fcn1.snp <- file.path("/Users/Rodrigo/Documents/Pos_doc/IdeaS/fine_mapping/exploratory", dir(pattern='.snp$'))

FCN1.snp <- read.delim(file = fcn1.snp, sep = "",header = F) 
# Digita o sep para separar os espacos da coluna 
#  e colocar o nome das colunas
head(FCN1.snp)

##glioma.snp <- file.path(system.file('extdata', package='FunciSNP'), dir(system.file('extdata',package='FunciSNP'), pattern='.snp$'));
##gsnp <- read.delim(file=glioma.snp,sep=" ",header=FALSE);
#gsnp

FCN1.bio <- ("/Users/Rodrigo/Documents/Pos_doc/IdeaS/fine_mapping/raw_data/BEDs")
FCN1.bio
as.matrix(list.files(FCN1.bio, pattern=".bed$"))

##Funci
FCN1 <- getFSNPs(snp.regions.file=FCN1.snp, bio.features.loc=FCN1.bio)
FCN1
head(FCN1)





#head(eQTLsFCN1)



