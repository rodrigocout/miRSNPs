####Get the mRNA data of CD samples from GEO################
library(GEOquery)
gse <- getGEO("GSE66207", GSEMatrix =TRUE)
if (length(gse) > 1) idx <- grep("GPL16791", attr(gse, "names")) else idx <- 1
gse <- gse[[idx]]
show(gse)
exp <- exprs(gse)
#show(pData(phenoData(gset)))
#head(exprs(gset), n = 40)
#exprs(gset)

gsmplatforms <- lapply(GSMList(gse),function(x) {Meta(x)$platform})
head(gsmplatforms)

