setwd('/Users/Rodrigo/Documents/Pos_doc/Projeto_epistasis_interacao/Circus/')
library(RCircos)
data(UCSC.HG19.Human.CytoBandIdeogram)
data(RCircos.Line.Data)

link.data <- read.delim('Circus_link.txt')
snps.label <- read.delim('SNPs_circus.txt')
gene.label <- read.delim('Genes_circus.txt')

chr.exclude <- c("chrX", "chrY", "chr1", "chr2", "chr4", "chr5", "chr6",
                 "chr7", "chr8", "chr9", "chr10", "chr11", "chr13", "chr14",'chr15', "chr16", "chr18", "chr20", "chr21")

cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
tracks.inside <- 5
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude,
                            tracks.inside, tracks.outside)

pdf(file="RCircos_Interactions_comSNPs.pdf", height=8, width=8)

RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()

RCircos.Gene.Name.Plot(gene.data=gene.label, name.col=4, 
                       track.num=2, side="in")

RCircos.Gene.Connector.Plot(genomic.data=gene.label, 
                            track.num=1, side="in")

RCircos.Gene.Name.Plot(gene.data=snps.label, name.col=4, 
                  track.num=4, side="in")

RCircos.Link.Plot(link.data=link.data, track.num=7, 
                  by.chromosome=FALSE)

RCircos.Line.Plot(line.data=cyto.info, data.col=4, 
                  	track.num=1, side="in")
dev.off()

pdf(file="RCircos_Interactions_semSNPs.pdf", height=8, width=8)

RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()

RCircos.Gene.Name.Plot(gene.data=gene.label, name.col=4, 
                       track.num=2, side="in")

RCircos.Gene.Connector.Plot(genomic.data=gene.label, 
                            track.num=1, side="in")


RCircos.Link.Plot(link.data=link.data, track.num=3, 
                  by.chromosome=FALSE)

dev.off()

data(RCircos.Histogram.Data)
data(RCircos.Gene.Label.Data)
data(RCircos.Ribbon.Data)
