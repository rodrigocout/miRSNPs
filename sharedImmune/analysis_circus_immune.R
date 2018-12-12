######Circus plot for shared mirsnps between autoimmune diseases##########
setwd('~/Dropbox/shared_immune_mirsnps/scripts/circos_data/')
library(RCircos)

#data(UCSC.HG19.Human.CytoBandIdeogram)
#head(UCSC.HG19.Human.CytoBandIdeogram)
data(RCircos.Line.Data)
head(RCircos.Line.Data)
write.table(RCircos.Line.Data, 'line.data.txt', sep = '\t', quote = F, row.names = F)
#line.data
line.data <- read.delim('line.data.txt')
head(line.data)

##make a link data
link.data <- read.delim('cicus_link_test.txt')
head(link.data, n = 10)

#shared genes names and diseases
gene.label <- read.delim('gene_label_mirsnps.txt')
head(gene.label)

#mirnas names
mirnas.label <- read.delim('mirnas_circus.txt')
head(mirnas.label)

#test
cyto.info <- read.delim('chr_circos_hg19.txt')
head(cyto.info)
chr.exclude <- c( "chrY","chr12", "chr13", "chr14",'chr15', "chr16", "chr17","chr19", "chr18", "chr20", "chr21", "chr22")
tracks.inside <- 5
tracks.outside <- 0

RCircos.Set.Core.Components(cyto.info, chr.exclude,
                            tracks.inside, tracks.outside)

###make a plot
pdf(file="RCircos.pdf", height=8, width=8)

RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

RCircos.Line.Plot(line.data=gene.label, data.col=5, 
                  track.num=4, side="in")


RCircos.Gene.Name.Plot(gene.data=gene.label, name.col=4, 
                       track.num=1, side="in")

RCircos.Link.Plot(link.data=link.data, track.num=3, 
                  by.chromosome=FALSE)

dev.off()

#RCircos.Gene.Connector.Plot(genomic.data=gene.label, 
#                            track.num=1, side="in")


RCircos.Link.Plot(link.data=link.data, track.num=3, 
                  by.chromosome=FALSE)

RCircos.Line.Plot(line.data=snps.label, data.col=5, 
                  track.num=6, side="in")


dev.off()



##############################################################
rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.cyto <- RCircos.Get.Plot.Ideogram()
rcircos.position <- RCircos.Get.Plot.Positions()
RCircos.List.Parameters()

##Modifying
rcircos.params <- RCircos.Get.Plot.Parameters();
rcircos.params$base.per.unit <- 30000;
RCircos.Reset.Plot.Parameters(rcircos.params);
RCircos.List.Parameters();

##Making a plot
out.file <- "RCircosDemoHumanGenome.pdf"
pdf(file=out.file, height=8, width=8, compress=TRUE)

RCircos.Set.Plot.Area()

par(mai=c(0.25, 0.25, 0.25, 0.25))
plot.new()
plot.window(c(-2.5,2.5), c(-2.5, 2.5))

RCircos.Chromosome.Ideogram.Plot()

#data(RCircos.Gene.Label.Data)
gene.label <- read.delim('gene.label.txt')
name.col <- 4;
side <- "in";
track.num <- 1;

RCircos.Gene.Connector.Plot(gene.label,
                            track.num, side)
track.num <- 2;

RCircos.Gene.Name.Plot(gene.label,
                       name.col,track.num, side)


#The code below will draw data tracks of heatmap, scatter, line, histogram,
#and tile plots.

data(RCircos.Heatmap.Data)
data.col <- 6

track.num <- 5
side <- "in"
RCircos.Heatmap.Plot(RCircos.Heatmap.Data, data.col,
                     track.num, side)

data(RCircos.Scatter.Data)
data.col <- 5
track.num <- 6
side <- "in"
by.fold <- 1
RCircos.Scatter.Plot(RCircos.Scatter.Data, data.col,
                     track.num, side, by.fold)
data(RCircos.Line.Data)
data.col <- 5
track.num <- 7
side <- "in"
RCircos.Line.Plot(RCircos.Line.Data, data.col,
                  track.num, side)

data(RCircos.Tile.Data)
track.num <- 9
side <- "in";
RCircos.Tile.Plot(RCircos.Tile.Data, track.num, side)

dev.off()




