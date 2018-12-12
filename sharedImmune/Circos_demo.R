This demo draw human chromosome ideogram and data tracks for:
  > #
  > #  1.	Connectors
  > #	2.	Gene lables
  > #	3.	Heatmap
  > #	4.	Scatterplot 
  > #	5.	Line plot
  > #	6.	Histogram
  > #	7.	Tile plot
  > #	8.	Link lines
  > #
  > # ________________________________________________________________________________________
  > # <><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><>
  > 
  > 
  > 
  > 	#	Load RCircos package and defined parameters
  > 	#  	_________________________________________________________________
  > 	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  > 	
  > 	library(RCircos);

> 	#	Load human cytoband data 
  > 	#  	_________________________________________________________________
  > 	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  > 
  > 	data(UCSC.HG19.Human.CytoBandIdeogram);

> 	hg19.cyto <- UCSC.HG19.Human.CytoBandIdeogram;

> 	#	Setup RCircos core components:
  > 	#
  > 	#	1. Chromosome ideogram plot information
  > 	#	2. x and y coordinates for a circular line and degrees of the
  > 	#		text rotation at each point
  > 	#	3. Plot parameters for image plot control
  > 	#  
  > 	#	These components will be stored in RCircos environment
  > 	#
  > 	#	Function arguments:
  > 	#
  > 	#	Chromosome ideogram data loaded above
  > 	#	Chromosomes need be excluded from cytoinfo
  > 	#	Number of tracks inside chromosome ideogram
  > 	#	Number of tracks inside chromosome ideogram
  > 	#  	_________________________________________________________________
  > 	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  > 
  > 	RCircos.Set.Core.Components(cyto.info=hg19.cyto, chr.exclude=NULL, 
                                 + 			tracks.inside=10, tracks.outside=0);

RCircos.Core.Components initialized.
Type ?RCircos.Reset.Plot.Parameters to see how to modify the core components.


> 	#	Open the graphic device (here a pdf file)
  > 	#
  > 	# 	png(file="RCircos.Demo.Human.png", height=8, width=8, unit="in", 
  > 	#		type="cairo", res=300);
  > 	#
  >  	#	tiff(file="RCircos.Demo.Human.tif", height=8, width=8, unit="in", 
  > 	#		type="cairo", res=300);
  > 	#  	_________________________________________________________________
  > 	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  > 
  > 	cat("Open graphic device and start plot ...\n");
Open graphic device and start plot ...

> 	pdf(file="RCircos.Demo.Human.pdf", height=8, width=8);

> 	RCircos.Set.Plot.Area();
Hit <Return> to see next plot: 
  
  > 	title("RCircos 2D Track Plot with Human Genome");

> 	#	Draw chromosome ideogram
  > 	#  	_________________________________________________________________
  > 	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  > 
  > 	cat("Draw chromosome ideogram ...\n");
Draw chromosome ideogram ...

> 	RCircos.Chromosome.Ideogram.Plot();

> 	#	Connectors in first track and gene names in the second track. 
  > 	#  	_________________________________________________________________
  > 	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  > 
  > 	cat("Add Gene and connector tracks ...\n");
Add Gene and connector tracks ...

> 	data(RCircos.Gene.Label.Data);

> 	RCircos.Gene.Connector.Plot(genomic.data=RCircos.Gene.Label.Data, 
                               + 			track.num=1, side="in");

Maximum lables for  chr1 is 16 . Extra ones are ignored.
Maximum lables for  chr3 is 13 . Extra ones are ignored.
Maximum lables for  chr7 is 10 . Extra ones are ignored.
Maximum lables for  chr9 is 9 . Extra ones are ignored.
Maximum lables for  chr10 is 9 . Extra ones are ignored.
Maximum lables for  chr11 is 9 . Extra ones are ignored.
Maximum lables for  chr12 is 8 . Extra ones are ignored.
Maximum lables for  chr16 is 6 . Extra ones are ignored.
Maximum lables for  chr17 is 5 . Extra ones are ignored.
Maximum lables for  chr18 is 5 . Extra ones are ignored.
Maximum lables for  chr19 is 3 . Extra ones are ignored.
Maximum lables for  chr20 is 4 . Extra ones are ignored.


> 	RCircos.Gene.Name.Plot(gene.data=RCircos.Gene.Label.Data, name.col=4, 
                          + 			track.num=2, side="in");

Maximum lables for  chr1 is 16 . Extra ones are ignored.
Maximum lables for  chr3 is 13 . Extra ones are ignored.
Maximum lables for  chr7 is 10 . Extra ones are ignored.
Maximum lables for  chr9 is 9 . Extra ones are ignored.
Maximum lables for  chr10 is 9 . Extra ones are ignored.
Maximum lables for  chr11 is 9 . Extra ones are ignored.
Maximum lables for  chr12 is 8 . Extra ones are ignored.
Maximum lables for  chr16 is 6 . Extra ones are ignored.
Maximum lables for  chr17 is 5 . Extra ones are ignored.
Maximum lables for  chr18 is 5 . Extra ones are ignored.
Maximum lables for  chr19 is 3 . Extra ones are ignored.
Maximum lables for  chr20 is 4 . Extra ones are ignored.


> 	#	Heatmap plot.  Since some gene names plotted above are longer 
  > 	#	than one track height, we skip two tracks 
  > 	#  	_________________________________________________________________
  > 	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  > 
  > 	cat("Add heatmap track ...\n");
Add heatmap track ...

> 	data(RCircos.Heatmap.Data);

> 	RCircos.Heatmap.Plot(heatmap.data=RCircos.Heatmap.Data, data.col=6, 
                        + 			track.num=5, side="in");

> 	#	Scatterplot
  > 	#  	_________________________________________________________________
  > 	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  > 
  > 	cat("Add scatterplot track ...\n");
Add scatterplot track ...

> 	data(RCircos.Scatter.Data)

> 	RCircos.Scatter.Plot(scatter.data=RCircos.Scatter.Data, data.col=5, 
                        + 			track.num=6, side="in", by.fold=1);

> 	#	Line plot.
  > 	#  	_________________________________________________________________
  > 	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  > 
  > 	cat("Add line plot track ...\n");
Add line plot track ...

> 	data(RCircos.Line.Data);

> 	RCircos.Line.Plot(line.data=RCircos.Line.Data, data.col=5, 
                     + 			track.num=7, side="in");

> 	#	Histogram plot
  > 	#  	_________________________________________________________________
  > 	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  > 
  > 	cat("Add histogram track ...\n");
Add histogram track ...

> 	data(RCircos.Histogram.Data);

> 	RCircos.Histogram.Plot(hist.data=RCircos.Histogram.Data, data.col=4, 
                          + 			track.num=8, side="in");

> 	#	Tile plot. Note: tile plot data have chromosome locations and each
  > 	#	data file is for one track
  > 	#  	_________________________________________________________________
  > 	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  > 
  > 	cat("Add tile track ...\n");
Add tile track ...

> 	data(RCircos.Tile.Data);

> 	RCircos.Tile.Plot(tile.data=RCircos.Tile.Data, track.num=9, side="in");

> 	#	Link lines. Link data has only paired chromosome locations in
  > 	#	each row and link lines are always drawn inside of chromosome 
  > 	#	ideogram.
  > 	#  	_________________________________________________________________
  > 	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  > 
  > 	cat("Add link track ...\n");
Add link track ...

> 	data(RCircos.Link.Data);

> 	RCircos.Link.Plot(link.data=RCircos.Link.Data, track.num=11, 
                     + 				by.chromosome=FALSE);

> 	#	Add ribbon link to the center of plot area (link lines).
  > 	#  	_________________________________________________________________
  > 	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  > 
  > 	cat("Add ribbons to link track ...\n");
Add ribbons to link track ...

> 	data(RCircos.Ribbon.Data);

> 	RCircos.Ribbon.Plot(ribbon.data=RCircos.Ribbon.Data, track.num=11, 
                       + 		by.chromosome=FALSE, twist=FALSE);

> 	#	Close the graphic device and clear memory
  > 	#  	_________________________________________________________________
  > 	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  > 
  > 	dev.off();
RStudioGD 
2 

> 	cat("R Circos Demo Done ...\n\n");
R Circos Demo Done ...
