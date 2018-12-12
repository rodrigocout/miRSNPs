###Analysis mirSNPs association per disease ###############
setwd("/Users/Rodrigo/Documents/Pos_doc/IdeaS/Immunobase/proceed_data/")

library(circlize)

mirs <- read.table("circlize_disease_overlapp_mirnas.txt", skip = 1, stringsAsFactors = FALSE)

chordDiagram(mirs)
circos.info()
circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), mean(ylim), sector.index, col = "white", cex = 0.6, facing = "inside", niceFacing = TRUE)
}, bg.border = NA)

#########Alternative#########################
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  if(abs(xplot[2] - xplot[1]) > 20) {
    circos.lines(xlim, c(mean(ylim), mean(ylim)), lty = 3) # dotted line
    for(p in seq(0.2, 1, by = 0.2)) {
      circos.text(p*(xlim[2] - xlim[1]) + xlim[1], mean(ylim) + 0.1,
                  p, cex = 0.3, adj = c(0.5, 0), niceFacing = TRUE)
    }
  }
  circos.text(mean(xlim), 1, sector.name, niceFacing = TRUE, adj = c(0.5, 0))
}, bg.border = NA)
######################################################################

chordDiagram(mirs, 
             preAllocateTracks = list(track.height = 0.3))

circos.text(3, 9, "Autoimmune Diseases", facing = "inside", cex = 0.8)


chordDiagram(mirs, annotationTrack = c("name", "grid"),
             annotationTrackHeight = c(0.09, 0.01))

circos.clear()

