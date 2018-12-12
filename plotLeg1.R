##plot legends
plotleg1<-function(leg,title="", fname="leg",plot=FALSE, width=5,height=2.5){
  if(!plot)pdf(file=paste(fname,".pdf",sep=""), width=width, height=height)
  par(mar = c(5,1,2.5,1))
  image(matrix(1:length(leg)),col=leg,axes=FALSE)
  if(!is.null(names(leg))){
    xlab<-names(leg)
    axis(1, at=seq(0, 1, length.out = length(xlab)), tcl = -0.2,
         labels=xlab, las=3,cex.axis=1, mgp=c(1,0.4,0),lwd=1)
  }
  box()
  mtext(title,adj=0,line=1,cex=1.2)
  if(!plot)dev.off()
  if(!plot)cat(paste("Legend '", fname,".pdf' generated!\n\n", sep=""))
}
