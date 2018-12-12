data<-read.table("India2.cd3.assoc.logistic", header=T)
data2<-read.table("NL.cd3.assoc.logistic", header=T)
#add=which(data$TEST=="ADD")
#data_r2<-read.table("IND.cd3.r0.8.assoc.logistic", header=T)
#data_trans<-read.table("IND.cd3.transfer.assoc.logistic", header=T)
index_bp<-181716045
index_p<-0.201
transfer_bp<-181849926
transfer_p<-0.001297
locus<- 'ITGA4/UBE2E3'
gene_info<-read.table("cd3genes.txt",header=T)
ncrna_info<-read.table("cd3ncrna.txt",header=T)
tiff("ITGA4-UBE2E3.regioassoc3.tiff",w=6,h=6,unit="in", res=400)
plot(data2$BP, -log10(data2$P), col=rgb(181,181,181,200,maxColorValue=255), pch=19, lwd=0.1, bty="n",xlim=range(181500000,182000000), ylim=range(0,5), xlab ="BP", ylab ="-LOG(p-val)",main=locus)
points(data$BP,-log10(data$P), col=rgb(0,0,139,150,maxColorValue=255), pch=19)
points(data_trans$BP,-log10(data_trans$P), col=rgb(255, 165, 0,255, maxColorValue=255), pch=19)
points(data_r2$BP,-log10(data_r2$P), col=rgb(0, 205, 0,110,maxColorValue=255), pch=19)
points(index_bp, -log10(index_p), col=rgb(238,0,0,238,maxColorValue=255), pch=19)#SNP-index
points(transfer_bp, -log10(transfer_p), col=rgb(0, 100, 0,255,maxColorValue=255), pch=19)
abline(h=-log10(0.05), lty="dashed")
abline(h=-log10(0.01), lty="dashed")
#abline(h=-log10(0.0004098), lty="dashed")
legend(181500000, 5, c("Top transfer-SNP", "r2>0.8 transfer-SNP", "All transferable SNPs", "Index-SNP"), cex=0.7, bty="n", pch=19, col= c(rgb(0, 100, 0,255,maxColorValue=255), rgb(0, 205, 0,110, maxColorValue=255) , rgb(255, 165, 0,255, maxColorValue=255), rgb(238,0,0,238,maxColorValue=255)))

for (i in 1:nrow(gene_info)){
  lines(range(gene_info[i,3],gene_info[i,5]), range(3.5,3.5), type="l", lwd=4,  col="turquoise 4", pch=19) #range values to specify the y-value of where to plot the lines for genes
  text(gene_info[i,4], range(3.60,3.60), gene_info[i,1], pos=4,col="black",cex=0.6)
}
for (i in 1:nrow(ncrna_info)){
  lines(range(ncrna_info[i,3],ncrna_info[i,5]), range(3.5,3.5), type="l", lwd=2,  col="indian red", pch=12) #range values to specify the y-value of where to plot the lines for ncrnas
  text(ncrna_info[i,4], range(3.60,3.60), ncrna_info[i,1], pos=4,col="black",cex=0.6)
} 
dev.off()

