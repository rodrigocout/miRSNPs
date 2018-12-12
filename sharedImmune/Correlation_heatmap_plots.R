###Correlation##########
library (RColorBrewer)    
source ('/Users/Rodrigo/Dropbox/BED_files_miRSNPS/TCGA/script/Sfunctions.R')

gxrna <- read.delim("~/Dropbox/shared_immune_mirsnps/data/Gene_names_mirsnps_geuvadis.txt")
gxrna <- gxrna[rowMeans(gxrna) > 1,]

gxmir <- read.delim("~/Dropbox/shared_immune_mirsnps/data/geuvadis_miRNAS_matched_rna.txt")
gxmir <- gxmir[rowMeans(gxmir) > 1,]


lista <- read.table("~/Dropbox/shared_immune_mirsnps/tables/results/exploratory/lista_microRNAs.txt")

miRInter <- intersect(rownames(gxmir), lista$V1)
length(miRInter)

targets <- read.table("~/Dropbox/shared_immune_mirsnps/tables/results/exploratory/Lista_target_genes_mirsnps.txt")
gensInter <- intersect(rownames(gxrna), targets$V1)

##----gexp intersections
ids <- colnames (gxmir)
gxrna <- gxrna[,ids]
gexp <- rbind (gxrna[gensInter,], gxmir[miRInter, ])


gexp1 <- log2(gexp + 1)

# By default calculates the distance between rows
dist1 = dist(t(gexp1))

colramp = colorRampPalette(c(3,"white",2))(9)
heatmap(as.matrix(dist1),col=colramp,Colv=NA,Rowv=NA)

##----correlation matrix
gexp <- t (gexp)
cormat <- cor (gexp[, miRInter], gexp[, gensInter], method = "spearman", use =  "complete.obs")
#colnames (cormat) <- names (gensInter)
#heatmap.cormat(cormat)
gexp1 <- t(gexp1)

cormat <- cor (gexp1[, miRInter], gexp1[, gensInter], method = "spearman", use =  "complete.obs")


colramp = colorRampPalette(c("red", "white", "blue"))
heatmap(cormat,col=colramp,Colv=NA,Rowv=NA)

library(gplots)
col3 <- colorRampPalette(c("red", "white", "blue"))

pdf("~/Dropbox/shared_immune_mirsnps/figs/correlation/Correlation_miRNAs_mRNAS.pdf")
heatmap.2(cormat,col=col3,Rowv=NA,Colv=NA,
          dendrogram="none", scale="row",trace="none",cexRow = 0.9, 
          cexCol = 0.9,
          margin=c(5,10), density.info=c("none"),  key.title = TRUE) 
dev.off()

cormat2 <- cor (gexp[, miRInter], gexp[, gensInter], method = "pearson", use =  "complete.obs")

pdf("~/Dropbox/shared_immune_mirsnps/figs/correlation/PearsonCorrelation_miRNAs_mRNAS.pdf")
heatmap.2(cormat2,col=col3,Rowv=NA,Colv=NA,
          dendrogram="none", scale="row",trace="none",cexRow = 0.9, 
          cexCol = 0.9,
          margin=c(5,10), main = "Pearson Correlation") 

dev.off()

mirs_list <- c("hsa-miR-326", "hsa-miR-3661", "hsa-miR-4463", "hsa-miR-4497", "hsa-miR-4518", 
               "hsa-miR-455-3p", "hsa-miR-4741", "hsa-miR-628-5p", "hsa-miR-660-5p")
length(mirs_list)

gene_list <- c("IKZF3", "SNAPC4", "FUT2", "IL18RAP", "UBE2L3","COG6", "MANBA" )

cormat2 <- cor (gexp[, mirs_list], gexp[, gene_list], method = "spearman", use =  "complete.obs")
               
#pdf("~/Dropbox/shared_immune_mirsnps/figs/correlation/Correlation_miRNAs_mRNAS_matched.pdf")

tiff("~/Dropbox/shared_immune_mirsnps/figs/correlation/Fig3_Correlation_miRNAs_mRNAS_matched.tiff", width = 593, height = 382, res = 400 )
graphics.off()
par("mar")
par(mar=c(1,1,1,1))

heatmap.2(cormat2,col=col3,Rowv=NA,Colv=NA,
          dendrogram="none", scale="row",trace="none",cexRow = 0.9, 
          cexCol = 0.9,
          margin=c(5,10),density.info=c("none"),  key.title = TRUE)
dev.off()

install.packages("ggcorrplot")
library(ggcorrplot)
               
ggcorrplot(cormat, outline.color = "black")

ggcorrplot(cormat2, hc.order = TRUE, 
            type = "lower",
           lab = TRUE)

# Get the lower triangle
ggcorrplot(cormat, hc.order = TRUE, type = "lower",
           outline.col = "white")

#########################
cormat <- cor (gexp1[, miRInter], gexp1[, gensInter], method = "pearson", use =  "complete.obs")

#install.packages("reshape2")
library(reshape2)
melted_cormat <- melt(cormat)
head(melted_cormat)

library(ggplot2)

pdf("~/Dropbox/shared_immune_mirsnps/figs/correlation/Pearson_Correlation_miRNAs_targets.pdf")

ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+   
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 9, hjust = 1))+
  coord_fixed()

dev.off()
#-----------------
cormat2 <- cor (gexp[, mirs_list], gexp[, gene_list], method = "pearson", use =  "complete.obs")

melted_cormat2 <- melt(cormat2)

library(ggplot2)

pdf("~/Dropbox/shared_immune_mirsnps/figs/correlation/Pearson_Correlation_miRNAs_targets_matched.pdf")

ggplot(data = melted_cormat2, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+   
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 9, hjust = 1))+
  coord_fixed()

dev.off()


# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(cormat)
upper_tri

melted_cormat <- melt(upper_tri, na.rm = TRUE)

#--------------------------
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
# Print the heatmap
print(ggheatmap)

##############
ggheatmap + 
  geom_text(aes(Var2, Var1, label = ""), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))


#####################
library(corrplot)
tiff("~/Dropbox/shared_immune_mirsnps/figs/correlationplot.tiff")
#col3 <- colorRampPalette(c("red", "white", "blue"))

tcormat <- t(cormat2)

corrplot(tcormat, method="color", tl.col = "black", tl.cex = 0.8)

corrplot(cormat, method="color", tl.col = "black", tl.cex = 0.8)

dev.off()

################################


