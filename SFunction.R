source("~/Dropbox/BED_files_miRSNPS/TCGA/script/plotleg1.R")
##############################################
####Functions of miRSNPs project
##############################################
naFilter <- function (gexp, maxPercentNas = 15, byNa = TRUE)
{
  ##----checks
  if (class (gexp) != "matrix")
    stop ("'gexp' should be an expression matrix!")
  if (class (maxPercentNas) != "numeric")
    stop ("'maxPercentNas' should be numeric!")
  
  ##----percent avaliation
  if (byNa)
  {
    ids <- apply (gexp, 1, function (nasum)
      {
        val <- (sum (is.na (nasum)) / length (nasum)) * 100
    })
  } else
  {
    ids <- apply (gexp, 1, function (nasum)
    {
      val <- (sum (nasum == 0) / length (nasum)) * 100
    })
  }
  idx <- which (ids <= maxPercentNas)
  gexp <- gexp[idx, ]
  return (gexp)
}

##############################################
####Nas barplot
##############################################
barPlotNas <- function (gexp, miRNavy = NULL)
{
  ##----checks
  if (class (gexp) != "matrix")
    stop ("'gexp' should be an expression matrix!")
  ##----percent avaliation
  percent <- list ()
  for (i in seq (0, 100, by = 10))
  {
    ids <- apply (gexp, 1, function (nasum)
    {
      val <- (sum (is.na (nasum)) / length (nasum)) * 100
    })
    #ids[ids == 0] <- 1e-5
    idx <- which (ids <= i)
    percent[[paste (i, "%", sep = "")]] <- ids[idx]
  }
  if (!is.null (miRNavy))
  {
    intersect <- list ()
    for (i in names (percent))
    {
      nms <- names (percent[[i]])
      inter <- intersect (nms, miRNavy)
      intersect[[i]] <- inter
    }
    bl<-brewer.pal(9,'Greens')[5]
    #barplot(rep(1,length (bl)), yaxt="n", col=bl)
    ylim <- c(0, 1.1 * max(unlist(lapply(intersect, length))))
    
    pdf(file=paste('intersect','.pdf',sep=''), width = 8, height = 5)
    
    par(mgp = c(2.5,0.5,0))
    xx <- barplot(unlist(lapply(intersect, length)),xlab = "Limite m??ximo de NAs (%)", 
                  ylab = "N??mero de miRNAs", axes = F,
                  axisnames = F, ylim = ylim, col = bl[1], border = NA)
    text(x = xx, y = unlist(lapply(intersect, length)) , label = unlist(lapply(intersect, length)),
         pos = 3, cex = 0.8, col = "darkslategrey", font = 2)
    axis(1, at=xx, labels=names(intersect), tick=T, las=1, line=0.5, cex.axis=1)
    axis(2 , las = 1)
    
    dev.off()
  } else
  {
    bl<-brewer.pal(9,'Blues')[4]
    #barplot(rep(1,length (bl)), yaxt="n", col=bl)
    ylim <- c(0, 1.1 * max(unlist(lapply(percent, length))))
    
    pdf(file=paste('percent','.pdf',sep=''), width = 8, height = 5)
    
    par(mgp = c(2.5,0.5,0))
    xx <- barplot(unlist(lapply(percent, length)),xlab = "Limite m??ximo de NAs (%)", 
                  ylab = "N??mero de miRNAs", axes = F,
                  axisnames = F, ylim = ylim, col = bl[1], border = NA)
    text(x = xx, y = unlist(lapply(percent, length)) , label = unlist(lapply(percent, length)),
         pos = 3, cex = 0.8, col = "darkblue", font = 2)
    axis(1, at=xx, labels=names(percent), tick=T, las=1, line=0.5, cex.axis=1)
    axis(2 , las = 1)
    
    dev.off()
  }
  
}


# ##############################################
# ####Order the gexp (more expression per miR and Sample)
# ##############################################
# expOrder <- function (gexp, decreasing = TRUE, bySample = TRUE)
# {
#   ##----checks
#   if (class (gexp) != "matrix")
#     stop ("'gexp' should be an expression matrix!")
#   if (is.null(colnames (gexp)) || is.null(rownames (gexp)))
#     stop ("'gexp' should have rownames and colnames !")
#   ##----get indices
#   idsR <- apply (gexp, 1, function (mRow)
#     {
#       meanRows <- mean (mRow, na.rm = TRUE)
#       ##coeficiente de variacao (sd/meadian)
#   })
#   idxR <- names (sort (idsR, decreasing = decreasing))
#   gexp <- gexp[idxR, ]
#   if (bySample)
#   {
#     idsC <- apply (gexp, 2, function (mCol)
#     {
#       meanCols <- mean (mCol, na.rm = TRUE)
#     })
#     idxC <- names (sort (idsC, decreasing = decreasing))
#     gexp <- gexp[, idxC]
#   }
#   return (gexp)
# }

##############################################
####characteristics extration (sd, std.error, var)
##############################################
gxmirChar <- function (gexp)
{
  sd <- apply (gexp, 1, function (x)
  {
    var <- sd (x, na.rm = T)
  })
  #coeficiente de var: sd/median
  std.error <- apply (gexp, 1, function (x)
    {
      var <- sd (x, na.rm = T) / sqrt (length(x))
  })
  
  variance <- apply (gexp, 1, function (x)
    {
      variance <- var (x, na.rm = T)
  })
  
  return (cbind (sd, std.error, variance))
}

# ##############################################
# ####Heatmap gexp
# ##############################################
# source ("~/Dropbox/reportGordon/tablesAndPlots/plotleg1.R")
# heatmap.gexp <- function(gexp)
# {
#   gexp[is.na(gexp)] <- 0
#   gexpz <- t (apply (gexp, 1, function (x)
#   {
#     z <- (x - median (x)) / sd(x)
#   }))
#   
#   #---
#   bl<-brewer.pal(7,'Blues')[7:4]
#   #barplot(rep(1,length (bl)), yaxt="n", col=bl)
#   rd <- brewer.pal(7,'Reds')[3:7]
#   #barplot(rep(1,length (rd)), yaxt="n", col=rd)
#   #hgexp<-t(rbind(regn,regp))
#   gexpz[gexpz<(-7)] <- -7
#   gexpz[gexpz>7] <- 7
#   #---
#   pdf(file=paste('Gexp_Heat','.pdf',sep=''), width = 5, height = 5)
#   
#   heatmap(gexpz, Rowv = NA, Colv = NA, breaks = quantile (gexpz, probs = seq (0,1,0.1)), 
#           col=c(bl,"white",rd), labRow=NA, labCol=NA, margins = c(1,1))
#   
#   dev.off()
#   leg <- c(bl[c(1,3,4)],rd[c(1,3,4)]); names(leg) <- names(quantile (gexpz, probs = seq (0,1,0.1))[c(1,2,4,6,9,11)])
#   plotleg1 (leg = leg)
#   
# }

##############################################
####Heatmap cormat
##############################################
heatmap.cormat <- function(cormat)
{

  #---
  bl<-brewer.pal(9,'Blues')[9:6]
  #barplot(rep(1,length (bl)), yaxt="n", col=bl)
  rd <- brewer.pal(9,'Reds')[6:9]
  #barplot(rep(1,length (c(bl,rd))), yaxt="n", col=c(bl,rd))
  #hgexp<-t(rbind(regn,regp))
  cormat[cormat < -0.4] <- -0.4
  cormat[cormat > 0.4] <- 0.4
  #---
  pdf(file=paste('Cormat_Heat','.pdf',sep=''), width = 6, height = 5)

  heatmap(cormat, scale = "none", breaks = seq (-0.4,0.4,0.1),
          col=c(bl,rd))

  dev.off()
  leg <- c(bl,rd); names(leg) <-seq (-0.4,0.4,0.1)[-5]
  plotleg1 (leg = leg)

}
