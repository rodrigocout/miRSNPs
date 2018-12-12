#####MH Plot#####

#####################
##Crohn's disease####
load(file = "/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/tables/cd_mirsnps.RData")
library(manhattanly)

head(cd_mirsnps)
cd_mirsnps <- subset(cd_mirsnps, select = c("Chr", "Position", "PValue", "Marker"))
head(cd_mirsnps)
colnames(cd_mirsnps) <- c("CHR", "BP", "P", "SNP")
cd < cd [as.numeric(cd$CHR),]
write.table(cd_mirsnps,"/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/tables/cd_mirsnps_all.txt", sep = '\t', row.names = F, quote = F)

cd_mirsnps <- read.table("/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/tables/cd_mirsnps_all.txt", header = T)


manhattanly(subset(cd_mirsnps, CHR %in% 1:22), snp = "SNP")

tiff("/Users/lgmh/Dropbox/shared_immune_mirsnps/figs/MH_plot_mirsnps_Crohn.tif", width = 890, height = 675, pointsize = 19)
manhattanly(cd_mirsnps, suggestiveline = -log10(5e-06),  genomewideline = -log10(5e-08),  title = "MirSNPs Crohn's disease")
dev.off()

###Alternative#####
library(qqman)

tiff("/Users/lgmh/Dropbox/shared_immune_mirsnps/figs/manhatan/MH_plot_mirsnps_Crohn.tif", width = 890, height = 675, pointsize = 19)
manhattan(cd_mirsnps,main = "miRSNPs Crohn's disease",cex = 0.9,cex.axis = 0.9,ylim = c(0, 45),
          col = c("blue4", "deepskyblue"), suggestiveline = -log10(5e-06), genomewideline = -log10(5e-08))
dev.off()

###########Non HLA########################################
#GRCh37
#chr6:28477797-33448354 

hla_cd <- cd_mirsnps[ cd_mirsnps$CHR == "6",]
hla_cd <- hla_cd[hla_cd$BP > 25780582 & hla_cd$BP < 33448354,] 

nonhla_cd <- cd_mirsnps[!(cd_mirsnps$SNP %in% hla_cd$SNP),]

tiff("/Users/lgmh/Dropbox/shared_immune_mirsnps/figs/manhatan/mh_plot_nonhla_mirnas_cd.tif",  width = 890, height = 675, pointsize = 19)
manhattan(nonhla_cd,main = "Non-HLA miRSNPs Crohn's disease",cex = 0.9,cex.axis = 0.9, ylim = c(0, 45),
          col = c("blue4", "deepskyblue"), suggestiveline = -log10(5e-06), genomewideline = -log10(5e-08))
dev.off()

####Write table##############
mirsnps_cd <- cd_mirsnps[cd_mirsnps$P < 5e-06,]
write.table(mirsnps_cd, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/sign_mirsnps_cd.txt", sep = '\t', 
            quote = F, row.names = F)


nonhla_mirsnps_cd <- nonhla_cd[nonhla_cd$P < 5e-06,]
write.table(nonhla_mirsnps_cd, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/sign_nonhla_mirsnps_cd.txt", sep = '\t', 
            quote = F, row.names = F)

##########
##ATD#####
load(file = "/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/tables/atd_mirsnps.Rdata")

atd_mirsnp <- read.table("/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/tables/atd_mirsnps.txt", header=T, stringsAsFactors=F)
head(atd_mirsnp)

atd_mirsnp <- subset(atd_mirsnp, select = c("Chr", "Position", "PValue", "Marker"))
head(atd_mirsnp)
colnames(atd_mirsnp) <- c("CHR", "BP", "P", "SNP")

#Get rid off CHR X and Y
atd_mirsnp <- atd_mirsnp[!atd_mirsnp$CHR == "X",]
atd_mirsnp <- atd_mirsnp[!atd_mirsnp$CHR == "Y",]
#Check dataframe
lapply(atd_mirsnp, data.class)

#Convert column to numeric
atd_mirsnp[, c(1)] <- sapply(atd_mirsnp[, c(1)], as.numeric)

tiff("/Users/lgmh/Dropbox/shared_immune_mirsnps/figs/manhatan/mh_plot_mirnas_atd.tif",  width = 890, height = 675, pointsize = 19)
manhattan(atd_mirsnp,main = "miRSNPs ATD",cex = 0.9,cex.axis = 0.9,ylim =c(0,50),
          col = c("blue4", "deepskyblue"), suggestiveline = -log10(5e-06), genomewideline = -log10(5e-08))
dev.off()
###########With no HLA########################################
#GRCh37
#chr6:28477797-33448354 

hla_atd <- atd_mirsnp[ atd_mirsnp$CHR == "6",]
hla_atd <- hla_atd[hla_atd$BP > 28477797,] 
hla_atd <- hla_atd[hla_atd$BP < 33448354,]

nonhla_atd <- atd_mirsnp[!(atd_mirsnp$SNP %in% hla_atd$SNP),]

tiff("/Users/lgmh/Dropbox/shared_immune_mirsnps/figs/manhatan/mh_plot_nohla_mirnas_atd.tif",  width = 890, height = 675, pointsize = 19)
manhattan(nonhla_atd,main = "Non-HLA miRSNPs ATD",cex = 0.9,cex.axis = 0.9, ylim = c(0, 25),
          col = c("blue4", "deepskyblue"), suggestiveline = -log10(5e-06), genomewideline = -log10(5e-08))
dev.off()

####Write table##############
mirsnps_atd <- atd_mirsnp[atd_mirsnp$P < 5e-06,]

write.table(mirsnps_atd, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/sign_mirsnps_atd.txt", sep = '\t', 
            quote = F, row.names = F)

nonhla_mirsnps_atd <- nonhla_atd[nonhla_atd$P < 5e-06,]

write.table(nonhla_mirsnps_atd, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/sign_nonhla_mirsnps_atd.txt", sep = '\t', 
            quote = F, row.names = F)

###########################################################################################################
library(manhattanly)
manhattanly(atd_mirsnp, suggestiveline = -log10(5e-06),  genomewideline = -log10(5e-08),  title = "mirSNPs ATD")

#################################################################################
##Celiac#######

cel_mirs <- read.table("/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/tables/cel_mirsnps.txt", header = T, stringsAsFactors = FALSE)
head(cel_mirs)

cel_mirs <- subset(cel_mirs, select = c("Chr", "Position", "PValue", "Marker"))
head(cel_mirs)
colnames(cel_mirs) <- c("CHR", "BP", "P", "SNP")

#Get rid off CHR X and Y
cel_mirs <- cel_mirs[!cel_mirs$CHR == "X",]
cel_mirs <- cel_mirs[!cel_mirs$CHR == "Y",]
#Check dataft1dme
lapply(cel_mirs, data.class)

#Convert column to numeric
cel_mirs$CHR = as.numeric(as.character(cel_mirs$CHR))

tif("/Users/lgmh/Dropbox/shared_immune_mirsnps/figs/manhatan/mh_plot_mirnas_cel.tif",  width = 890, height = 675, pointsize = 19)
manhattan(cel_mirs,main = "miRSNPs Cel",cex = 0.9,cex.axis = 0.9,ylim =c(0,330),
          col = c("blue4", "deepskyblue"), suggestiveline = -log10(5e-06), genomewideline = -log10(5e-08))
dev.off()

###########Non HLA########################################
#GRCh37
#chr6:28477797-33448354 

hla_cel <- cel_mirs[ cel_mirs$CHR == "6",]
hla_cel <- hla_cel[hla_cel$BP > 25780582,] 
hla_cel <- hla_cel[hla_cel$BP < 33448354,]

nonhla_cel <- cel_mirs[!(cel_mirs$SNP %in% hla_cel$SNP),]

tiff("/Users/lgmh/Dropbox/shared_immune_mirsnps/figs/manhatan/mh_plot_nohla_mirnas_cel.tif",  width = 890, height = 675, pointsize = 19)
manhattan(nonhla_cel,main = "Non-HLA miRSNPs CEL",cex = 0.9,cex.axis = 0.9, ylim = c(0, 20), 
          col = c("blue4", "deepskyblue"), suggestiveline = -log10(5e-06), genomewideline = -log10(5e-08))
dev.off()


####Write table ##########################
mirsnps_cel <- cel_mirs[cel_mirs$P < 5e-06,]

write.table(mirsnps_cel, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/sign_mirsnps_cel.txt", sep = '\t', 
            quote = F, row.names = F)

nonhla_mirsnps_cel <- nonhla_cel[nonhla_cel$P < 5e-06,]

write.table(nonhla_mirsnps_cel, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/sign_nonhla_mirsnps_cel.txt", sep = '\t', 
            quote = F, row.names = F)

####################################################################################
##JIA############

jia_mirs <- read.table("/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/tables/jia_mirsnps.txt", header = T)

jia_mirs <- subset(jia_mirs, select = c("Chr", "Position", "PValue", "Marker"))
head(jia_mirs)
colnames(jia_mirs) <- c("CHR", "BP", "P", "SNP")

#Get rid off CHR X and Y
jia_mirs <- jia_mirs[!jia_mirs$CHR == "X",]
jia_mirs <- jia_mirs[!jia_mirs$CHR == "Y",]

#Check dataft1dme
lapply(jia_mirs, data.class)


#Convert column to numeric
jia_mirs$CHR = as.numeric(as.character(jia_mirs$CHR))

tiff("/Users/lgmh/Dropbox/shared_immune_mirsnps/figs/manhatan/mh_plot_mirnas_jia.tif", width = 890, height = 675, pointsize = 19)
manhattan(jia_mirs,main = "miRSNPs JIA",cex = 0.9,cex.axis = 0.9,
          col = c("blue4", "deepskyblue"), suggestiveline = -log10(5e-06), genomewideline = -log10(5e-08))
dev.off()

###########Non HLA########################################
#GRCh37
#chr6:28477797-33448354 

hla_jia <- jia_mirs[ jia_mirs$CHR == "6",]
hla_jia <- hla_jia[hla_jia$BP > 25780580 & hla_jia$BP < 33448354,] 

nonhla_jia <- jia_mirs[!(jia_mirs$SNP %in% hla_jia$SNP),]

tiff("/Users/lgmh/Dropbox/shared_immune_mirsnps/figs/manhatan/mh_plot_nonhla_mirnas_jia.tif", width = 890, height = 675, pointsize = 19)
manhattan(nonhla_jia,main = "Non-HLA miRSNPs JIA",cex = 0.9,cex.axis = 0.9,
          col = c("blue4", "deepskyblue"), suggestiveline = -log10(5e-06), genomewideline = -log10(5e-08))
dev.off()


####Write table ##########################
mirsnps_jia <- jia_mirs[jia_mirs$P < 5e-06,]

write.table(mirsnps_jia, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/sign_mirsnps_jia.txt", sep = '\t', 
            quote = F, row.names = F)

nonhla_mirsnps_jia <- nonhla_jia[nonhla_jia$P < 5e-06,]

write.table(nonhla_mirsnps_jia, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/sign_nonhla_mirsnps_jia.txt", sep = '\t', 
            quote = F, row.names = F)


#############################################################################
####Narcolepsy####################################################

nar_mirs <- read.table("/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/tables/nar_mirsnps.txt", header = T)
head(nar_mirs)

nar_mirs <- subset(nar_mirs, select = c("Chr", "Position", "PValue", "Marker"))
head(nar_mirs)
colnames(nar_mirs) <- c("CHR", "BP", "P", "SNP")

#Get rid off CHR X and Y
nar_mirs <- nar_mirs[!nar_mirs$CHR == "X",]
nar_mirs <- nar_mirs[!nar_mirs$CHR == "Y",]

#Check dataft1dme
lapply(nar_mirs, data.class)

#Convert column to numeric
nar_mirs$CHR = as.numeric(as.character(nar_mirs$CHR))

tiff("/Users/lgmh/Dropbox/shared_immune_mirsnps/figs/manhatan/mh_plot_mirnas_nar.tif", width = 890, height = 675, pointsize = 19)
manhattan(nar_mirs,main = "miRSNPs Narcolepsy",cex = 0.9,cex.axis = 0.8,
          col = c("blue4", "deepskyblue"), suggestiveline = -log10(5e-06), genomewideline = -log10(5e-08))
dev.off()


##############################
####PBC#####

pbc_mirs <- read.table("/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/tables/pbc_mirsnps.txt", header = T)
head(pbc_mirs)

pbc_mirs <- subset(pbc_mirs, select = c("Chr", "Position", "PValue", "Marker"))
head(pbc_mirs)
colnames(pbc_mirs) <- c("CHR", "BP", "P", "SNP")

#Get rid off CHR X and Y
pbc_mirs <- pbc_mirs[!pbc_mirs$CHR == "X",]
pbc_mirs <- pbc_mirs[!pbc_mirs$CHR == "Y",]

#Check dataft1dme
lapply(pbc_mirs, data.class)

#Convert column to numeric
pbc_mirs$CHR = as.numeric(as.chat1dcter(pbc_mirs$CHR))

tiff("/Users/lgmh/Dropbox/shared_immune_mirsnps/figs/manhatan/mh_plot_mirnas_pbc.tif", width = 890, height = 675, pointsize = 19)
manhattan(pbc_mirs,main = "miRSNPs PBC",cex = 0.9,cex.axis = 0.8,ylim = c(0,50),
          col = c("blue4", "deepskyblue"), suggestiveline = -log10(5e-06), genomewideline = -log10(5e-08))
dev.off()
###########Non HLA########################################
#GRCh37
#chr6:28477797-33448354 

hla_pbc <- pbc_mirs[ pbc_mirs$CHR == "6",]
hla_pbc <- hla_pbc[hla_pbc$BP > 25780582 & hla_pbc$BP < 33448354,] 

nonhla_pbc <- pbc_mirs[!(pbc_mirs$SNP %in% hla_pbc$SNP),]

tiff("/Users/lgmh/Dropbox/shared_immune_mirsnps/figs/manhatan/mh_plot_nohla_mirnas_pbc.tif",  width = 890, height = 675, pointsize = 19)
manhattan(nonhla_pbc,main = "Non-HLA miRSNPs PBC",cex = 0.9,cex.axis = 0.9, ylim = c(0, 25), 
          col = c("blue4", "deepskyblue"), suggestiveline = -log10(5e-06), genomewideline = -log10(5e-08))
dev.off()


####Write table ##########################
mirsnps_pbc <- pbc_mirs[pbc_mirs$P < 5e-06,]

write.table(mirsnps_pbc, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/sign_mirsnps_pbc.txt", sep = '\t', 
            quote = F, row.names = F)

nonhla_mirsnps_pbc <- nonhla_pbc[nonhla_pbc$P < 5e-06,]

write.table(nonhla_mirsnps_pbc, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/sign_nonhla_mirsnps_pbc.txt", sep = '\t', 
            quote = F, row.names = F)

#############################################
####PSO#################

pso_mirs <- read.table("/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/tables/pso_mirsnps.txt", header = T)
head(pso_mirs)

pso_mirs <- subset(pso_mirs, select = c("Chr", "Position", "PValue", "Marker"))
head(pso_mirs)
colnames(pso_mirs) <- c("CHR", "BP", "P", "SNP")

#Get rid off CHR X and Y
pso_mirs <- pso_mirs[!pso_mirs$CHR == "X",]
pso_mirs <- pso_mirs[!pso_mirs$CHR == "Y",]

#Check dataft1dme
lapply(pso_mirs, data.class)

#Convert column to numeric
pso_mirs$CHR = as.numeric(as.chat1dcter(pso_mirs$CHR))

tiff("/Users/lgmh/Dropbox/shared_immune_mirsnps/figs/manhatan/mh_plot_mirnas_pso.tif",width = 890, height = 675, pointsize = 19 )
manhattan(pso_mirs,main = "miRSNPs PSO",cex = 0.9,cex.axis = 0.8,
          col = c("blue4", "deepskyblue"), suggestiveline = -log10(5e-06), genomewideline = -log10(5e-08))
dev.off()

###########Non HLA########################################
#GRCh37
#chr6:28477797-33448354 

hla_pso <- pso_mirs[ pso_mirs$CHR == "6",]
hla_pso <- hla_pso[hla_pso$BP > 25780582 & hla_pso$BP < 33448354,] 

nonhla_pso <- pso_mirs[!(pso_mirs$SNP %in% hla_pso$SNP),]

tiff("/Users/lgmh/Dropbox/shared_immune_mirsnps/figs/manhatan/mh_plot_nohla_mirnas_pso.tif",  width = 890, height = 675, pointsize = 19)
manhattan(nonhla_pso,main = "Non-HLA miRSNPs PSO",cex = 0.9,cex.axis = 0.9, ylim = c(0, 20),
          col = c("blue4", "deepskyblue"), suggestiveline = -log10(5e-06), genomewideline = -log10(5e-08))
dev.off()


####Write table ##########################
mirsnps_pso <- pso_mirs[pso_mirs$P < 5e-06,]

write.table(mirsnps_pso, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/sign_mirsnps_pso.txt", sep = '\t', 
            quote = F, row.names = F)

nonhla_mirsnps_pso <- nonhla_pso[nonhla_pso$P < 5e-06,]

write.table(nonhla_mirsnps_pso, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/sign_nonhla_mirsnps_pso.txt", sep = '\t', 
            quote = F, row.names = F)


###############################################
######T1D###############

t1d_mirs <- read.table("/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/tables/t1d_mirsnps.txt", header = T)
head(t1d_mirs)

t1d_mirs <- subset(t1d_mirs, select = c("Chr", "Position", "PValue", "Marker"))
head(t1d_mirs)
colnames(t1d_mirs) <- c("CHR", "BP", "P", "SNP")

#Get rid off CHR X and Y
t1d_mirs <- t1d_mirs[!t1d_mirs$CHR == "X",]
t1d_mirs <- t1d_mirs[!t1d_mirs$CHR == "Y",]

#Check dataft1dme
lapply(t1d_mirs, data.class)

#Convert column to numeric
t1d_mirs$CHR = as.numeric(as.chat1dcter(t1d_mirs$CHR))

tiff("/Users/lgmh/Dropbox/shared_immune_mirsnps/figs/manhatan/mh_plot_mirnas_t1d.tif", width = 890, height = 675, pointsize = 19)
manhattan(t1d_mirs,main = "miRSNPs T1D",cex = 0.9,cex.axis = 0.8, ylim  = c(0,50),
          col = c("blue4", "deepskyblue"), suggestiveline = -log10(5e-06), genomewideline = -log10(5e-08))
dev.off()

###########Non HLA########################################
#GRCh37
#chr6:28477797-33448354 

hla_t1d <- t1d_mirs[ t1d_mirs$CHR == "6",]
hla_t1d <- hla_t1d[hla_t1d$BP > 25780582 & hla_t1d$BP < 33448354,] 

nonhla_t1d <- t1d_mirs[!(t1d_mirs$SNP %in% hla_t1d$SNP),]

tiff("/Users/lgmh/Dropbox/shared_immune_mirsnps/figs/manhatan/mh_plot_nohla_mirnas_t1d.tif",  width = 890, height = 675, pointsize = 19)
manhattan(nonhla_t1d,main = "Non-HLA miRSNPs T1D",cex = 0.9,cex.axis = 0.9, ylim = c(0,50),
          col = c("blue4", "deepskyblue"), suggestiveline = -log10(5e-06), genomewideline = -log10(5e-08))
dev.off()

####Write table ##########################
mirsnps_t1d <- t1d_mirs[t1d_mirs$P < 5e-06,]

write.table(mirsnps_t1d, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/sign_mirsnps_t1d.txt", sep = '\t', 
            quote = F, row.names = F)

nonhla_mirsnps_t1d <- nonhla_t1d[nonhla_t1d$P < 5e-06,]

write.table(nonhla_mirsnps_t1d, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/sign_nonhla_mirsnps_t1d.txt", sep = '\t', 
            quote = F, row.names = F)

############################################
##########SLE###############

sle_mirs <- read.table("/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/tables/sle_mirsnps.txt", header = T)
head(sle_mirs)

sle_mirs <- subset(sle_mirs, select = c("Chr", "Position", "PValue", "Marker"))
head(sle_mirs)
colnames(sle_mirs) <- c("CHR", "BP", "P", "SNP")

#Get rid off CHR X and Y
sle_mirs <- sle_mirs[!sle_mirs$CHR == "X",]
sle_mirs <- sle_mirs[!sle_mirs$CHR == "Y",]

#Check datafsleme
lapply(sle_mirs, data.class)

#Convert column to numeric
sle_mirs$CHR = as.numeric(as.character(sle_mirs$CHR))

tiff("/Users/lgmh/Dropbox/shared_immune_mirsnps/figs/manhatan/mh_plot_mirnas_sle.tif",width = 890, height = 675, pointsize = 19 )
manhattan(sle_mirs,main = "miRSNPs SLE",cex = 0.9,cex.axis = 0.8, ylim = c(0, 60),
          col = c("blue4", "deepskyblue"), suggestiveline = -log10(5e-06), genomewideline = -log10(5e-08))
dev.off()

###########Non HLA########################################
#GRCh37
#chr6:28477797-33448354 

hla_sle <- sle_mirs[ sle_mirs$CHR == "6",]
hla_sle <- hla_sle[hla_sle$BP > 25780582 & hla_sle$BP < 33448354,] 

nonhla_sle <- sle_mirs[!(sle_mirs$SNP %in% hla_sle$SNP),]

tiff("/Users/lgmh/Dropbox/shared_immune_mirsnps/figs/manhatan/mh_plot_nohla_mirnas_sle.tif",  width = 890, height = 675, pointsize = 19)
manhattan(nonhla_sle,main = "Non-HLA miRSNPs SLE",cex = 0.9,cex.axis = 0.9, ylim = c(0,50),
          col = c("blue4", "deepskyblue"), suggestiveline = -log10(5e-06), genomewideline = -log10(5e-08))
dev.off()

####Write table ##########################
mirsnps_sle <- sle_mirs[sle_mirs$P < 5e-06,]

write.table(mirsnps_sle, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/sign_mirsnps_sle.txt", sep = '\t', 
            quote = F, row.names = F)

nonhla_mirsnps_sle <- nonhla_sle[nonhla_sle$P < 5e-06,]

write.table(nonhla_mirsnps_sle, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/sign_nonhla_mirsnps_sle.txt", sep = '\t', 
            quote = F, row.names = F)


########################################
####UC##################
uc_mirs <- read.table("/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/tables/uc_mirsnps.txt", header = T)
head(uc_mirs)

uc_mirs <- subset(uc_mirs, select = c("Chr", "Position", "PValue", "Marker"))
head(uc_mirs)
colnames(uc_mirs) <- c("CHR", "BP", "P", "SNP")

#Get rid off CHR X and Y
uc_mirs <- uc_mirs[!uc_mirs$CHR == "X",]
uc_mirs <- uc_mirs[!uc_mirs$CHR == "Y",]

#Check dataframe
lapply(uc_mirs, data.class)

#Convert column to numeric
uc_mirs$CHR = as.numeric(as.character(uc_mirs$CHR))

tiff("/Users/lgmh/Dropbox/shared_immune_mirsnps/figs/manhatan/mh_plot_mirnas_uc.tif", width = 890, height = 675, pointsize = 19)
manhattan(uc_mirs,main = "miRSNPs UC",cex = 0.9,cex.axis = 0.8,ylim = c(0, 35),
          col = c("blue4", "deepskyblue"), suggestiveline = -log10(5e-06), genomewideline = -log10(5e-08))
dev.off()

###########Non HLA########################################
#GRCh37
#chr6:28477797-33448354 

hla_uc <- uc_mirs[ uc_mirs$CHR == "6",]
hla_uc <- hla_uc[hla_uc$BP > 25780582 & hla_uc$BP < 33448354,] 

nonhla_uc <- uc_mirs[!(uc_mirs$SNP %in% hla_uc$SNP),]

tiff("/Users/lgmh/Dropbox/shared_immune_mirsnps/figs/manhatan/mh_plot_nonhla_mirnas_uc.tif",  width = 890, height = 675, pointsize = 19)
manhattan(nonhla_uc,main = "Non-HLA miRSNPs UC",cex = 0.9,cex.axis = 0.9, ylim = c(0, 35), 
          col = c("blue4", "deepskyblue"), suggestiveline = -log10(5e-06), genomewideline = -log10(5e-08))
dev.off()

####Write table ##########################
mirsnps_uc <- uc_mirs[uc_mirs$P < 5e-06,]

write.table(mirsnps_uc, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/sign_mirsnps_uc.txt", sep = '\t', 
            quote = F, row.names = F)

nonhla_mirsnps_uc <- nonhla_uc[nonhla_uc$P < 5e-06,]

write.table(nonhla_mirsnps_uc, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/sign_nonhla_mirsnps_uc.txt", sep = '\t', 
            quote = F, row.names = F)

########################################
####MS##################
ms_mirs <- read.table("/Users/lgmh/Documents/rodrigo/IdeaS/Immunobase/tables/ms_mirsnps.txt", header = T)
head(ms_mirs)

ms_mirs <- subset(ms_mirs, select = c("Chr", "Position", "PValue", "Marker"))
head(ms_mirs)
colnames(ms_mirs) <- c("CHR", "BP", "P", "SNP")

#Get rid off CHR X and Y
ms_mirs <- ms_mirs[!ms_mirs$CHR == "X",]
ms_mirs <- ms_mirs[!ms_mirs$CHR == "Y",]

#Check dataframe
lapply(ms_mirs, data.class)

#Convert column to numeric
ms_mirs$CHR = as.numeric(as.character(ms_mirs$CHR))

tiff("/Users/lgmh/Dropbox/shared_immune_mirsnps/figs/manhatan/mh_plot_mirnas_ms.tif", width = 890, height = 675, pointsize = 19)
manhattan(ms_mirs,main = "miRSNPs MS",cex = 0.9,cex.axis = 0.8,ylim = c(0, 290),
          col = c("blue4", "deepskyblue"), suggestiveline = -log10(5e-06), genomewideline = -log10(5e-08))
dev.off()

###########Non HLA########################################
#GRCh37
#chr6:28477797-33448354 

hla_ms <- ms_mirs[ ms_mirs$CHR == "6",]
hla_ms <- hla_ms[hla_ms$BP > 25780582 & hla_ms$BP < 33448354,] 

nonhla_ms <- ms_mirs[!(ms_mirs$SNP %in% hla_ms$SNP),]

tiff("/Users/lgmh/Dropbox/shared_immune_mirsnps/figs/manhatan/mh_plot_nonhla_mirnas_ms.tif",  width = 890, height = 675, pointsize = 19)
manhattan(nonhla_ms,main = "Non-HLA miRSNPs MS",cex = 0.9,cex.axis = 0.9, ylim = c(0, 18), 
          col = c("blue4", "deepskyblue"), suggestiveline = -log10(5e-06), genomewideline = -log10(5e-08))
dev.off()

####Write table ##########################
mirsnps_ms <- ms_mirs[ms_mirs$P < 5e-06,]

write.table(mirsnps_ms, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/sign_mirsnps_ms.txt", sep = '\t', 
            quote = F, row.names = F)

nonhla_mirsnps_ms <- nonhla_ms[nonhla_ms$P < 5e-06,]

write.table(nonhla_mirsnps_ms, "/Users/lgmh/Dropbox/shared_immune_mirsnps/tables/top_mirsnps/sign_nonhla_mirsnps_ms.txt", sep = '\t', 
            quote = F, row.names = F)

