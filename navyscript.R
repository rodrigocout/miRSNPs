#########################
### Load
#########################
#---data
load (file = "~/Dropbox/BED_files_miRSNPS/databases/navyData.RData")
#---source
source ("~/Dropbox/BED_files_miRSNPS/scripts/source_stringdb.R")
source ("~/Dropbox/BED_files_miRSNPS/scripts/runboxcox.R")
#########################
### Process
#########################

#--- normalize in z
#--- Z score
datanavy$mirsnpscore <- (datanavy$mirsnpscore - median (datanavy$mirsnpscore, na.rm = TRUE) ) / sd (datanavy$mirsnpscore, na.rm = TRUE)
datanavy$mirsnptarget <- (datanavy$mirsnptarget - median (datanavy$mirsnptarget, na.rm = TRUE)) / sd (datanavy$mirsnptarget, na.rm = TRUE)
datanavy$polymirts <- (datanavy$polymirts - median (datanavy$polymirts, na.rm = TRUE)) / sd (datanavy$polymirts, na.rm = TRUE)

#--- boxcox transformation
datanavy$mirsnpscore<-runboxcox(datanavy$mirsnpscore)
datanavy$mirsnptarget<-runboxcox(datanavy$mirsnptarget)
datanavy$polymirts<-runboxcox(datanavy$polymirts)

#---
plot(density(datanavy$mirsnpscore,na.rm = TRUE))
plot(density(datanavy$mirsnptarget,na.rm = TRUE))
plot(density(datanavy$polymirts,na.rm = TRUE))

# #---pnorm
datanavy$mirsnpscore <- 1 - pnorm(datanavy$mirsnpscore, mean=mean(datanavy$mirsnpscore, na.rm=TRUE), sd=sd(datanavy$mirsnpscore, na.rm = TRUE), lower.tail = FALSE)
datanavy$mirsnptarget <- 1 - pnorm(datanavy$mirsnptarget, mean=mean(datanavy$mirsnptarget, na.rm=TRUE), sd=sd(datanavy$mirsnptarget, na.rm = TRUE), lower.tail = FALSE)
datanavy$polymirts <- pnorm(datanavy$polymirts, mean=mean(datanavy$polymirts, na.rm=TRUE), sd=sd(datanavy$polymirts, na.rm = TRUE), lower.tail = FALSE)

#--- NA = 0
datanavy [is.na (datanavy)] <- 0

#########################
### Navy
#########################
navyD <- combinescores(datanavy, evidences = c ("mirsnpscore", "mirsnptarget", "polymirts"))
  
save (navyD, file = "~/Dropbox/BED_files_miRSNPS/databases/navyResultsData.RData")
