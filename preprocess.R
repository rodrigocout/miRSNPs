###############################################
########Preprocess
###############################################

########################
########-mirsnpscore
########################

#---clean table
# system("cat databases/mirsnpscore.txt | awk '{print $1,$2,$3,$4,$5,$6, $7,$8,$9,$10,$11}' > ~/Dropbox/BED_files_miRSNPS/databases/newmirsnpscore.txt")
# 
# #---read and process table
# mirsnpscore<-read.table(file = "databases//newmirsnpscore.txt",sep="",header =T,stringsAsFactors = FALSE)
# colnames(mirsnpscore)<-c("chr","Gene", "RefGene","Pos","Strand","miRNA","DeltaS","sc1","sc2","RsID","Pos")
# save(mirsnpscore,file="databases//newmirsnpscore.RData")

########################
########-mirsnpscore REV
########################

mirsnpscore<-read.table(file = "~/Dropbox/BED_files_miRSNPS/databases/mirsnpscore.txt",sep="#",header=F,stringsAsFactors=FALSE)
dtmir<-NULL
dtmirsnps<-list()
for(i in 1:nrow(mirsnpscore)){
  tp<-strsplit(mirsnpscore[i,],"\t")[[1]]
  #---
  tp1<-tp[1:9]
  mir<-tp1[6]
  #---
  tp2<-tp[10:length(tp)]
  vecmir<-NULL
  for(j in 1:length(tp2)){
    tp3<-strsplit(tp2[j]," ")[[1]]
    vecmir<-c(vecmir,tp3[1])
  }
  dtmir<-rbind(dtmir,c(tp1,vecmir[1]))
  dtmirsnps[[i]]<-vecmir
}
#---
rownames(dtmir)<-NULL
colnames(dtmir)<-c("chr","Gene", "RefGene","Pos","Strand","miRNA","DeltaS","sc1","sc2","SNP")
dtmir<-data.frame(dtmir,stringsAsFactors = FALSE)
for(i in 7:9)dtmir[,i]<-as.numeric(dtmir[,i])
#---
save(dtmir,dtmirsnps,file="~/Dropbox/BED_files_miRSNPS/databases/newmirsnpscoreREV.RData")

########################
########-mirSnpTarget
########################
mirsnptarget<-read.table(file = "~/Dropbox/BED_files_miRSNPS/databases/MirSNPInTarget.txt",sep="",header =F,stringsAsFactors = FALSE)
tmp<-paste("rs",mirsnptarget[,1],sep="")
mirsnptarget[,1]<-tmp
mirsnptarget<-mirsnptarget[,-c(5:9,12,13,15:20)]
colnames(mirsnptarget)<-c("SNP","Gene","miRNA","transcript","score","energy","conservation")
save(mirsnptarget,file="~/Dropbox/BED_files_miRSNPS/databases/mirsnptarget.RData")
#---
annot<-paste(mirsnptarget$SNP,mirsnptarget$Gene,mirsnptarget$miRNA,sep ="/")
mirsnptarget$annot<-annot
idx<-sort.list(mirsnptarget$score,decreasing = TRUE)
mirsnptarget<-mirsnptarget[idx,]
#---
annot<-unique(annot)
idx<-match(annot,mirsnptarget$annot)
mirsnptarget<-mirsnptarget[idx,]
all(mirsnptarget$annot==annot)
save(mirsnptarget,file="~/Dropbox/BED_files_miRSNPS/databases/mirsnptarget.RData")


########################
########-Polymirts
########################
polymirts<-read.table(file = "~/Dropbox/BED_files_miRSNPS/databases/Polymirts.txt", sep="\t", quote = '',header=TRUE,stringsAsFactors = FALSE)
polymirts<-polymirts[,c(1:11,28:29)]
polymirts<-polymirts[,c(3,10,11,12,13)]
#---
n=nrow(polymirts)
Allele1miR<-strsplit(polymirts[1:n,"Allele1miR"],"|",fixed=T)
cs_diff_Allele1Site<-strsplit(polymirts[1:n,"cs_diff_Allele1Site"],"|",fixed=T)
cs_diff_Allele1<-lapply(1:length(cs_diff_Allele1Site), function(i){
  tp<-as.numeric(cs_diff_Allele1Site[[i]])
  names(tp)<-Allele1miR[[i]]
  tp
})
#---
Allele2miR<-strsplit(polymirts[1:n,"Allele2miR"],"|",fixed=T)
cs_diff_Allele2Site<-strsplit(polymirts[1:n,"cs_diff_Allele2Site"],"|",fixed=T)
cs_diff_Allele2<-lapply(1:length(cs_diff_Allele2Site), function(i){
  tp<-as.numeric(cs_diff_Allele2Site[[i]])
  names(tp)<-Allele2miR[[i]]
  tp
})
#---#---#---#---#---#---#---
cs_diff_Allele<-lapply(1:length(cs_diff_Allele1), function(i){
  tp1<-cs_diff_Allele1[[i]]
  tp2<-cs_diff_Allele2[[i]]
  tp<-sort(c(tp1,tp2),decreasing=FALSE)
  nms<-unique(names(tp))
  idx<-match(nms,names(tp))
  tp[idx]
})

#---#---#---#---#---#---#---
snpmap<-unlist(lapply(cs_diff_Allele,length))
names(snpmap)<-polymirts[,1]
snpmap<-rep(names(snpmap),times=snpmap)
#---
tp<-unlist(cs_diff_Allele)
polymirts_merge<-data.frame(SNPID=snpmap,miR=names(tp),cs_diff=tp)

#---#---#---#---#---#---#---
save(polymirts_merge,file="~/Dropbox/BED_files_miRSNPS/databases/polymirts_merge.RData")
