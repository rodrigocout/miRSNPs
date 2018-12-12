setwd("/Users/Rodrigo/Dropbox/shared_immune_mirsnps/tables/files_for_circosplot")

mirs <- read.table("lista_microRNAs.txt")
colnames(mirs) <- c("mirnas")

bed <- read.table("mirnas_bed_file_hg38.txt")
colnames(bed) <- c("chr", "start", "end", "mirnas")


mir_bed <- merge(bed, mirs, by = c("mirnas"))
mir_bed <- mir_bed[order(mir_bed$chr),]
mir_bed <- mir_bed[,c(2,3,4,1)]

write.table(mir_bed, "mirnas_bed_hg38_merged.txt", sep = '\t', row.names = F, quote = F)

mirs$match <- match(mirs$mirnas, bed$mirnas, nomatch=0)

non.matched <- mirs[mirs$match == "0",]

write.table(non.matched, "non_matched_mirnas_bed.txt", sep = '\t', quote = F, row.names = F)
