#################################
#########Analises Max###########
###############################
setwd('/Users/Rodrigo/Documents/Pos_doc/IdeaS/fine_mapping/exploratory/')

##Carregar os dados
eqtl <- read.delim('Scores_eQTLs.FCN1.txt', header = T) # 188 eQTL mais LD SNPs regiao FCN1
non_eqtl <- read.delim('scores_SNPs_totais.txt', header = T) # Todos os SNPs nao eQTL encontrados na regiao FCN1

hist(eqtl$total)
hist(non_eqtl$total)
class(eqtl)
class(non_eqtl)
dim(eqtl)
dim(non_eqtl)

#Primeiro calcular a media 
eqtl_mean <- mean(eqtl$total)

non_eqtl_mean_chrm_st <- mean(non_eqtl$total)


#Depois calcular o desvio padrao (sd)
eqtl_sd <- sd(eqtl$total)*sqrt((length(eqtl$total)-1)/(length(eqtl$total)))

non_eqtl_sd <- sd(non_eqtl$total)*sqrt((length(non_eqtl$total)-1)/(length(non_eqtl$total)))


#z-score calculation

z <- scale(eqtl$total, center = T, scale = T) #usando a func scale

#ou sem a func scale
z_eqtl <- (eqtl$total - eqtl_mean)/eqtl_sd 

#colocar a coluna com os z-scores na matrix
eqtl_z <- cbind(eqtl,z) 

#distribution
hist(z)

#salvar em um arquivo txt
write.table(eqtl_z, "eQTL_Zscore.txt", sep = '\t')



