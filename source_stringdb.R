###################################################
### Combine scores code based on
### STRING naive Bayesian fashion
### von Mering et al., Nucleic Acids Res. 2005 Jan 1; 
### 33(Database issue): D433–D437.
###################################################
# OBS: scores should be transformed to 0-1 scale!
combinescores<-function(dat, evidences="all", confLevel=0.7){
  if(evidences[1]=="all"){
    edat<-dat[,-c(1,2,ncol(dat))]
  } else {
    if(!all(evidences%in%colnames(dat))){
      stop("NOTE: one or more 'evidences' not listed in 'dat' colnames!")
    }
    edat<-dat[,evidences]
  }
  edat<-1-edat
  sc<-sapply(1:nrow(edat),function(i){
    tp<-edat[i,]
    1-prod(tp)
  })
  dat<-cbind(dat[,c(1,2)],combined_score=sc)
  idx<-dat$combined_score>confLevel
  dat<-dat[idx,]
  return(dat)
}
