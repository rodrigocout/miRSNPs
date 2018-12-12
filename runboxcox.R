#----powerTransform
library(car)
runboxcox<-function(obsdist){
  if(shtest(obsdist)){
    minval<-min(obsdist, na.rm = TRUE)
    minval<-ifelse(minval<=0,abs(minval)+1,minval)
    nullm<-obsdist+minval
    l<-coef(powerTransform(nullm), round=TRUE)
    ptdat<-bcPower(nullm,l)
    ptdist<-(ptdat-median(ptdat, na.rm = TRUE))/sd(ptdat,na.rm = TRUE)
  } else {
    ptdist<-obsdist
  }
  #pvals<-pnorm(ptdist, lower.tail=FALSE)
  ptdist
}
#-------------------------------------------------------------------------
shtest<-function(obsdist, nd=5000){
  nnull<-length(obsdist)
  qt<-quantile(obsdist,probs=seq(0,1,length.out=nd),names=FALSE, na.rm = TRUE)
  shapiro.test(qt)$p.value<0.05
}
