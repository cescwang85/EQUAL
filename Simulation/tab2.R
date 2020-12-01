rm(list = ls())
setwd('~/Dropbox/Share/Cheng&Yan/matinv/Simulation/')
library('MASS')
library('stargazer')
set.seed(123)

TT<-array(0,dim=c(300,15,3));
for (id in 1:3)
{
  TT[,,id]=readRDS(paste('error',id,'-tab2','.rds',sep =''))[,1:15]
}


Re<-matrix(0,27,5)
dig=3
ff<-function(obj){
  m<-format(round(apply(obj,2,mean,na.rm=TRUE),dig),nsmall =dig)
  s<-format(round(apply(obj,2,sd,na.rm=TRUE),dig),nsmall =dig)
  l<-rep('(',length(m));
  r<-rep(')',length(m));
  R1=paste(m,l,s,r,sep='')
  return(R1)
}

for (k in 1:3)
{Re[1:3+(k-1)*3,]=matrix(ff(TT[1:100+(k-1)*100,,1]),ncol =5,byrow =TRUE);
Re[1:3+(k-1)*3+9,]=matrix(ff(TT[1:100+(k-1)*100,,2]),ncol =5,byrow =TRUE);
Re[1:3+(k-1)*3+18,]=matrix(ff(TT[1:100+(k-1)*100,,3]),ncol =5,byrow =TRUE);
}

colnames(Re)<-c('Loss1','Loss2','Loss3','Loss4','min-Eigen');
Re<-cbind(rep(c('EQUAL','EQUALs','glasso'),9),Re)
Re
stargazer(Re)
