rm(list = ls())
set.seed(1985)
dig=3
setwd('~/Dropbox/Share/Cheng&Yan/matinv/Simulation/')
library('stargazer')

TT<-array(0,dim=c(8,500,3));
for (id in 1:3)
{
TT[,,id]=readRDS(paste('time',id,'-tab1','.rds',sep =''))[-(1:5),]
}

Re<-matrix(0,24,5)
ff<-function(obj){
  m<-format(round(apply(obj,1,mean,na.rm=TRUE),dig),nsmall =dig)
  s<-format(round(apply(obj,1,sd,na.rm=TRUE),dig),nsmall =dig)
  l<-rep('(',length(m));
  r<-rep(')',length(m));
  R1=paste(m,l,s,r,sep='')
  return(R1)
}
for (k in 1:5)
{Re[1:8,k]=ff(TT[,1:100+(k-1)*100,1]);
Re[9:16,k]=ff(TT[,1:100+(k-1)*100,2]);
Re[17:24,k]=ff(TT[,1:100+(k-1)*100,3]);
}
colnames(Re)<-c('p=100','p=200','p=400','p=800','p=1600');
Re<-cbind(rep(c('CLIME','glasso','BigQuic','glasso-ADMM','SCIO','EQUAL','D-trace','	EQUALs'),3),Re)
stargazer(Re)

