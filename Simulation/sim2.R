rm(list = ls())
library(doParallel)
library('doRNG',quietly =TRUE)
setwd('~/Dropbox/Share/Cheng&Yan/matinv/Simulation/')

cl <- makeCluster(4)
registerDoParallel(cl)
n=200
id=3; ## Select Omega
ff<-function(p){
  source('ffun.R')
  Omega<-array(0,dim=c(p,p,3));
  Omega[,,1]<-toeplitz(0.5^(1:p-1))
  Omega[,,2]<-solve(toeplitz(0.5^(1:p-1)));
  for (k in 1:(p/5))
  {Omega[1:5+5*(k-1),1:5+5*(k-1),3]=runif(1,0.5,5)*(matrix(0.5,5,5)+0.5*diag(5))}
  Omega[,,3]=Omega[,,3]/mean(diag(Omega[,,3]))
  X=mvrnorm(n,rep(0,p),solve(Omega[,,id]))
  t1<-ttime(Ome1<-CVEQUAL(X,lambda.min=sqrt(log(p)/n),type =FALSE)$Omega)
  t2<-ttime(Ome2<-CVEQUAL(X,lambda.min=sqrt(log(p)/n),type =TRUE)$Omega)
  t3<-ttime(Ome3<-cvglasso(X,lambda.min=sqrt(log(p)/n)))
  return(c(loss(Ome1,Omega[,,id]),loss(Ome2,Omega[,,id]),loss(Ome3,Omega[,,id]),t1,t2,t3))
}
Re<-foreach(p=c(rep(500,100),rep(1000,100),rep(2000,100)),.combine ='rbind',.options.RNG=123)%dorng% ff(p)
stopCluster(cl)
saveRDS(Re,file=paste('error',id,'-tab2','.rds',sep =''))
