rm(list = ls())
set.seed(1985)
setwd('~/Dropbox/Share/Cheng&Yan/matinv/Simulation/')
library('MASS')
library('Rcpp')
library('Matrix')
library('fastclime')
library('glasso')
library('BigQuic')
library('scio')
library('EQUAL')
sourceCpp('hp.cpp')
id=3; ## Select Omega
ttime<-function(obj)
{t1<-proc.time();
obj
t2<-as.vector(proc.time()-t1);
return(sum(t2[1:2]))}

ff<-function(n,p){
  Omega<-array(0,dim=c(p,p,3));
  Omega[,,1]<-toeplitz(0.5^(1:p-1))
  Omega[,,2]<-solve(toeplitz(0.5^(1:p-1)));
  for (k in 1:(p/5))
  {Omega[1:5+5*(k-1),1:5+5*(k-1),3]=runif(1,0.5,5)*(matrix(0.5,5,5)+0.5*diag(5))}
  Omega[,,3]=Omega[,,3]/mean(diag(Omega[,,3]))
  X=mvrnorm(n,rep(0,p),solve(Omega[,,id]))
  Sn<-cov(X);
  nlambda=50
  lambda<-exp(seq(0,log(sqrt(log(p)/n)),length.out =nlambda))*max(abs(Sn))
  tt<-rep(0,8)
  tt[1]<-ttime(obj0<-fastclime(Sn,lambda.min=0.1, nlambda =50))
  tt[2]<-ttime(obj1<-glassopath(Sn,rholist=lambda,penalize.diagonal=FALSE,trace =0))
  tt[3]<-ttime(obj2<-BigQuic(X=X,lambda =lambda,use_ram =TRUE))
  tt[4]<-ttime(obj3<-gOmega(X,lambda =lambda))
  tt[5]<-ttime(obj4<-sciopath(Sn,lambdalist =lambda,pen.diag =FALSE))
  tt[6]<-ttime(obj5<-EQUAL(X,lambda =lambda,type=FALSE))
  tt[7]<-ttime(obj6<-dtrace(X,lambda =lambda))
  tt[8]<-ttime(obj7<-EQUAL(X,lambda =lambda,type=TRUE))
  print(c(n,p,sum(obj5$niter),sum(obj7$niter),sum(obj3$niter),tt))}
nn<-200;
pp<-c(100,200,400,800,1600);
Re<-mapply(ff,nn,rep(pp,each=100))
saveRDS(Re,file=paste('time',id,'-tab1','.rds',sep =''))
