rm(list = ls())

setwd('~/Dropbox/Share/Cheng&Yan/matinv/Simulation/')

n=200
id=3; ## Select Omega
p=1000
  source('ffun.R')
  Omega<-array(0,dim=c(p,p,3));
  Omega[,,1]<-toeplitz(0.5^(1:p-1))
  Omega[,,2]<-solve(toeplitz(0.5^(1:p-1)));
  for (k in 1:(p/5))
  {Omega[1:5+5*(k-1),1:5+5*(k-1),3]=runif(1,1,5)*(matrix(0.5,5,5)+0.5*diag(5))}
  X=mvrnorm(n,rep(0,p),solve(Omega[,,id]))
  AA<-EQUAL(X)
re<-matrix(0,5,50)
for (k in 1:50)
{re[,k]=loss(AA$Omega[[k]],Omega[,,3])}
plot(AA$lambda, re[3,])
