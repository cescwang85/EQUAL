rm(list = ls())
set.seed(123)
library('MASS')
library('Rcpp')
library('Matrix')
library('EQUAL')
n=200
p=100
Omega<-toeplitz(0.5^(1:p-1))
X=mvrnorm(n,rep(0,p),solve(Omega))
aa<-EQUAL(X)
bb<-EQUAL(X,type=FALSE);

obj1<-CVEQUAL(X);
obj2<-CVEQUAL(X,type=FALSE)
obj1$Omega[1:10,1:10]
obj2$Omega[1:10,1:10]