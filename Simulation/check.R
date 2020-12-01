#Numerical verification of Proposition 1

rm(list = ls())
library('MASS')
setwd('~/Dropbox/Share/Cheng&Yan/matinv/Simulation/')
set.seed(1985)
n=5;
p=11;
m=min(n,p);
rho=0.38;
X<-matrix(rnorm(n*p),nrow =n,ncol=p);
S<-t(X)%*%X/n;

obj<-eigen(S);
la<-obj$values[1:m];
U<-obj$vectors[,1:m];
Lambda<-diag(la);

Lambda1<-diag(la/(la+rho));
Lambda2<-diag(la/(la+2*rho));
Lambda3<-matrix(0,m,m)
for (i in 1:m){
  for (j in 1:m)
Lambda3[i,j]=la[i]*la[j]*(la[i]+la[j]+4*rho)/(la[i]+2*rho)/(la[j]+2*rho)/(la[i]+la[j]+2*rho)}

dim(U)
dim(Lambda1)
dim(Lambda2)
dim(Lambda3)
max(abs(S-U%*%Lambda%*%t(U)))


###Check 
A<-solve(S+rho*diag(p))
B<-solve(kronecker(S,diag(p))/2+kronecker(diag(p),S)/2+rho*diag(p^2))
A1<-diag(p)/rho-U%*%Lambda1%*%t(U)/rho;
B1<-diag(p^2)/rho-kronecker(U%*%Lambda2%*%t(U),diag(p))/rho-kronecker(diag(p),U%*%Lambda2%*%t(U))/rho+1/rho*kronecker(U,U)%*%diag(as.vector(Lambda3))%*%kronecker(t(U),t(U));

#Verification of formula (12)
max(abs(A-A1))
#Verification of formula (13)
max(abs(B-B1))



  