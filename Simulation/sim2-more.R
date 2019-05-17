rm(list = ls())
library('MASS')
library('Rcpp')
library('Matrix')
library('EQUAL')
setwd('~/Dropbox/Share/Cheng&Yan/matinv/Simulation/')
set.seed(1985)
n=200;
p=500;
Omega<-array(0,dim=c(p,p,3));
Omega[,,1]<-toeplitz(0.5^(1:p-1))
Omega[,,2]<-solve(toeplitz(0.5^(1:p-1)));
for (k in 1:(p/5))
{Omega[1:5+5*(k-1),1:5+5*(k-1),3]=runif(1,0.5,5)*(matrix(0.5,5,5)+0.5*diag(5))}
Omega[,,3]=Omega[,,3]/mean(diag(Omega[,,3]))

for (id in 1:3)## Select Omega
{X=mvrnorm(n,rep(0,p),solve(Omega[,,id]))
Ome1<-EQUAL(X,lambda.min=sqrt(log(p)/n),type =FALSE)
Ome2<-EQUAL(X,lambda.min=sqrt(log(p)/n),type =TRUE)
lambda<-Ome1$lambda;
Re<-matrix(0,2,length(lambda));
for (k in 1:length(lambda))
{
  Re[1,k]=min(eigen(Ome1$Omega[[k]])$values);
    Re[2,k]=min(eigen(Ome2$Omega[[k]])$values)
}
pdf(paste('fig/','sim2-',id,'.pdf',sep=''))
plot(lambda,Re[1,],ylim=c(min(Re)*0.95,max(Re)*1.05),ylab='minEigen',type='o',pch=1)
lines(lambda,Re[2,],type='o',pch=8)
legend('bottomright',legend =c('EQUAL', 'EQUALs'),pch=c(1,8),lty=c(1,1))
dev.off()
}

