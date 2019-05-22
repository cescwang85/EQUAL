rm(list = ls())
set.seed(1985)
setwd('~/Dropbox/Share/Cheng&Yan/matinv/simulation/')
library('MASS')
library('Rcpp')
library('Matrix')
sourceCpp('hp.cpp')
library('EQUAL')
soft0<-function(x,lambda){(x>lambda)*(x-lambda)+(x<(-lambda))*(x+lambda)}
mcp<-function(x,lambda){soft0(2*lambda,abs(x))/2}
scad<-function(x,lambda){lambda*(abs(x)<=lambda)+(abs(x)>lambda)*soft0(3.7*lambda,abs(x))/2.7}
loss<-function(A,B)
{aa<-eigen(solve(B,A))$values;
p=ncol(B);
a1<-norm(A-B,type='F');
a2<-norm(A-B,type='2');
a3<-mean(aa-log(aa)-1);
a4<-sum(diag(t(A)%*%solve(B,A)))/2-sum(diag(A))+sum(diag(B))/2;
return(c(a1/sqrt(p),a2,sqrt(a3),sqrt(a4/p)))}

p=200
n=200
id=3; ## Select Omega


Omega<-array(0,dim=c(p,p,3));
Omega[,,1]<-toeplitz(0.5^(1:p-1))
Omega[,,2]<-solve(toeplitz(0.5^(1:p-1)));
for (k in 1:(p/5))
{Omega[1:5+5*(k-1),1:5+5*(k-1),3]=runif(1,0.5,5)*(matrix(0.5,5,5)+0.5*diag(5))}
Omega[,,3]=Omega[,,3]/mean(diag(Omega[,,3]))
X=mvrnorm(n,rep(0,p),solve(Omega[,,id]))
Sn<-cov(X);
max(cov(Sn))
nlambda=50
lambda<-exp(seq(log(0.5),log(0.1),length.out =nlambda)) ##for Omega 1 and Omega2
#lambda<-exp(seq(log(0.4),log(0.08),length.out =nlambda)) ##for Omega 3
#lambda<-exp(seq(log(0.1),0,length.out =nlambda))
#lambda<-exp(seq(log(2*log(p)/n),log(0.5),length.out=NN))

obj1<-EQUAL(X,lambda =lambda)
lam1<-array(0,dim=c(p,p,nlambda));
lam2<-array(0,dim=c(p,p,nlambda));
for (k in 1:nlambda){
lam1[,,k]=mcp(as.matrix(obj1$Omega[[k]]),lambda[k]);
lam2[,,k]=scad(as.matrix(obj1$Omega[[k]]),lambda[k]);
}
obj2<-sOmega(X,lambda=lam1)
obj3<-sOmega(X,lambda=lam2)

Re<-matrix(0,nlambda,12)
for (k in 1:nlambda)
{Re[k,1:4]=loss(obj1$Omega[[k]],Omega[,,id]);
Re[k,5:8]=loss(obj2$Omega[[k]],Omega[,,id]);
Re[k,9:12]=loss(obj3$Omega[[k]],Omega[,,id]);
}
for (tt in 1:4){
pdf(paste('fig/',id,'-','loss',tt,'-sim3','.pdf',sep=''))
plot(lambda,Re[,tt+4],type='o',ylab=paste('loss',tt,sep =''),pch=16)
points(lambda,Re[,tt+8],type='o',pch=1)
points(lambda,Re[,tt],type='o',pch=8)
legend('bottomright',legend =c('MCP','SCAD','LASSO'),pch=c(16,1,8),lty=c(1,1,1))
dev.off()
}
#AA=matrix(apply(Re,2,min),nrow=5)
#rownames(AA)=c('loss1','loss2','loss3','loss4','min-Eigen')
#colnames(AA)=c('LASSO','SCAD','MCP')
##Optimal loss
#AA

