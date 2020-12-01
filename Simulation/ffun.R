library('MASS')
library('Rcpp')
library('Matrix')
library('fastclime')
library('glasso')
library('BigQuic')
library('scio')
library('EQUAL')
ttime<-function(obj)
{t1<-proc.time();
obj
t2<-as.vector(proc.time()-t1);
return(sum(t2[1:2]))}
loss<-function(A,B)
{aa<-eigen(solve(B,A))$values;
p=ncol(B);
  a1<-norm(A-B,type='F');
a2<-norm(A-B,type='2');
a3<-mean(aa-log(aa)-1);
a4<-sum(diag(t(A)%*%solve(B,A)))/2-sum(diag(A))+sum(diag(B))/2;
return(c(a1/sqrt(p),a2,sqrt(a3),sqrt(a4/p),min(eigen(A)$values)))}

fglasso<-function(X,lambda,err=1e-4,maxIter =1e4)
{
  S<-cov(X);
  obj<-glassopath(S, rholist=lambda,penalize.diagonal=FALSE,trace=0,thr=err,maxit =maxIter);
  nlambda<-length(lambda);
  Omega<-NULL;
  for (k in 1:nlambda)
  {Omega=c(Omega,list(as(obj$wi[,,k], "sparseMatrix")))}
  return(list(Omega=Omega,lambda=obj$rholist))
}
cvglasso<-function(X,K=5,lambda=NULL,lambda.min=0.1,nlambda=50)
{
p=ncol(X);
n=nrow(X);
Sn<-cov(X);
if (is.null(lambda)){lambda=exp(seq(log(lambda.min),0,length.out =nlambda))*max(abs(Sn))}
lambda<-sort(lambda,decreasing =TRUE)

fold<-split(1:n, rep(1:K, length =n));
CV<-matrix(0,nrow=K,ncol=nlambda);
for (k in 1:K){
  xtrain=X[-fold[[k]],];
  xtest=X[fold[[k]],];
  obj<-glassopath(cov(xtrain), rholist=lambda,penalize.diagonal=FALSE,trace=0)
  x1<-scale(xtest,center =TRUE,scale =FALSE)/sqrt(nrow(xtest)-1);
  for (i in 1:nlambda){
    hOme<-obj$wi[,,k];
    CV[k,i]=sum((x1%*%hOme)^2)/2-sum(diag(hOme));
  }
}
mcv<-apply(CV, 2, mean);
print(mcv)
lambda<-obj$lambda;
cvlambda<-lambda[which.min(mcv)];
fobj<-glassopath(Sn, rholist=cvlambda,penalize.diagonal=FALSE,trace=0);
return(fobj$wi[,,1])
}