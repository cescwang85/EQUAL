#' CVEQUAL where the tuning parameter is chosen by cross-validation for each column
#' @param X data matrix of dimension n*p.
#' @param K the number of folds. Default is 5.
#' @param type Should the loss function be symmetric? Default is TRUE.
#' @param sdiag Should diagonal of inverse covariance be penalized? Default is FALSE.
#' @param lambda user supplied tuning parameter; Default is NULL and the program compute its own
#' sequence based on \code{nlambda}.
#' @param  nlambda the length of the tuning parameter sequence which is available when lambda is NULL. Default is 50.
#' @param  lambda.min smallest value for lambda, as a fraction of lambda.max which is available when lambda is NULL. 
#' Default is sqrt(log(p)/n).
#' @param err the precision used to stop the convergence. Default is 1e-5. 
#' Iterations stop when average absolute parameter change is less than \code{err}.
#' @param maxIter Maximum number of iterations. Default is 1000.
#' @param rho step parameter for the ADMM. Default is 1.
#' @return A list with components
#' \item{Omega}{the estimated p*p precision matrix.}
#' \item{cvlambda}{the chosen lambda by cross-validation.}
#' \item{lambda}{the used lambda list for cross-validation.}
#' \item{CVloss}{the empirical loss of cross-validation related to lambda.}
col_CVEQUAL<-function(X,K=5,type=TRUE,sdiag=FALSE,lambda=NULL,lambda.min=sqrt(log(ncol(X))/nrow(X)),nlambda=50,err=10^(-5),maxIter =1000,rho=1)
{
p=ncol(X);
n=nrow(X);
Sn<-cov(X);
A<-abs(Sn-diag(diag(Sn)));
if (is.null(lambda)){lambda=exp(seq(log(lambda.min),0,length.out =nlambda))*max(A)}
lambda<-sort(lambda,decreasing =TRUE)

fold<-split(1:n, rep(1:K, length =n));
CV<-array(0,dim=c(p,K,nlambda));
for (k in 1:K){
  xtrain=X[-fold[[k]],];
  xtest=X[fold[[k]],];
  S1<-cov(xtest);
  obj<-EQUAL(xtrain,sdiag=sdiag,type=type,lambda =lambda,err=10*err,maxIter =100)
  for (i in 1:nlambda){
    hOme<-as.matrix(obj$Omega[[i]]);
    CV[,k,i]=diag(hOme%*%S1%*%t(hOme))/2-diag(hOme)
  }
}
cumCV<-apply(CV,c(1,3),sum)
id<-apply(cumCV, 1, which.min);

hOmega<-EQUAL(X,lambda =lambda,sdiag=sdiag,type=type,err=err,maxIter=maxIter,rho=rho)$Omega

Omega_final=matrix(0,p,p)
for (j in 1:p)
{Omega_final[,j]=hOmega[[id[j]]][,j]}
Omega_final=(Omega_final+t(Omega_final))/2
return(list(Omega=Omega_final,cvlambda=lambda[id],lambda=lambda,CVloss=cumCV))
}

