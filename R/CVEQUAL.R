#' CVEQUAL where the tuning parameter is chosen by cross-validation
#' @param X data matrix of dimension n*p.
#' @param K the number of folds. Default is 5.
#' @param type Should the loss function be symmetric? Default is TRUE.
#' @param sdiag Should diagonal of inverse covariance be penalized? Default is FALSE.
#' @param lambda user supplied tuning parameter; Default is NULL and the program compute its own
#' sequence based on \code{nlambda}.
#' @param  nlambda the length of the tuning parameter sequence which is available when lambda is NULL. Default is 50.
#' @param  lambda.min.ratio smallest value for lambda, as a fraction of lambda.max which is available when lambda is NULL. 
#' Default is sqrt(log(p)/n).
#' @param err the precision used to stop the convergence. Default is 1e-5.
#' @param maxIter Maximum number of iterations. Default is 1000.
#' @param rho step parameter for the ADMM. Default is 1.
#' @return A list with components
#' \item{Omega}{the estimated p*p precision matrix.}
#' \item{cvlambda}{the chosen lambda by cross-validation.}
#' \item{lambda}{the used lambda list for cross-validation.}
#' \item{CVloss}{the empirical loss of cross-validation related to lambda.}
CVEQUAL<-function(X,K=5,type=TRUE,sdiag=FALSE,lambda=NULL,lambda.min=sqrt(log(ncol(X))/nrow(X)),nlambda=50,err=10^(-5),maxIter =1000)
{p=ncol(X);
n=nrow(X);
Sn<-cov(X);
A<-abs(Sn-diag(diag(Sn)));
if (is.null(lambda)){lambda=exp(seq(log(lambda.min),0,length.out =nlambda))*max(A)}
lambda<-sort(lambda,decreasing =TRUE)

fold<-split(1:n, rep(1:K, length =n));
CV<-matrix(0,nrow=K,ncol=nlambda);
for (k in 1:K){
  xtrain=X[-fold[[k]],];
  xtest=X[fold[[k]],];
  x1<-scale(xtest,center =TRUE,scale =FALSE)/sqrt(nrow(xtest)-1);
  obj<-EQUAL(xtrain,sdiag=sdiag,type=type,lambda =lambda,err=10*err,maxIter =100)
  for (i in 1:nlambda){
    hOme<-as.matrix(obj$Omega[[i]]);
    CV[k,i]=sum((x1%*%hOme)^2)/2-sum(diag(hOme));
  }
}
mcv<-apply(CV, 2, mean);
lambda<-obj$lambda;
cvlambda<-lambda[which.min(mcv)];
hOmega<-EQUAL(X,lambda =cvlambda,sdiag=sdiag,type=type,err=err,maxIter=maxIter)
return(list(Omega=as.matrix(hOmega$Omega[[1]]),cvlambda=cvlambda,lambda=lambda,CVloss=mcv))
}

