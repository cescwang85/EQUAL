#' Efficient admm algorithm via the QUAdratic Loss (EQUAL) for precision matrix estimation
#' @param X data matrix of dimension n*p.
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
#' \item{Omega}{a list of sparse p*p matrices corresponding to lambda.}
#' \item{lambda}{the used lambda for the solution path.}
#' \item{niter}{the number of iterations for each element of lambda.}
EQUAL<-function(X,type=TRUE,sdiag=FALSE,lambda=NULL,lambda.min=sqrt(log(ncol(X))/nrow(X)),nlambda=50,err=10^(-5),maxIter =1000,rho=1)
{p=ncol(X);
n=nrow(X);
Sn<-cov(X);
A<-abs(Sn-diag(diag(Sn)));
if (is.null(lambda)){lambda=exp(seq(log(lambda.min),0,length.out =nlambda))*max(A)}
lambda<-sort(lambda,decreasing =TRUE)

if (type){obj<-equal2(X,lambda =lambda,diag=as.numeric(sdiag),err=err,maxIter =maxIter,rho=rho)}
else {obj<-equal1(X,lambda =lambda,diag=as.numeric(sdiag),err=err,maxIter =maxIter,rho=rho)}
return(obj)
}  