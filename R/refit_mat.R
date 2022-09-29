#' Estimate the precision matrix with the given support
#' @param \code{X} data matrix of dimension \eqn{n\times p}.
#' @param \code{obj} a \eqn{p \times p} matrix whose support is used to estimate the precision matrix. 
#' @param \code{rho} step parameter for the ADMM. Default is 1.
#' @param \code{err} the precision used to stop the convergence. Default is 1e-3. 
#' @return A \eqn{p \times p} precision matrix whose support is consistent with the one of \code{obj}.

refit_mat<-function(X,obj,rho=1,err=1e-3)
{ 
Omega=as.matrix(obj);
Sn=cov(X);
p=ncol(X);
n=nrow(X);
for (k in 1:p)
{  id=(Omega[,k]!=0);
  if (sum(id)<n-1){
  ek=rep(0,p);
  ek[k]=1
  Omega[id,k]=solve(Sn[id,id],ek[id])
  }
}
Omega=(Omega+t(Omega))/2
return(as(Omega, "sparseMatrix"))  
}

