#' Estimate the precision matrix with the given support
#' @param \code{X} data matrix of dimension \eqn{n\times p}.
#' @param \code{obj} a \eqn{p \times p} matrix whose support is used to estimate the precision matrix. 
#' @param \code{rho} step parameter for the ADMM. Default is 1.
#' @param \code{err} the precision used to stop the convergence. Default is 1e-3. 
#' @return A \eqn{p \times p} precision matrix whose support is consistent with the one of \code{obj}.

refit_mat<-function(X,obj,rho=1,err=1e-3)
{ 
aa=as(obj,"sparseMatrix")
p=ncol(X);
H=solve(cov(X)+rho*diag(p));
idx<-1+aa@i
idy<-rep(1:p,diff(aa@p))
aZ=diag(p)
aU=matrix(0,p,p)
##iterations
for (k in 1:100){
  aX=H%*%(diag(p)+rho*(aZ-aU));
  aZ=matrix(0,p,p)
  aZ[idx+(idy-1)*p]=(aX+aU)[idx+(idy-1)*p];
  aU=aU+aX-aZ
  ee=norm(aX-aZ,type='F')
  if (ee<err) break;
}
return(aZ)  
}

