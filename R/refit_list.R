#' Estimate the precision matrix with the given support from a list of sparse matrix
#' @param \code{X} data matrix of dimension \eqn{n\times p}.
#' @param \code{obj} a list of sparse \eqn{p \times p} matrices whose support is used to estimate the precision matrix. 
#' @param \code{rho} step parameter for the ADMM. Default is 1.
#' @param \code{err} the precision used to stop the convergence. Default is 1e-3. 
#' @return A  list of \eqn{p \times p} precision matrices whose support is consistent with the one of \code{obj}.

refit_list<-function(X,obj,rho=1,err=1e-3)
{
  p=ncol(X);
  re<-obj
  m=length(obj);
  H=solve(cov(X)+rho*diag(p));
for (i in 1:m){  
  aa=obj[[i]]
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
  re[[i]]=Matrix(aZ)}
return(re)  
}













