#' Estimate the precision matrix with the given support from a list of sparse matrix
#' @param \code{X} data matrix of dimension \eqn{n\times p}.
#' @param \code{obj} a list of sparse \eqn{p \times p} matrices whose support is used to estimate the precision matrix. 
#' @param \code{rho} step parameter for the ADMM. Default is 1.
#' @param \code{err} the precision used to stop the convergence. Default is 1e-3. 
#' @return A  list of \eqn{p \times p} precision matrices whose support is consistent with the one of \code{obj}.

refit_list<-function(X,obj,rho=1,err=1e-3)
{ m=length(obj);
  re=obj;
  Sn=cov(X);
  p=ncol(X);
  n=nrow(X);
for (i in 1:m){  
  Omega=as.matrix(obj[[i]]);
  for (k in 1:p)
  {  id=(Omega[,k]!=0);
  if (sum(id)<n-1){
    ek=rep(0,p);
    ek[k]=1
    Omega[id,k]=solve(Sn[id,id],ek[id])
  }}
  re[[i]]=as((Omega+t(Omega))/2, "sparseMatrix")
}
return(re)  
}













