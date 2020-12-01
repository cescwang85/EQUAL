#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
// [[Rcpp::export]]
mat soft(mat X,double lambda){
  vec X0=diagvec(X);
  mat Y=X-diagmat(X0);
  Y=(Y>=lambda)%(Y-lambda)+(Y<=(-lambda))%(Y+lambda);
  return Y+diagmat(X0);}

// [[Rcpp::export]]
Rcpp::List dtrace(mat X,vec lambda,double err=10^(-5),int maxIter=10^3,double rho=1){
  int p=X.n_cols;
  int nlambda=lambda.size();
  /*Centering int*/
  mat Sn=cov(X);
  mat U;
  vec eigd;
  eig_sym(eigd,U, Sn); 
  vec ld=exp(eigd);
  mat D=2/(log(ld*ld.t())+2*rho);
  Rcpp::List Omega_all(nlambda);
  vec niter=lambda;
  /*Intialization*/
  mat aZ=eye(p,p);
  mat aU=zeros(p,p);
  mat aX;
  mat L;
  mat Z1;
  double lam;
  double ee;
  for (int k=0;k<nlambda;++k) {
    lam=lambda(k);
    int i=0;
    while (((i<maxIter)&&(ee>err))|(i==0))
    { Z1=aZ;
      L=eye(p,p)+rho*(aZ-aU);
      aX=U*(D%(U.t()*L*U))*U.t();
      aZ=soft(aX+aU,lam/rho);
      aU=aU+aX-aZ;
      ee=mean(mean(abs(aZ-Z1)));
      i=i+1;
    }
    Omega_all(k)=sp_mat(aZ);
    niter(k)=i;
  }
  return Rcpp::List::create(Rcpp::Named("Omega") =Omega_all,
                            Rcpp::Named("lambda") =lambda,
                            Rcpp::Named("niter") =niter); }

// [[Rcpp::export]]
Rcpp::List gOmega(mat X,vec lambda,double err=10^(-5),int maxIter=10^3,double rho=1){
  int p=X.n_cols;
  int nlambda=lambda.size();
  /*Centering int*/
  mat Sn=cov(X);
  mat U;
  vec la;
  Rcpp::List Omega_all(nlambda);
  vec niter=lambda;
  /*Intialization*/
  mat aZ=eye(p,p);
  mat aU=zeros(p,p);
  mat aX;
  mat L;
  mat Z1;
  double lam;
  double ee;
  for (int k=0;k<nlambda;++k) {
    lam=lambda(k);
    int i=0;
    while (((i<maxIter)&&(ee>err))|(i==0))
    { Z1=aZ;
      eig_sym(la,U,rho*(aZ-aU)-Sn);
      aX=U*diagmat((la+sqrt(la%la+4*rho))/2/rho)*U.t();
      aZ=soft(aX+aU,lam/rho);
      aU=aU+aX-aZ;
      ee=mean(mean(abs(aZ-Z1)));
      i=i+1;
    }
    Omega_all(k)=sp_mat(aZ);
    niter(k)=i;
  }
  return Rcpp::List::create(Rcpp::Named("Omega") =Omega_all,
                            Rcpp::Named("lambda") =lambda,
                            Rcpp::Named("niter") =niter); }
// [[Rcpp::export]]
Rcpp::List sOmega(mat X,cube lambda,double err=10^(-5),int maxIter=10^3,double rho=1){
  int n=X.n_rows;
  int p=X.n_cols;
  int m=(p>n)*n+(p<=n)*p;
  int nlambda=lambda.n_slices;
  /*Centering int*/
  mat dX=(eye(n,n)-ones(n,n)/n)*X/sqrt(n-1);
  mat U;
  mat Uv;
  vec eigd;
  svd(U, eigd, Uv, dX.t()); 
  eigd=eigd%eigd;
  mat U1=U.cols(0,m-1)*diagmat(sqrt(eigd/(eigd+2*rho)));
  vec ld=exp(eigd);
  mat D=2*rho/(log(ld*ld.t())+2*rho)+1;
  Rcpp::List Omega_all(nlambda);
  vec niter=zeros<vec>(nlambda);
  /*Intialization*/
  mat aZ=eye(p,p);
  mat aU=zeros(p,p);
  mat aX;
  mat L;
  mat L1;
  mat L2;
  mat Z1;
  mat lam;
  double ee;
  for (int k=0;k<nlambda;++k) {
    lam=lambda.slice(k);
    lam.diag()=zeros<vec>(p);
    int i=0;
    while (((i<maxIter)&&(ee>err))|(i==0))
    { Z1=aZ;
      L=eye(p,p)/rho+aZ-aU;
      L=(L+L.t())/2;
      L1=L*U1;
      L2=L1*U1.t();
      aX=L-L2-L2.t()+U1*(D%(U1.t()*L1))*U1.t();
      aZ=(aX+aU>=lam/rho)%(aX+aU-lam/rho)+(aX+aU<=(-lam/rho))%(aX+aU+lam/rho);
      aU=aU+aX-aZ;
      ee=mean(mean(abs(aZ-Z1)));
      i=i+1;
    }
    Omega_all(k)=sp_mat(aZ);
    niter(k)=i;
  }
  return Rcpp::List::create(Rcpp::Named("Omega") =Omega_all,
                            Rcpp::Named("niter") =niter); }

