#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


//' @title Element-wise Soft Thresholding Function for Matrix
//' @description Soft thresholding function for Matrix
//'
//' @param A Matrix
//' @param a scalar
//' @param diag  Should diagonal elements of the Matrix be thresholding? Default is FALSE.
//' @return Matrix after threholding
//' 
//' @export
//'
// [[Rcpp::export]]
arma::mat soft(arma::mat A,double a, int diag=0){
  a=std::abs(a);
  arma::mat C=(A>=a)%(A-a)+(A<=(-a))%(A+a);
if (diag==0){
  arma::mat B=A-arma::diagmat(arma::diagvec(A));
  B=(B>=a)%(B-a)+(B<=(-a))%(B+a);
  B=B+arma::diagmat(arma::diagvec(A));
  C=B;}
  return C;}


//' @title Non-symmetric version of EQUAL
//' @description ADMM algorithm for high precision matrix estimation using asymmetic loss
//' @param X data matrix of dimension n*p.
//' @param lambda user supplied tuning parameter; 
//' @param err the precision used to stop the convergence. Default is 1e-5. 
//' Iterations stop when average absolute parameter change is less than \code{err}.
//' @param maxIter Maximum number of iterations. Default is 1000.
//' @param rho step parameter for the ADMM. Default is 1.
//' @param diag  Should diagonal elements of the Matrix be thresholding? Default is FALSE.
//' @return A list with components
//' \item{Omega}{a list of sparse p*p matrices corresponding to lambda.}
//' \item{lambda}{the used lambda for the solution path.}
//' \item{niter}{the number of iterations for each element of lambda.}
//' 
//' @export
//' 
// [[Rcpp::export]]
Rcpp::List equal1(arma::mat X,arma::vec lambda,double err=10^(-5),int maxIter=1000,double rho=1, int diag=0){
  int n=X.n_rows;
  int p=X.n_cols;
  int m=(p>n)*n+(p<=n)*p;
  int nlambda=lambda.size();
  /*Centering int*/
  arma::mat dX=(arma::eye(n,n)-arma::ones(n,n)/n)*X/sqrt(n-1);
  arma::mat U;
  arma::mat Uv;
  arma::vec eigd;
  svd(U, eigd, Uv, dX.t()); 
  U=U.cols(0,m-1);
  eigd=eigd%eigd;
  arma::mat D=U*arma::diagmat(eigd/(eigd+rho));
  Rcpp::List Omega_all(nlambda);
  arma::vec niter=lambda;
  /*Intialization*/
  arma::mat aZ=arma::eye(p,p);
  arma::mat aU=arma::zeros(p,p);
  arma::mat aX;
  arma::mat L;
  arma::mat Z1;
  double lam;
  double ee=1;
  for (int k=0;k<nlambda;++k) {
    lam=lambda(k);
    int i=0;
    while (((i<maxIter)&&(ee>err))||(i==0))
    { Z1=aZ;
      L=arma::eye(p,p)/rho+aZ-aU;
      aX=L-D*(U.t()*L);
      aZ=soft(aX+aU,lam/rho,diag);
      aU=aU+aX-aZ;
      ee=mean(mean(abs(aZ-Z1)));
      i=i+1;
    }
    Omega_all(k)=arma::sp_mat((abs(aZ)<abs(aZ.t()))%aZ+(abs(aZ)>=abs(aZ.t()))%aZ.t());
    niter(k)=i;
  }
  return Rcpp::List::create(Rcpp::Named("Omega") =Omega_all,
                            Rcpp::Named("lambda") =lambda,
                            Rcpp::Named("niter") =niter); }




//' @title Symmetric version of EQUAL
//' @description ADMM algorithm for high precision matrix estimation using symmetic loss
//' @param X data matrix of dimension n*p.
//' @param lambda user supplied tuning parameter; 
//' @param err the precision used to stop the convergence. Default is 1e-5. 
//' Iterations stop when average absolute parameter change is less than \code{err}.
//' @param maxIter Maximum number of iterations. Default is 1000.
//' @param rho step parameter for the ADMM. Default is 1.
//' @param diag  Should diagonal elements of the Matrix be thresholding? Default is FALSE.
//' @return A list with components
//' \item{Omega}{a list of sparse p*p matrices corresponding to lambda.}
//' \item{lambda}{the used lambda for the solution path.}
//' \item{niter}{the number of iterations for each element of lambda.}
//' 
//' @export
//' 
// [[Rcpp::export]]
Rcpp::List equal2(arma::mat X,arma::vec lambda,double err=10^(-5),int maxIter=1000,double rho=1,int diag=0){
  int n=X.n_rows;
  int p=X.n_cols;
  int m=(p>n)*n+(p<=n)*p;
  int nlambda=lambda.size();
  /*Centering int*/
  arma::mat dX=(arma::eye(n,n)-arma::ones(n,n)/n)*X/sqrt(n-1);
  arma::mat U;
  arma::mat Uv;
  arma::vec eigd;
  svd(U, eigd, Uv, dX.t()); 
  eigd=eigd%eigd;
  arma::mat U1=U.cols(0,m-1)*arma::diagmat(sqrt(eigd/(eigd+2*rho)));
 arma::vec ld=exp(eigd);
  arma::mat D=2*rho/(log(ld*ld.t())+2*rho)+1;
 Rcpp::List Omega_all(nlambda);
 arma::vec niter=lambda;
 /*Intialization*/
 arma::mat aZ=arma::eye(p,p);
 arma::mat aU=arma::zeros(p,p);
 arma::mat aX;
 arma::mat L;
 arma::mat L1;
 arma::mat L2;
 arma::mat Z1;
 double lam;
 double ee=1;
 for (int k=0;k<nlambda;++k) {
   lam=lambda(k);
   int i=0;
   while (((i<maxIter)&&(ee>err))||(i==0))
   { Z1=aZ;
     L=arma::eye(p,p)/rho+aZ-aU;
     L=(L+L.t())/2;
     L1=L*U1;
     L2=L1*U1.t();
     aX=L-L2-L2.t()+U1*(D%(U1.t()*L1))*U1.t();
     aZ=soft(aX+aU,lam/rho,diag);
     aU=aU+aX-aZ;
     ee=mean(mean(abs(aZ-Z1)));
     i=i+1;
   }
   Omega_all(k)=arma::sp_mat(aZ);
   niter(k)=i;
 }
 return Rcpp::List::create(Rcpp::Named("Omega") =Omega_all,
                           Rcpp::Named("lambda") =lambda,
                           Rcpp::Named("niter") =niter); }
