# Introduction for R package EQUAL
We develop an Efficient admm algorithm via the QUAdratic Loss (EQUAL) for precision matrix estimation. This is my first R package and welcome any comments or suggestions.

# References 
Cheng Wang and Binyan Jiang. "[An efficient ADMM algorithm for high dimensional precision matrix estimation via penalized quadratic loss.](https://arxiv.org/abs/1811.04545)" arXiv preprint arXiv:1811.04545 (2018).

## Getting Started
These instructions will give you a toy example for implementing the package.

### Prerequisites
What things you need to install the software and how to install them

```
install.packages("MASS")
install.packages("Rcpp")
install.packages("Matrix")
install.packages("devtools")
```
### Install EQUAL

```
library("devtools")
devtools::install_github("cescwang85/EQUAL")
```

### Toy example 

```
rm(list = ls())
set.seed(123)
library('MASS')
library('Rcpp')
library('Matrix')
library('EQUAL')
n=200
p=100
Omega<-toeplitz(0.5^(1:p-1))
X=mvrnorm(n,rep(0,p),solve(Omega))
aa<-EQUAL(X)
bb<-EQUAL(X,type=FALSE);

obj1<-CVEQUAL(X);
obj2<-CVEQUAL(X,type=FALSE)
obj1$Omega[1:10,1:10]
obj2$Omega[1:10,1:10]
```
Overall, the computation complexity for each iteration of the algorithm is linear in both the sample size (n) and the number of parameters (p^2).  
