# Introduction for the R package EQUAL
We develop an R package EQUAL for the paper
"Efficient admm algorithm via the QUAdratic Loss (EQUAL) for precision matrix estimation" [arxiv](https://arxiv.org/abs/1811.04545)  

## Welcome Any Comments or Suggestions
This is my first R package and welcome any comments or suggesiongs.

## Getting Started

These instructions will give you a toy Example for implement the package.

### Prerequisites

What things you need to install the software and how to install them

```
require("MASS")
require("Rcpp")
require("Matrix")
```
### Toy Example for EQUAL

```
rm(list = ls())
set.seed(520)
require("devtools")
devtools::install_github("cescwang85/EQUAL")
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
