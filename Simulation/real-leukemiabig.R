rm(list = ls())
set.seed(1985)
setwd('~/Dropbox/Share/Cheng&Yan/matinv/Simulation/')
library('Rcpp')
library('glasso')
library('Matrix')
library('EQUAL')
library('igraph')
ttime<-function(obj)
{t1<-proc.time();
obj
t2<-as.vector(proc.time()-t1);
return(sum(t2[1:2]))}
leukemia_big <- read.csv("http://web.stanford.edu/~hastie/CASI_files/DATA/leukemia_big.csv")
id<-'ALL'
for (k in 1:46)
{id<-c(id,paste('ALL.',k,sep=''))}
X<-t(leukemia_big[,id]);
dim(X)
n=nrow(X);
p=ncol(X);
X<-scale(X,center =TRUE,scale =TRUE)
nlambda=50
Sn<-cov(X);
An<-abs(Sn-diag(diag(Sn)));
lambda=exp(seq(log(0.5),0,length.out =nlambda))*max(An)
ttime(obj1<-EQUAL(X,lambda=lambda))
saveRDS(obj1,file='bigALL.RDS')


ttime(obj3<-glassopath(Sn,rholist=lambda,trace =0))
saveRDS(obj3,file='gbigALL.RDS')

obj1<-readRDS('ALL.RDS')

Re<-rep(0,nlambda)
for ( k in 1:nlambda)
{
  
  Re[k]=nnzero(obj1$Omega[[k]])/p;
}
pdf("fig/fig5.pdf")
plot(obj1$lambda,Re,type='o',xlab ='lambda',ylab='Sparsity level')
dev.off()

A<-as.matrix(obj1$Omega[[36]])
obj1$lambda[[36]]
nnzero(A)

lab1<-apply(A, 1, nnzero)
aa<-sort(lab1,decreasing =TRUE)[round(p*0.01)]
id1<-(lab1>=aa)
A=A[id1,id1]
dim(A)
g1<- graph_from_adjacency_matrix(A!=0,mode='undirected',diag=FALSE)
layout.grid = layout.fruchterman.reingold(g1)
pdf("fig/fig6.pdf")
#plot(g1, edge.color='gray50',vertex.color="gray50", vertex.size=2, vertex.label=rownames(leukemia_big)[id1],layout=layout.circle)
plot(g1, edge.color='gray50',vertex.color="gray50", vertex.size=2, vertex.label='',layout=layout.circle)
dev.off()


