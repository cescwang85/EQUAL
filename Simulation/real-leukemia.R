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
leukemia_small <- read.csv("http://web.stanford.edu/~hastie/CASI_files/DATA/leukemia_small.csv")
id<-'ALL'
for (k in 1:46)
{id<-c(id,paste('ALL.',k,sep=''))}
X<-t(leukemia_small[,id]);
dim(X)
n=nrow(X);
p=ncol(X);
X<-scale(X,center =TRUE,scale =TRUE)
nlambda=100
Sn<-cov(X);
An<-abs(Sn-diag(diag(Sn)));
lambda=exp(seq(log(0.5),0,length.out =nlambda))*max(An)
ttime(obj1<-EQUAL(X,lambda=lambda))
saveRDS(obj1,file='ALL.RDS')


ttime(obj3<-glassopath(Sn,rholist=lambda,trace =0))
saveRDS(obj3,file='gALL.RDS')

#obj1<-readRDS('ALL.RDS')

Re<-rep(0,nlambda)
for ( k in 1:nlambda)
{
  
  Re[k]=nnzero(obj1$Omega[[k]])/p;
}
pdf("fig/fig5.pdf")
plot(obj1$lambda,Re,type='o',xlab ='lambda',ylab='Sparsity level')
dev.off()

A<-as.matrix(obj1$Omega[[63]])
obj1$lambda[[63]]
nnzero(A)

lab1<-apply(A, 1, nnzero)
aa<-sort(lab1,decreasing =TRUE)[round(p*0.02)]
id1<-(lab1>=aa)
A=A[id1,id1]
dim(A)
g1<- graph_from_adjacency_matrix(A!=0,mode='undirected',diag=FALSE)
layout.grid = layout.fruchterman.reingold(g1)
pdf("fig/fig6.pdf")
plot(g1, edge.color='gray50',vertex.color="gray50", vertex.size=2, vertex.label=NA,layout=layout.circle)
dev.off()


