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

prostmat <- read.csv("http://web.stanford.edu/~hastie/CASI_files/DATA/prostmat.csv")
#prostmat<-read.csv('prostmat.csv')
X<-t(prostmat[,1:50]);
Y<-t(prostmat[,51:102])
X<-scale(X,center =TRUE,scale =TRUE)
Y<-scale(Y,center =TRUE,scale =TRUE)
nlambda=50
S1<-cov(X);
A1<-abs(S1-diag(diag(S1)));
lambda1=exp(seq(log(0.5),0,length.out =nlambda))*max(A1)
lambda1<-sort(lambda1,decreasing =TRUE)


S2<-cov(Y);
A2<-abs(S2-diag(diag(S2)));
lambda2=exp(seq(log(0.5),0,length.out =nlambda))*max(A2)
lambda2<-sort(lambda2,decreasing =TRUE)

#ttime(a1<-EQUAL(X,lambda=lambda1,type =FALSE))
#saveRDS(a1,file='b1.RDS')
#ttime(a2<-EQUAL(Y,lambda=lambda2,type=FALSE))
#saveRDS(a2,file='b2.RDS')
ttime(a3<-EQUAL(X,lambda=lambda1))
saveRDS(a3,file='b3.RDS')
ttime(a4<-EQUAL(Y,lambda=lambda2))
saveRDS(a4,file='b4.RDS')

#ttime(a5<-glassopath(S1,rholist=lambda1,trace =0))
#saveRDS(a5,file='b5.RDS')
#ttime(a6<-glassopath(S2,rholist=lambda2,trace =0))
#saveRDS(a6,file='b6.RDS')


obj1<-readRDS('b3.Rds')
obj2<-readRDS('b4.Rds')
Re<-matrix(0,nlambda,2)
for ( k in 1:nlambda)
{
  
  Re[k,1]=nnzero(obj1$Omega[[k]])/6033;
  Re[k,2]=nnzero(obj2$Omega[[k]])/6033;
}
pdf("fig/fig1.pdf")
plot(lambda1,Re[,1],type='o',xlab ='lambda',ylab='Sparsity level')
dev.off()
pdf("fig/fig2.pdf")
plot(lambda2,Re[,2],type='o',xlab ='lambda',ylab='Sparsity level')
dev.off()
A<-obj1$Omega[[21]];
B<-obj2$Omega[[21]]
lab1<-apply(A, 1, nnzero)
lab2<-apply(B, 1, nnzero)
id<-(lab1>=2)+(lab2>=2)
id=(id>=2)
sum(id)
A=A[id,id]
B=B[id,id]
g1<- graph_from_adjacency_matrix(A!=0,mode='undirected',diag=FALSE)
g2<- graph_from_adjacency_matrix(B!=0,mode='undirected',diag=FALSE)
layout.grid = layout.fruchterman.reingold(g2)
pdf("fig/fig3.pdf")
plot(g1, layout=layout.grid, edge.color='gray50',vertex.color="gray50", vertex.size=2, vertex.label=NA)
dev.off()
pdf("fig/fig4.pdf")
plot(g2, layout=layout.grid, edge.color='gray50',vertex.color="gray50", vertex.size=2, vertex.label=NA)
dev.off()


