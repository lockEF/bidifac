set.seed(4122020)
#data
p.vec = c(50,50)
n.vec = c(100,100)
p = sum(p.vec)
n = sum(n.vec)
p.ind = list(c(1:50),c(51:100))
n.ind = list(c(1:100),c(101:200))
#initialize data X as a noisy matrix
X = matrix(rnorm(p*n),nrow=p,ncol=n)
#add rank-1 global structure 
J = rnorm(p)%*%t(rnorm(n))
X = X+J
#add rank-1 structure shared for row set 1 and both column sets
R = rnorm(p.vec[1])%*%t(rnorm(n))
X[p.ind[[1]],] = X[p.ind[[1]],]+R
#add rank-1 structure shared for both row sets and column set 1
C = rnorm(p)%*%t(rnorm(n.vec[1]))
X[,n.ind[[1]]] = X[,n.ind[[1]]]+C
#add rank-1 structure to row set 1 and column set 1
I = rnorm(p.vec[1])%*%t(rnorm(n.vec[1]))
X[p.ind[[1]],n.ind[[1]]] = X[p.ind[[1]],n.ind[[1]]]+I

#install package RSpectra
#library(RSpectra)

#run bidifac+
res <- bidifac.plus(X0=X,p.ind=p.ind,n.ind=n.ind)
#res$Sums[i,j,k] gives sum of squares for module i in submatrix (j,k)
#res$p.ind.list and res$n.ind.list are lists that give row and column sets, respectively, for each modules

#test accuracy 
#here first module, res$S[[1]], corresponds to global structure
sum(res$S[[1]]-J)^2/sum(J^2)
#res$S[[2]] corresponds to R
sum(res$S[[2]][res$p.ind.list[[2]],res$n.ind.list[[2]]]-R)^2/sum(R^2)
#res$S[[3]] corresponds to C
sum(res$S[[3]][p.ind.list[[3]],res$n.ind.list[[3]]]-C)^2/sum(C^2)
#res$S[[4]] corresponds to I
sum(res$S[[4]][p.ind.list[[4]],res$n.ind.list[[4]]]-I)^2/sum(I^2)
