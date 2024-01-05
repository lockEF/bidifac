#load functions
source('EV_BIDIFAC_Functions.R')

##EVB low-rank estimation for a single matrix: 
#Simulate rank 10 structure and noise
D=1000
N=100
R=10

E=array(data=rnorm(D*N),dim=c(D,N)) # noise
A=matrix(rnorm(D*R),ncol=R)%*%t(matrix(rnorm(N*R),ncol=R)) #signal
X=A+E
#Estimate signal using EVB method
A.est=eb.glob.o(X)
#check accuracy: RSE
sum((A.est-A)^2)/sum(A^2)

###EVB missing data imputation
E=array(data=rnorm(D*N),dim=c(D,N)) # noise
A=matrix(rnorm(D*R),ncol=R)%*%t(matrix(rnorm(N*R),ncol=R)) #signal
X=A+E
#set 20% of values to missing
X.na=X
missing=sample(length(X),length(X)*0.2) 
X.na[missing]=NA
#estimate matrix signal using EVB approach
A.est=eb.glob.imp(X.na)
#check accuracy of signal estimation for missing values
sum((A.est[missing]-A[missing])^2)/sum(A[missing]^2)

##bidimensional decomposition
#Generate matrix with to row sets (p=50 each) and two column sets (n=100 each)
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
#p.ind.list and n.ind.list give modules of shared/individual structure defined on row and column sets
#here, they enumeate all possibilities:
p.ind.list=list(c(1:100),c(1:100),c(1:100),p.ind[[1]],p.ind[[2]],p.ind[[1]],p.ind[[1]],p.ind[[2]],p.ind[[2]])
n.ind.list=list(c(1:200),n.ind[[1]],n.ind[[2]],c(1:200),c(1:200),n.ind[[1]],n.ind[[2]],n.ind[[1]],n.ind[[2]])
#estimate underlying structure and decomposition
res <- glob.bidi.o(X,p.ind,n.ind,p.ind.list,n.ind.list)
#res$S is a list of estimates for identified modules in the decomposition.  Sum to get overall structure:
A.est=Reduce('+',res$S)

#here first module, res$S[[1]], corresponds to global structure
sum(res$S[[1]]-J)^2/sum(J^2)
#res$S[[2]] corresponds to C
sum(res$S[[2]][p.ind.list[[2]],n.ind.list[[2]]]-C)^2/sum(C^2)
#res$S[[4]] corresponds to R
sum(res$S[[4]][p.ind.list[[4]],n.ind.list[[4]]]-R)^2/sum(R^2)
#res$S[[6]] corresponds to I
sum(res$S[[6]][p.ind.list[[6]],n.ind.list[[6]]]-I)^2/sum(I^2)

#####Missing data imputation

#Randomly set 10% of values to missing
X.miss=X
all.miss=sample(c(1:prod(dim(X))),prod(dim(X))/10)
X.miss[all.miss]=NA
#impute results, simultaneously estimatying all possible modules (given by p.ind.list and n.ind.list)
p.ind.list=list(c(1:100),c(1:50),c(1:50),c(51:100),c(51:100),c(51:100),c(1:50),c(1:100),c(1:100))
n.ind.list=list(c(1:200),c(1:100),c(101:200),c(1:100),c(101:200),c(1:200),c(1:200),c(1:100),c(101:200))
res.impute<-glob.bidi.imp(X.miss,p.ind=p.ind,n.ind=n.ind,p.ind.list,n.ind.list)
#get overall signal estimation
A.est=Reduce('+',res.impute$S)

#check relative error of imputed values
sum(A.est[all.miss]-X[all.miss])^2/sum(X[all.miss]^2)

