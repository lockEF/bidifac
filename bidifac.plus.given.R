bidifac.plus.given <- function(X0,p.ind,n.ind,p.ind.list,n.ind.list, S = list(), pen=c(), given.inits = FALSE,max.iter=500,conv.thresh=0.001){

  n.source <- length(p.ind)
  n.type <- length(n.ind)
if(given.inits==FALSE){
  for(i in 1:max.comb){
    S[[i]]=0
    pen[i]=0}
}
obj.vec <- c(sum(X0^2))
max.comb=length(p.ind.list)
for(jj in 1:max.iter){
  print(jj)
  for(k in 1:max.comb){
    print(paste('*',k))
    X0.temp <- X0.resid+S[[k]]
    cur.p.ind <- p.ind.list[[k]] 
    cur.n.ind <- n.ind.list[[k]] 
    a <- svd(X0.temp[cur.p.ind,cur.n.ind],nu=0,nv=0)$d
    nc <- sum(a>(lambda*(sqrt(length(cur.p.ind))+sqrt(length(cur.n.ind)))))
    SVD <- svd(X0.temp[cur.p.ind,cur.n.ind], nu=nc,nv=nc)
    s.vals <- pmax(SVD$d[1:nc]-lambda*(sqrt(length(cur.p.ind))+sqrt(length(cur.n.ind))),0)
    Diag <- diag(s.vals,nrow=nc,ncol=nc)
    Est <- SVD$u%*%Diag%*%t(SVD$v)
    S[[k]] <- array(rep(0,prod(dim(X0))),dim=dim(X0))
    S[[k]][cur.p.ind,cur.n.ind] <- Est
    X0.resid[cur.p.ind,cur.n.ind] <- X0.temp[cur.p.ind,cur.n.ind]-Est
    X0.resid = X0-Reduce('+',S)
    pen[k] <- lambda*(sqrt(length(cur.n.ind))+sqrt(length(cur.p.ind)))*(sum(s.vals))
    obj.vec <- c(obj.vec,sum(X0.resid^2)+2*sum(pen))
  }
  obj.cur <- sum(X0.resid^2)+2*sum(pen)
  if(abs(obj.cur-obj.prev)<conv.thresh) break
  obj.prev <- sum(X0.resid^2)+2*sum(pen)
}

Sums <- array(dim=c(max.comb,n.source,n.type))
for(kk in 1:max.comb){for(j in 1:n.source){ for(i in 1:n.type){
  Sums[kk,j,i] = sum(S[[kk]][p.ind[[j]],n.ind[[i]]]^2)
}}}

return(list(S=S,Sums=Sums,obj.vec=obj.vec))
}





