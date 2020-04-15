bidifac.plus <- function(X0,p.ind,n.ind,max.comb=20,num.comp=20,max.iter=500,conv.thresh=0.001,temp.iter=100){
# X0=X0; p.ind = p.ind; n.ind=n.ind; max.comb=20; num.comp=20; max.iter=500; conv.thresh=1; temp.iter=100
  n.source <- length(p.ind)
  n.type <- length(n.ind)
  X0.resid <- X0
  S <- list()
  pen <- c()
    n.ind.list <- list()
    p.ind.list <- list()
  for(i in 1:max.comb){
    S[[i]]=0
    pen[i]=0
    n.ind.list[[i]]<-c(1:dim(X0)[2])
    p.ind.list[[i]]<-c(1:dim(X0)[1])}
  obj.vec <- c(sum(X0^2))
  #obj.vec <- c()
  temp.fac <- svd(X0,nu=0,nv=0)$d[1]/(sum(sqrt(dim(X0))))-1
  obj.prev <- 10^10
  for(jj in 1:max.iter){
    print(jj)
    if(jj < temp.iter){
      lambda <- 1+(temp.iter-1-jj)/temp.iter*temp.fac
    }
    
   for(k in 1:max.comb){
  #    print(paste(k,'*'))
      X0.temp <- X0.resid+S[[k]]
      cur.n.ind <- n.ind.list[[k]]#c(1:dim(X0.temp)[2])
      cur.p.ind <- c()
      prev.p.ind <- p.ind.list[[k]]#c(1:dim(X0.temp)[1])
      prev.n.ind <- c()
      for(ii in 1:10){
        X0.temp.p <- X0.temp[,cur.n.ind]
        res.p <- bidifac.cycle(X0.temp=X0.temp,x.ind=p.ind,num.comp=num.comp, lambda=lambda)
        cur.p.ind <- sort(res.p$cur.x.ind)
        if(isTRUE(all.equal(prev.p.ind,cur.p.ind))&&isTRUE(all.equal(prev.n.ind,cur.n.ind))) break
          X0.temp.n <- t(X0.temp[cur.p.ind,])
        res.n<- bidifac.cycle(X0.temp=X0.temp.n,x.ind=n.ind,num.comp=num.comp,lambda=lambda)
        cur.n.ind <- sort(res.n$cur.x.ind)
        if(isTRUE(all.equal(prev.p.ind,cur.p.ind))&&isTRUE(all.equal(prev.n.ind,cur.n.ind))) break
        prev.p.ind <- cur.p.ind
        prev.n.ind <- cur.n.ind
      }
      CurPen <- res.n$CurPen
      if(!isTRUE(all.equal(p.ind.list[[k]],cur.p.ind))|!isTRUE(all.equal(n.ind.list[[k]],cur.n.ind))){
        Old.a <- svd(X0.temp[p.ind.list[[k]],n.ind.list[[k]]],nu=0,nv=0)$d
        New.a <- svd(X0.temp[cur.p.ind,cur.n.ind],nu=0,nv=0)$d
        OldPen <- sum(pmax(Old.a-lambda*(sqrt(length(p.ind.list[[k]]))+sqrt(length(n.ind.list[[k]]))),0)^2)
        NewPen <- sum(pmax(New.a-lambda*(sqrt(length(p.ind.list[[k]]))+sqrt(length(n.ind.list[[k]]))),0)^2)
        if(OldPen>NewPen){ 
          CurPen <- OldPen  
          cur.p.ind = p.ind.list[[k]]
          cur.n.ind = n.ind.list[[k]]
        }
      }
      Repeat=FALSE
      if(k>1){
      for(j in 1:(k-1)){
        if(isTRUE(all.equal(p.ind.list[[j]],cur.p.ind))&&isTRUE(all.equal(n.ind.list[[j]],cur.n.ind))){
          Repeat=TRUE
        }
      }}
      if(CurPen<=0.001|Repeat){
        pen = pen[-k]
        S[[k]] <- NULL
        pen[max.comb]=0
        S[[max.comb]]=0
        X0.resid=X0.temp
        break
      }
      p.ind.list[[k]] <- cur.p.ind
      n.ind.list[[k]] <- cur.n.ind
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
  num.modules=k
  if(S[[max.comb]]==0) num.modules=k-1
  Sums <- array(dim=c(num.modules,n.source,n.type))
  for(kk in 1:num.modules){for(j in 1:n.source){ for(i in 1:n.type){
    Sums[kk,j,i] = sum(S[[kk]][p.ind[[j]],n.ind[[i]]]^2)
  }}}
  p.ind.list=p.ind.list[1:num.modules]
  n.ind.list=n.ind.list[1:num.modules]
  return(list(S=S,p.ind.list=p.ind.list,n.ind.list=n.ind.list,Sums=Sums,obj.vec=obj.vec))
}

bidifac.cycle <- function(X0.temp,x.ind,num.comp,lambda){
  n <- dim(X0.temp)[2]
  n.x <- length(x.ind)
  temp.x <- c(1:n.x)
  X0.c <- tcrossprod(X0.temp)
  X0.rl <- list()
  for(i in 1:n.x) X0.rl[[i]] <- crossprod(X0.temp[x.ind[[i]],])
  cur.x.ind <- c()
  sse.pen.tot <- c()
  SSE.pen.temp <- c()
  cur.x.ind.list <- list()
  cur.n.ind <- c()
  cur_cross=0
  for(i in 1:n.x){
    SSE.pen.temp <- c()
    for(j in temp.x){
      temp.x.ind <- c(cur.x.ind,x.ind[[j]]) 
      pl <- length(temp.x.ind)
      if(pl<=n){
        a <- sqrt(eigs_sym(X0.c[temp.x.ind,temp.x.ind], num.comp, which = "LM",retvec=FALSE)$values)
      }  
      if(pl>n){
        a <- sqrt(eigs_sym(cur_cross+X0.rl[[j]], num.comp, which = "LM",retvec=FALSE)$values)
      }  
      SSE.pen.temp[j] <- sum(pmax(a-lambda*(sqrt(pl)+sqrt(n)),0)^2)
    }    
    indmin <- which.max(SSE.pen.temp)
    sse.pen.tot[i] <- SSE.pen.temp[indmin]
    cur.x.ind <- c(cur.x.ind,x.ind[[indmin]])
    cur.x.ind.list[[i]] <- cur.x.ind
    cur_cross <- cur_cross+X0.rl[[indmin]]
    temp.x <- temp.x[temp.x!=indmin]
  }
  Order = order(sse.pen.tot, decreasing=TRUE)
  i=1
 # for(i in 1:length(sse.pen.tot)){
    cur.x.ind <- cur.x.ind.list[[Order[i]]]
    CurPen <- sse.pen.tot[Order[i]]
 #   if(!(list(sort(cur.p.ind))%in%p.ind.list)) break
#    print('GGG!')
#  }
    return(list(cur.x.ind=cur.x.ind,CurPen=CurPen))
}

