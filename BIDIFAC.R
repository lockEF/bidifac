library(gplots)
library(denoiseR)
library(Matrix)

frob=function(X){ sum(X^2,na.rm=T) }

diag2=function(x){
  if(length(x)==1) return(as.matrix(x))
  else if(length(x)>1) return(diag(x))
}

sample2=function(x){
  if (length(x)==1) return(x)
  else if (length(x)>1) return(sample(x))
}

sigma.rmt=function(X){ estim_sigma(X,method="MAD") }

fill.matrix=function(X){
  na.ind.original=na.ind=which(is.na(X),arr.ind = T)
  bool=T
  while (bool){
    impute.X=rep(NA, nrow(na.ind))
    for (j in 1:nrow(na.ind)){
      colmean=mean(X[,na.ind[j,2]], na.rm=T)
      rowmean=mean(X[na.ind[j,1],], na.rm=T)
      impute.X[j]=mean(c(colmean,rowmean), na.rm=T)
    }
    X[na.ind]=impute.X
    
    na.ind=which(is.na(X),arr.ind = T)
    if (length(na.ind)==0) bool=F
  }

  return(list(X.fill=X, na.ind=na.ind.original))
}

softSVD=function(X, lambda){
  svdX=svd(X)
  nuc=pmax(svdX$d-lambda,0)
  out=tcrossprod(svdX$u, tcrossprod( svdX$v,diag(nuc) ))
  return(list(out=out, nuc=sum(nuc)))
}

data.rearrange=function(data,rmt=F,sigma=NULL){
  out=NULL
  p=nrow(data)
  q=ncol(data)
  
  m.vec=rep(NA,p)
  n.vec=do.call(c, lapply(data[1,], ncol))

  if (is.null(sigma)) sigma=matrix(1,p,q)

  for (i in 1:p){
    dimm=do.call(cbind, lapply(data[i,],dim))
    m1=unique(dimm[1,])
    if (length(m1)==1 ){m.vec[i]=m1 } 
    else{ stop("the number of rows do not match.") }
    if (!all(dimm[2,], n.vec)){ stop("the number of columns do not match")}
    
    for (j in 1:q){
      if (rmt) sigma[i,j]=sigma.rmt(data[[i,j]])
      data[[i,j]]=data[[i,j]]/sigma[i,j]
    }
    
    out=rbind(out,do.call(cbind,data[i,]))
  }
  
  return(list(out=out, nrows=m.vec, ncols=n.vec, sigma.mat=sigma))
}

BIDsim=function(m.vec, n.vec, 
                  rkG=1, rkC=1, rkR=1, rkI=1, SNR=1){
  p=length(m.vec); q=length(n.vec)
  
  m.vec.end=cumsum(m.vec)
  m.vec.start=c(1,m.vec.end[-p]+1)
  n.vec.end=cumsum(n.vec)
  n.vec.start=c(1,n.vec.end[-q]+1)
  
  X=S=G=R=C=I=matrix(list(), p,q)
  
  if (length(rkR)==1){rkR=rep(rkR, p)}
  if (length(rkC)==1){rkC=rep(rkC, q)}
  if (length(rkI)==1){rkI=matrix(rkI, p,q)}
  vec.rk=c(rkG, rkR, rkC, rkI)
  sum.rk=sum(vec.rk)
  
  rk.ind.end=cumsum(vec.rk)
  rk.ind.start=c(1,rk.ind.end[-((p+1)*(q+1))]+1)
  
  sigma.mat=1/sqrt(tcrossprod(m.vec,n.vec))/SNR
  
  m=sum(m.vec); n=sum(n.vec)
  svdX=svd(matrix(rnorm(m*n),m,n))
  min.mn=min(m,n)
  
  if (sum.rk>min.mn) stop("the rank of the signal is too high")
  
  ind=sample(min.mn, sum.rk)
  svdX.u=svdX$u[,ind]
  svdX.d=svdX$d[1:sum.rk]/sqrt(frob(svdX$d[1:sum.rk]))
  svdX.v=svdX$v[,ind]
  
  #generate G
  ind.G=rk.ind.start[1]:rk.ind.end[1]
  for (i in 1:p){
    U.g=svdX.u[m.vec.start[i]:m.vec.end[i],ind.G]
    for (j in 1:q){
      D.g=diag2(sample2(svdX$d[ind.G]))
      V.g=svdX.v[n.vec.start[j]:n.vec.end[j], ind.G]
      G[[i,j]]=tcrossprod(U.g, tcrossprod(V.g,D.g))
    }
  }
  rk.ind.start=rk.ind.start[-1]
  rk.ind.end=rk.ind.end[-1]
  
  #generate R
  for (i in 1:p){
    ind.R=rk.ind.start[i]:rk.ind.end[i]
    U.r=svdX.u[m.vec.start[i]:m.vec.end[i],ind.R]
    for (j in 1:q){
      D.r=diag2(sample2(svdX$d[ind.R]))
      V.r=svdX.v[n.vec.start[j]:n.vec.end[j], ind.R]
      R[[i,j]]=tcrossprod(U.r, tcrossprod(V.r,D.r))
    }
  }
  rk.ind.start=rk.ind.start[-(1:p)]
  rk.ind.end=rk.ind.end[-(1:p)]
  
  #generate C
  for (j in 1:q){
    ind.C=rk.ind.start[j]:rk.ind.end[j]
    V.c=svdX.v[n.vec.start[j]:n.vec.end[j], ind.C]
    for (i in 1:p){
      D.c=diag2(sample2(svdX$d[ind.C]))
      U.c=svdX.u[m.vec.start[i]:m.vec.end[i],ind.C]
      C[[i,j]]=tcrossprod(U.c, tcrossprod(V.c,D.c))
    }
  }
  rk.ind.start=matrix(rk.ind.start[-(1:q)],p)
  rk.ind.end=matrix(rk.ind.end[-(1:q)],p)
    
  #generate I
  for (i in 1:p){
    for (j in 1:q){
      ind.I=rk.ind.start[i,j]:rk.ind.end[i,j]
      U.i=svdX.u[m.vec.start[i]:m.vec.end[i],ind.I]
      D.i=diag2(sample2(svdX$d[ind.I]))
      V.i=svdX.v[n.vec.start[j]:n.vec.end[j], ind.I]
      I[[i,j]]=tcrossprod(U.i, tcrossprod(V.i,D.i))
    }
  }

  for (i in 1:p){
    for (j in 1:q){
      d=svd(G[[i,j]]+R[[i,j]]+C[[i,j]]+I[[i,j]])
      G[[i,j]]=G[[i,j]]/sqrt(frob(d$d))
      R[[i,j]]=R[[i,j]]/sqrt(frob(d$d))
      C[[i,j]]=C[[i,j]]/sqrt(frob(d$d))
      I[[i,j]]=I[[i,j]]/sqrt(frob(d$d))
      
      d$d=d$d/sqrt(frob(d$d))
      S[[i,j]]=d$u%*%diag(d$d)%*%t(d$v)
      X[[i,j]]=S[[i,j]]+replicate(n.vec[j],rnorm(m.vec[i],0,sigma.mat[i,j] ))
    }
  }
  
  return(list(X=X,S=S,G=G,R=R,C=C,I=I,sigma=sigma.mat))
}

BIDIFAC=function(data,rmt=T, sigma=NULL,
                 start=NULL, out=FALSE,
                 eps=1e-3, max.iter=1000, pbar=TRUE, seed=NULL, ...){
  if (!is.null(seed)){set.seed(seed)}
  if (!rmt & class(sigma)!="matrix") stop("sigma must be a matrix.")
  
  fit=data.rearrange(data, rmt, sigma)
  sigma.mat=fit$sigma.mat
  X00=fit$out
  
  mvec=fit$nrows; nvec=fit$ncols
  p=length(mvec); q=length(nvec)
  rm(fit)
  
  start.ind.m=c(1, cumsum(mvec)[1:(p-1)]+1)
  end.ind.m=cumsum(mvec)
  
  start.ind.n=c(1, cumsum(nvec)[1:(q-1)]+1)
  end.ind.n=cumsum(nvec)
  
  lambda.G=sqrt(sum(mvec))+sqrt(sum(nvec))
  lambda.R=sqrt(mvec)+sqrt(sum(nvec))
  lambda.C=sqrt(sum(mvec))+sqrt(nvec)
  lambda.I=tcrossprod(sqrt(mvec), rep(1, length(nvec)))+
    tcrossprod(rep(1, length(mvec)),sqrt(nvec))

  if (!is.null(start)){
    G00=start[[1]]; R00=start[[2]]
    C00=start[[3]]; I00=start[[4]]
  } else {
    G00=replicate(sum(nvec),rnorm(sum(mvec)))
    R00=replicate(sum(nvec),rnorm(sum(mvec)))
    C00=replicate(sum(nvec),rnorm(sum(mvec)))
    I00=replicate(sum(nvec),rnorm(sum(mvec)))
  }
  
  G00.nuc=NA; R00.nuc=rep(NA, p)
  C00.nuc=rep(NA, q); I00.nuc=matrix(NA,p,q)

  bool=TRUE
  count=1;crit0=0
  if (pbar) pb = txtProgressBar(min = 0, max=max.iter, initial=0, char="-", style = 3)
  while (bool){
    if (pbar){  setTxtProgressBar(pb, count)  }
    crit0.old = crit0
    
    #Update G
    fit1=softSVD(X00-R00-C00-I00,lambda.G)
    G00=fit1$out; G00.nuc=fit1$nuc
    
    #update R
    for (i in 1:p){
      ind=start.ind.m[i]:end.ind.m[i]
      fit1=softSVD(X00[ind,]-G00[ind,]-C00[ind,]-I00[ind,], lambda.R[i])
      R00[ind,]=fit1$out; R00.nuc[i]=fit1$nuc
    }
    
    for (j in 1:q){
      ind=start.ind.n[j]:end.ind.n[j]
      fit1=softSVD(X00[,ind]-G00[,ind]-R00[,ind]-I00[,ind], lambda.C[j])
      C00[,ind]=fit1$out; C00.nuc[j]=fit1$nuc
    }
    
    for (i in 1:p){
      for (j in 1:q){
        ind1= start.ind.m[i]:end.ind.m[i]
        ind2=start.ind.n[j]:end.ind.n[j]
        fit1=softSVD(X00[ind1,ind2]-G00[ind1,ind2]-R00[ind1,ind2]-C00[ind1,ind2], lambda.I[i,j])
        I00[ind1,ind2]=fit1$out; I00.nuc[i,j]=fit1$nuc
      }
    }
    
    crit0 = frob(X00-G00-R00-C00-I00)+
      2*lambda.G*G00.nuc+2*sum(lambda.R*R00.nuc)+
      2*sum(lambda.C*C00.nuc)+2*sum(lambda.I*I00.nuc)

    if (abs(crit0.old-crit0)<eps){ bool=FALSE }
    else if (count==max.iter){ bool=FALSE}
    else{ count = count+1 }
  }
  
  S00.mat=G00.mat=R00.mat=C00.mat=I00.mat=data
  for (i in 1:p){
    ind1= start.ind.m[i]:end.ind.m[i]
    for (j in 1:q){
      ind2=start.ind.n[j]:end.ind.n[j]
      G00.mat[[i,j]]=G00[ind1,ind2]*sigma.mat[i,j]
      R00.mat[[i,j]]=R00[ind1,ind2]*sigma.mat[i,j]
      C00.mat[[i,j]]=C00[ind1,ind2]*sigma.mat[i,j]
      I00.mat[[i,j]]=I00[ind1,ind2]*sigma.mat[i,j]
      S00.mat[[i,j]]=G00.mat[[i,j]]+R00.mat[[i,j]]+C00.mat[[i,j]]+I00.mat[[i,j]]
    }
  }
  
  return(list(X=data, S=S00.mat,
              G=G00.mat, R=R00.mat, C=C00.mat, I=I00.mat,
              sigma.mat=sigma.mat, n.vec=nvec,m.vec=mvec))
}

impute.BIDIFAC=function(data, 
                        rmt=T, sigma=NULL,
                        pbar=TRUE,
                        start=NULL, max.iter.impute=20, 
                        eps.impute=1e-3,...){
  dim.data=dim(data)
  p=dim.data[1]; q=dim.data[2]
  
  dim.list=do.call(cbind,lapply(data, dim))
  mvec=apply(matrix(dim.list[1,],p),1,unique)
  nvec=apply(matrix(dim.list[2,],p),2,unique)
  if (class(mvec)=="list" ) stop("the number of rows do not match")
  if (class(nvec)=="list" ) stop("the number of columns do not match")
  
  impute.index=matrix(list(), nrow = p, ncol=q)
  if (is.null(sigma)) sigma=matrix(1,p,q)

  for (i in 1:p){
    for (j in 1:q){
      fillmat=fill.matrix(data[[i,j]])
      impute.index[[i,j]]=fillmat$na.ind
      if (rmt) sigma[i,j]=sigma.rmt(fillmat$X.fill)
      data[[i,j]]=fillmat$X.fill/sigma[i,j]
    }
  }
  
  bool2=TRUE
  count2=1; impute.vec=0
  if (pbar) pb = txtProgressBar(min = 0, max=max.iter.impute, initial=0, char="-", style = 3)
  while (bool2){
    if (pbar){  setTxtProgressBar(pb, count2)  }
    impute.vec.old=impute.vec
    fit=BIDIFAC(data, rmt=F, sigma=matrix(1,p,q),start=start, pbar = F, ...)
    
    start=list(
      data.rearrange(fit$G)$out, data.rearrange(fit$R)$out,
      data.rearrange(fit$C)$out, data.rearrange(fit$I)$out)
    
    impute.vec=NULL
    for (i in 1:p){
      for (j in 1:q){
        imp=fit$S[[i,j]][impute.index[[i,j]]]
        data[[i,j]][impute.index[[i,j]]]=imp
        impute.vec=c(impute.vec,imp)
      }
    }
    
    crit2=frob(impute.vec.old-impute.vec)/length(impute.vec)
    if (crit2<eps.impute){ bool2=FALSE }
    else if (count2==max.iter.impute){ 
      bool2=FALSE 
      cat("\n The algorithm did not converge. \n")
      cat("Try to increase max.iter.impute or relax eps/eps.impute.")
      }
    else{ count2=count2+1 }
  }
  
  for (i in 1:p){
    for (j in 1:q){
      fit$X[[i,j]]=fit$X[[i,j]]*sigma[i,j]
      fit$S[[i,j]]=fit$S[[i,j]]*sigma[i,j]
      fit$G[[i,j]]=fit$G[[i,j]]*sigma[i,j]
      fit$R[[i,j]]=fit$R[[i,j]]*sigma[i,j]
      fit$C[[i,j]]=fit$C[[i,j]]*sigma[i,j]
      fit$I[[i,j]]=fit$I[[i,j]]*sigma[i,j]
    }
  }
  fit$sigma.mat=sigma
  
  return(fit)
}

summary.BIDIFAC=function(fit){
  X00=fit$X
  dimm=dim(X00)
  p=dimm[1]; q=dimm[2]
  out=NULL
  rnames=NULL
  for (i in 1:p){
    for (j in 1:q){
      vec=rep(NA,8)
      rnames=c(rnames, paste0("X",i,",",j))      
      frob.Xij=frob(fit$X[[i,j]])
      vec[1]=1-frob(fit$G[[i,j]]-fit$X[[i,j]])/frob.Xij
      vec[2]=1-frob(fit$R[[i,j]]-fit$X[[i,j]])/frob.Xij
      vec[3]=1-frob(fit$C[[i,j]]-fit$X[[i,j]])/frob.Xij
      vec[4]=1-frob(fit$I[[i,j]]-fit$X[[i,j]])/frob.Xij
      vec[5]=1-frob(fit$G[[i,j]]+fit$R[[i,j]]-fit$X[[i,j]])/frob.Xij
      vec[6]=1-frob(fit$G[[i,j]]+fit$C[[i,j]]-fit$X[[i,j]])/frob.Xij
      vec[7]=1-frob(fit$G[[i,j]]+fit$R[[i,j]]+fit$C[[i,j]]-fit$X[[i,j]])/frob.Xij
      vec[8]=1-frob(fit$S[[i,j]]-fit$X[[i,j]])/frob.Xij
      out=rbind(out, vec)
    }
  }
  
  rownames(out)=  rnames
  colnames(out)= c("G","R","C","I","GR","GC", "GRC", "GRCI")
  return(out)
}

show.image = function(Image,sd=NULL,ylab=''){
  if (is.null(sd)){s=sd(Image)}
  else{s=sd}
  
  lower = mean(Image)-3*s
  upper = mean(Image)+3*s
  Image[Image<lower] = lower
  Image[Image>upper] = upper
  image(x=1:dim(Image)[2], y=1:dim(Image)[1], z=t(Image), zlim = c(lower,upper),axes=FALSE,col=bluered(100),xlab="",ylab=ylab)
} 

