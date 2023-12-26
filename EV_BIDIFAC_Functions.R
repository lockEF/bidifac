
##EB.glob.o finds EVB solution for a single matrix
eb.glob.o<-function(X,sig=NULL,kap=NULL){
  #thm 2 in https://proceedings.neurips.cc/paper_files/paper/2012/file/26337353b7962f533d78c762373b3318-Paper.pdf
  D=dim(X)[1]
  N=dim(X)[2]
  if(is.null(sig)) sig=get_sigma(X)
  if(is.null(kap)) kap=get_kappa(X)
  gambar=sqrt(sig)*sqrt(D+N+sqrt(D*N)*(kap+1/kap))
  svd.x=svd(X)
  ds=svd.x$d
  ind1=ds>gambar
  dvec1=rep(0,length(ds))
  dvec1[ind1]=(ds[ind1]/2)*(1-(D+N)*sig/ds[ind1]^2+sqrt((1-(D+N)*sig/ds[ind1]^2)-(4*D*N*sig^2)/ds[ind1]^4))
  A= svd.x$u%*%diag(dvec1)%*%t(svd.x$v)
  #####unpack solution
  return(A)
}

##eb.glob.imp finds EVB solution for a single matrix, with missing data imputation
eb.glob.imp<-function(X,kap=NULL,max.iter=100,conv.thresh=0.1){
  if(is.null(kap)) kap=get_kappa(X)
  isna=is.na(X)
  X[isna]=0
  sig.est=mean(X^2)
  A.prev=0
  A.est=nn.app.imp(X)
  for(i in 1:max.iter){
    A=sqrt(sig.est)*nn.app(X/sqrt(sig.est))
    sig.est=mean((X-A)^2) 
    X[isna]=A[isna]
    if(sum((A.prev-A)^2)<conv.thresh) break
    print(sum((A.prev-A)^2))
    A.prev=A
  }
  for(i in 1:max.iter){
    A=eb.glob.o2(X,sig=sig.est)
    X[isna]=A[isna]
    sig.est=get_sigma(X)*length(X)/(length(X)-sum(isna))
    if(sum((A.prev-A)^2)<conv.thresh) break
    print(sum((A.prev-A)^2))
    A.prev=A
  }
  return(A)
}

#glob.bidi.o identifies decomposition for bidimensionally linked matrices, without missingness
glob.bidi.o <- function(X,p.ind,n.ind,p.ind.list,n.ind.list,kappa, conv.thresh=0.1,max.iter=1000){ 
  I=length(p.ind)
  J=length(n.ind)
  K = length(p.ind.list) 
  S=list()
  for(k in 1:K) S[[k]]=array(rep(0,prod(dim(X))),dim=dim(X))
  sigma2mat=matrix(rep(1,I*J),nrow=I,ncol=J)
  Xs=X
  for(i in 1:I){ for(j in 1:J){
    sigma2mat[i,j]=get_sigma(X[p.ind[[i]],n.ind[[j]]])
    Xs[p.ind[[i]],n.ind[[j]]]=X[p.ind[[i]],n.ind[[j]]]/sqrt(sigma2mat[i,j])
  }}
  for(iter in 1:max.iter){
    S.prev=S
    for(k in 1:K){
      Resid=Xs-Reduce('+',S[-k])
      S[[k]][p.ind.list[[k]],n.ind.list[[k]]]=nn.app(Resid[p.ind.list[[k]],n.ind.list[[k]]])
    }
    if(sum((unlist(S)-unlist(S.prev))^2)<conv.thresh){
      break}
    print(iter)
    print(sum((unlist(S)-unlist(S.prev))^2))
  }
  for(iter in 1:max.iter){
    S.prev=S
    for(k in 1:K){
      Resid=Xs-Reduce('+',S[-k])
      S[[k]][p.ind.list[[k]],n.ind.list[[k]]]=eb.glob.o(Resid[p.ind.list[[k]],n.ind.list[[k]]],sig=1)
    }
    if(sum((unlist(S)-unlist(S.prev))^2)<conv.thresh){
      break}
    print(iter)
    print(sum((unlist(S)-unlist(S.prev))^2))
  }
  for(i in 1:I){ for(j in 1:J){
    for(k in 1:K) S[[k]][p.ind[[i]],n.ind[[j]]]=sqrt(sigma2mat[i,j])*S[[k]][p.ind[[i]],n.ind[[j]]]
  }}
  return(list(S=S,sigma2mat=sigma2mat))
}

#glob.bidi.imp identifies decomposition for bidimensionally linked matrices, with missing data imputation
glob.bidi.imp <- function(X,p.ind,n.ind,p.ind.list,n.ind.list, conv.thresh=0.1,max.iter=1000){
  I=length(p.ind)
  J=length(n.ind)
  K = length(p.ind.list) 
  isna=is.na(X)
  S=list()
  for(k in 1:K) S[[k]]=array(rep(0,prod(dim(X))),dim=dim(X))
  sigma2mat=matrix(rep(1,I*J),nrow=I,ncol=J)
  S.scale=S
  X[isna]=0
  Xs=X
  for(i in 1:I){ for(j in 1:J){
    sigma2mat[i,j]=mean(X[p.ind[[i]],n.ind[[j]]]^2)
    Xs[p.ind[[i]],n.ind[[j]]]=X[p.ind[[i]],n.ind[[j]]]/sqrt(sigma2mat[i,j])
  }}
  for(iter in 1:max.iter){
    S.prev=S
    for(k in 1:K){
      Resid=Xs-Reduce('+',S[-k])
      S[[k]][p.ind.list[[k]],n.ind.list[[k]]]=nn.app(Resid[p.ind.list[[k]],n.ind.list[[k]]])
      #    S[[k]][p.ind.list[[k]],n.ind.list[[k]]][,colSums(isna[p.ind.list[[k]],n.ind.list[[k]]])==length(p.ind.list[[k]])]=0
      #    S[[k]][p.ind.list[[k]],n.ind.list[[k]]][rowSums(isna[p.ind.list[[k]],n.ind.list[[k]]])==length(n.ind.list[[k]]),]=0
    }
    sig=Reduce('+',S)
    Xs[isna]=sig[isna]
    for(i in 1:I){ for(j in 1:J){
      sprev=sigma2mat[i,j]
      resid=X[p.ind[[i]],n.ind[[j]]]-sqrt(sigma2mat[i,j])*sig[p.ind[[i]],n.ind[[j]]]
      sigma2mat[i,j]=mean(resid^2)#mean(resid[!isna[p.ind[[i]],n.ind[[j]]]]^2)#mean(resid^2) #update sigma2
      X[p.ind[[i]],n.ind[[j]]][isna[p.ind[[i]],n.ind[[j]]]]=sqrt(sigma2mat[i,j])*sig[p.ind[[i]],n.ind[[j]]][isna[p.ind[[i]],n.ind[[j]]]] ##need to do this after calculating variance, to account for variability in different solutions 
      Xs[p.ind[[i]],n.ind[[j]]]=X[p.ind[[i]],n.ind[[j]]]/sqrt(sigma2mat[i,j])
      for(k in 1:K){
        S[[k]][p.ind[[i]],n.ind[[j]]]=sqrt(sprev/sigma2mat[i,j])*S[[k]][p.ind[[i]],n.ind[[j]]]
        S.scale[[k]][p.ind[[i]],n.ind[[j]]]=sqrt(sigma2mat[i,j])*S[[k]][p.ind[[i]],n.ind[[j]]]
      }}}
    if(sum((unlist(S)-unlist(S.prev))^2)<conv.thresh){
      break}
    print(iter)
    print(sum((unlist(S)-unlist(S.prev))^2))
  }
  for(iter in 1:max.iter){
    S.prev=S
    S.scale.prev=S.scale
    for(k in 1:K){
      Resid=Xs-Reduce('+',S[-k])
      S[[k]][p.ind.list[[k]],n.ind.list[[k]]]=eb.glob.o2(Resid[p.ind.list[[k]],n.ind.list[[k]]],sig=1)
      #set to 0 anything that is missing entire row or column
      S[[k]][p.ind.list[[k]],n.ind.list[[k]]][,colSums(isna[p.ind.list[[k]],n.ind.list[[k]]])==length(p.ind.list[[k]])]=0
      S[[k]][p.ind.list[[k]],n.ind.list[[k]]][rowSums(isna[p.ind.list[[k]],n.ind.list[[k]]])==length(n.ind.list[[k]]),]=0
    }
    sig=Reduce('+',S)
    Xs[isna]=sig[isna]
    sig_unnorm=sig
    for(i in 1:I){ for(j in 1:J){
      sig_unnorm[p.ind[[i]],n.ind[[j]]]=sqrt(sigma2mat[i,j])*sig[p.ind[[i]],n.ind[[j]]] }}
    X[isna]=sig_unnorm[isna]
    for(i in 1:I){ for(j in 1:J){
      sprev=sigma2mat[i,j]
      cor_fac=length(X[p.ind[[i]],n.ind[[j]]])/(length(X[p.ind[[i]],n.ind[[j]]])-sum(isna[p.ind[[i]],n.ind[[j]]]))
      sigma2mat[i,j]=get_sigma(X[p.ind[[i]],n.ind[[j]]])*cor_fac #update sigma2
      X[p.ind[[i]],n.ind[[j]]][isna[p.ind[[i]],n.ind[[j]]]]=sqrt(sigma2mat[i,j])*sig[p.ind[[i]],n.ind[[j]]][isna[p.ind[[i]],n.ind[[j]]]] ##need to do this after calculating variance, to account for variability in different solutions 
      Xs[p.ind[[i]],n.ind[[j]]]=X[p.ind[[i]],n.ind[[j]]]/sqrt(sigma2mat[i,j])
      for(k in 1:K){
        S[[k]][p.ind[[i]],n.ind[[j]]]=sqrt(sprev/sigma2mat[i,j])*S[[k]][p.ind[[i]],n.ind[[j]]]
        S.scale[[k]][p.ind[[i]],n.ind[[j]]]=sqrt(sigma2mat[i,j])*S[[k]][p.ind[[i]],n.ind[[j]]]
      }}}
    if(sum((unlist(S.scale)-unlist(S.scale.prev))^2)<conv.thresh){
      break}
    print(iter)
    print(sum((unlist(S.scale)-unlist(S.scale.prev))^2))
    print(sum((Reduce('+',S)-Reduce('+',S.prev))^2))
    sig.prev=sig
    sigma2mat.prev=sigma2mat
  }
  for(i in 1:I){ for(j in 1:J){
    for(k in 1:K) S[[k]][p.ind[[i]],n.ind[[j]]]=sqrt(sigma2mat[i,j])*S[[k]][p.ind[[i]],n.ind[[j]]]
  }}
  return(list(S=S,sigma2mat=sigma2mat))
}

#helper functions!
get_kappa <- function(X,lb=0.5,ub=4)
{
  N=ncol(X)
  D=nrow(X)
  alpha=D/N 
  k.vec=seq(lb,ub,length.out=50000)
  EE=log(sqrt(alpha)*k.vec+1)/(sqrt(alpha)*k.vec)+log(k.vec/sqrt(alpha)+1)/(k.vec/sqrt(alpha))-1
  kappa=k.vec[which.min(EE^2)]
  return(kappa)
}

psi<-function(x,alpha,kappa,H){
  xbar=1+alpha+sqrt(alpha)*(kappa+1/kappa)
  res=x-log(x)
  ind=c(1:H)[x[1:H]>=xbar]
  k.func=(1/(2*sqrt(alpha)))*((x[ind]-(1+alpha))+sqrt((x[ind]-(1+alpha))^2-4*alpha))
  res[ind]=x[ind]-log(x[ind])+log(sqrt(alpha)*k.func+1)+alpha*log(k.func/sqrt(alpha)+1)-sqrt(alpha)*k.func
  return(res)
}

get_sigma <- function(X,lb=0.1,ub=NULL,kappa=NULL){
  if(is.null(ub)) ub = mean(X^2)
  if(is.null(kappa)) kappa=get_kappa(X)
  D=max(dim(X))
  N=min(dim(X))
  alpha=N/D
  svd.x=svd(X)
  sig.vec=seq(lb,ub,length.out=5000)
  om=rep(0,length(sig.vec))
  for(i in 1:length(sig.vec)){
 #   om[i]=sum(psi(svd.x$d^2/(D*sig.vec[i]),alpha,kappa,H=min(D,N)))##should be D (bigger), not N?
    om[i]=sum(psi(svd.x$d^2/(D*sig.vec[i]),alpha,kappa,H=sum(svd.x$d>0.001)))##should be D (bigger), not N?
  }
  sig=sig.vec[which.min(om)]
  return(sig)
}
 
 
