library(gplots)
library(MASS)
library(Matrix)
#library(gpuR)
library(denoiseR)

frob<-function(X){
  #Also works for vclMatrix
  sum(X^2,na.rm=T)
}


findNoiseSD.RMT<-function(X){
  estim_sigma(X,method="MAD")
}


generateRandomMatrix<-function(nrows, ncols, gpu, svd=F){
  out=matrix(rnorm(nrows*ncols),nrows)
  if (svd){out=svd(out)$u}
  if (gpu){out=vclMatrix(out) }
  return(out)
}

myridge<-function(U,Y,lambda, transpose=FALSE){
  X=crossprod(U)+lambda*diag(ncol(U))
  X.solve=solve(X)
  if (transpose){
    out= crossprod(Y,U)%*%X.solve
  }
  else{
    out=X.solve%*% crossprod(U,Y)
  }
  return(out)
}

BIDIFAC<-function(X1,X2, noise.est="RMT", sigma1=NULL, sigma2=NULL,
                   U.start=NULL,V.start=NULL, out=FALSE,
                   eps=1e-5, max.iter=1000, gpu=FALSE, pbar=TRUE, seed=NULL, ...){
  #It fits BID.jive to the bidimensionally linked matrices X1=(X11,X12), X2=(X21,X22)
  #noise.RMT: if TRUE, BID.jive() detects noise variance by random matrix theory
  #sigma1,sigma2: a vector, where each cell is the standard deviation of noise matrices X1 and X2
  #       only works when noise.RMT=FALSE
  #       If noise.RMT=FALSE and sigma=NULL, BID.jive uses sample standard deviation.
  #eps: threshold to stop the iteration
  #max.iter: total number of iterations
  #pbar: progress bar
  #gpu: wheter or not gpu is used to handle big matries. The computer or system must have GPU installed.
  
  if (!is.null(seed)){set.seed(seed)}
  
  X11.original = X1[[1]]; X12.original = X1[[2]]; 
  X21.original = X2[[1]]; X22.original = X2[[2]];
  
  n1 = nrow(X11.original); n2 = nrow(X21.original)
  m1 = ncol(X11.original); m2 = ncol(X12.original)
  rkG<-min(n1,n2,m1,m2)
  rkC1<-min(n1,n2,m1)
  rkC2<-min(n1,n2,m2)
  rkR1<-min(n1,m1,m2)
  rkR2<-min(n2,m1,m2)
  rkI11<-min(n1,m1)
  rkI12<-min(n1,m2)
  rkI21<-min(n2,m1)
  rkI22<-min(n2,m2)
  
  
  if (is.null(noise.est)){
    if (is.null(sigma1) || is.null(sigma2)){
      sigma11=sd(X11.original); sigma12=sd(X12.original); 
      sigma21=sd(X21.original); sigma22=sd(X22.original)
    }
    else{ sigma11=sigma1[1]; sigma12=sigma1[2]; sigma21=sigma2[1]; sigma22=sigma2[2] }
  }
  else if (noise.est=="RMT"){
    sigma11=findNoiseSD.RMT(X11.original); sigma12=findNoiseSD.RMT(X12.original); 
    sigma21=findNoiseSD.RMT(X21.original); sigma22=findNoiseSD.RMT(X22.original);
  } 
  
  
  
  X11 = X11.original/sigma11; X11.t = t(X11); 
  X12 = X12.original/sigma12; X12.t = t(X12);
  X21 = X21.original/sigma21; X21.t = t(X21); 
  X22 = X22.original/sigma22; X22.t = t(X22);
  
  lambda1.R = (sqrt(n1)+sqrt(m1+m2))
  lambda2.R = (sqrt(n2)+sqrt(m1+m2))
  lambda1.C = (sqrt(n1+n2)+sqrt(m1))
  lambda2.C = (sqrt(n1+n2)+sqrt(m2))
  lambda11.I = (sqrt(n1)+sqrt(m1))
  lambda12.I = (sqrt(n1)+sqrt(m2))
  lambda21.I = (sqrt(n2)+sqrt(m1))
  lambda22.I = (sqrt(n2)+sqrt(m2))
  
  bool<-TRUE
  crit<-NULL
  count<-1
  
  if (is.null(U.start) || is.null(V.start)){
    U1.G<-generateRandomMatrix(n1,rkG,gpu=gpu, ...) ;
    U1.R<-generateRandomMatrix(n1,rkR1,gpu=gpu, ...);
    U11.C<-generateRandomMatrix(n1,rkC1,gpu=gpu, ...);
    U12.C<-generateRandomMatrix(n1,rkC2,gpu=gpu, ...);
    U11.I<-generateRandomMatrix(n1,rkI11,gpu=gpu, ...);
    U12.I<-generateRandomMatrix(n1,rkI12,gpu=gpu, ...);
    
    U2.G<-generateRandomMatrix(n2,rkG,gpu=gpu, ...);
    U2.R<-generateRandomMatrix(n2,rkR2,gpu=gpu, ...);
    U21.C<-generateRandomMatrix(n2,rkC1,gpu=gpu, ...);
    U22.C<-generateRandomMatrix(n2,rkC2,gpu=gpu, ...);
    U21.I<-generateRandomMatrix(n2,rkI21,gpu=gpu, ...);
    U22.I<- generateRandomMatrix(n2,rkI22,gpu=gpu, ...);
    
    V1.G<-generateRandomMatrix(m1,rkG,gpu=gpu, ...);
    V1.C<-generateRandomMatrix(m1,rkC1,gpu=gpu, ...);
    V11.R<-generateRandomMatrix(m1,rkR1,gpu=gpu, ...);
    V21.R<-generateRandomMatrix(m1,rkR2,gpu=gpu, ...);
    V11.I<-generateRandomMatrix(m1,rkI11,gpu=gpu, ...);
    V21.I<-generateRandomMatrix(m1,rkI21,gpu=gpu, ...);
    
    V2.G<-generateRandomMatrix(m2,rkG,gpu=gpu, ...);
    V2.C<-generateRandomMatrix(m2,rkC2,gpu=gpu, ...);
    V12.R<-generateRandomMatrix(m2,rkR1,gpu=gpu, ...);
    V22.R<-generateRandomMatrix(m2,rkR2,gpu=gpu, ...);
    V12.I<-generateRandomMatrix(m2,rkI12,gpu=gpu, ...);
    V22.I<-generateRandomMatrix(m2,rkI22,gpu=gpu, ...);    
  }
  else{
    U1.G<-U.start[[1]];
    U1.R<-U.start[[2]];
    U11.C<-U.start[[3]];
    U12.C<-U.start[[4]];
    U11.I<-U.start[[5]];
    U12.I<-U.start[[6]];
    
    U2.G<-U.start[[7]];
    U2.R<-U.start[[8]];
    U21.C<-U.start[[9]];
    U22.C<-U.start[[10]];
    U21.I<-U.start[[11]];
    U22.I<-U.start[[12]];
    
    V1.G<-V.start[[1]];
    V11.R<-V.start[[2]];
    V21.R<-V.start[[3]];
    V1.C<-V.start[[4]];
    V11.I<-V.start[[5]];
    V21.I<-V.start[[6]];
    
    V2.G<-V.start[[7]];
    V12.R<-V.start[[8]];
    V22.R<-V.start[[9]];
    V2.C<-V.start[[10]];
    V12.I<-V.start[[11]];
    V22.I<-V.start[[12]];    
  }
  
  
  
  pred11.G=pred11.C=pred11.R=pred11.I=
    pred12.G=pred12.C=pred12.R=pred12.I=
    pred21.G=pred21.C=pred21.R=pred21.I=
    pred22.G=pred22.C=pred22.R=pred22.I=0
  
  crit0<-0
  if (pbar==TRUE) pb <- txtProgressBar(min = 0, max=max.iter, initial=0, char="-", style = 3)
  while (bool){
    if (pbar==TRUE){  setTxtProgressBar(pb, count)  }
    crit0.old = crit0
    U1.G.old = U1.G
    U2.G.old = U2.G
    V1.G.old = V1.G
    V2.G.old = V2.G
    
    V1.C.old = V1.C
    V2.C.old = V2.C
    U11.C.old = U11.C
    U21.C.old = U21.C
    U12.C.old = U12.C
    U22.C.old = U22.C
    
    U1.R.old = U1.R
    U2.R.old = U2.R
    V11.R.old = V11.R
    V12.R.old = V12.R
    V21.R.old = V21.R
    V22.R.old = V22.R
    
    U11.I.old = U11.I
    U12.I.old = U12.I
    U21.I.old = U21.I
    U22.I.old = U22.I
    
    V11.I.old = V11.I
    V12.I.old = V12.I
    V21.I.old = V21.I
    V22.I.old = V22.I
    
    pred11.G.old = pred11.G; pred11.R.old = pred11.R; pred11.C.old = pred11.C; pred11.I.old = pred11.I
    pred12.G.old = pred12.G; pred12.R.old = pred12.R; pred12.C.old = pred12.C; pred12.I.old = pred12.I
    pred21.G.old = pred21.G; pred21.R.old = pred21.R; pred21.C.old = pred21.C; pred21.I.old = pred21.I
    pred22.G.old = pred22.G; pred22.R.old = pred22.R; pred22.C.old = pred22.C; pred22.I.old = pred22.I
    
    
    #Update V1.G and V2.G
    Y = rbind(X11- tcrossprod(U11.C,V1.C) - tcrossprod(U1.R, V11.R)- tcrossprod(U11.I,V11.I),
              X21-tcrossprod(U21.C,V1.C) - tcrossprod(U2.R, V21.R)- tcrossprod(U21.I,V21.I))
    V1.G = myridge(rbind(U1.G, U2.G), Y, lambda1.C, transpose = T)
    
    Y = rbind(X12- tcrossprod(U12.C, V2.C) - tcrossprod(U1.R, V12.R) - tcrossprod(U12.I,V12.I),
              X22-tcrossprod(U22.C,V2.C) - tcrossprod(U2.R, V22.R) - tcrossprod(U22.I,V22.I))
    V2.G = myridge(rbind(U1.G, U2.G),Y, lambda2.C, transpose = T)
    
    #Update U1.J and U2.J
    Y = rbind(X11.t-tcrossprod(V1.C,U11.C)- tcrossprod(V11.R, U1.R) - tcrossprod(V11.I,U11.I),
              X12.t-tcrossprod(V2.C,U12.C)- tcrossprod(V12.R, U1.R) - tcrossprod(V12.I,U12.I))
    U1.G = myridge(rbind(V1.G, V2.G), Y, lambda1.R, transpose = T)
    
    Y = rbind(X21.t-tcrossprod(V1.C,U21.C)- tcrossprod(V21.R, U2.R) - tcrossprod(V21.I,U21.I),
              X22.t-tcrossprod(V2.C,U22.C)-tcrossprod(V22.R, U2.R) -tcrossprod(V22.I,U22.I))
    U2.G = myridge(rbind(V1.G, V2.G),Y, lambda2.R, transpose = T)
    
    #Update V1.C and V2.C
    Y = rbind(X11-tcrossprod(U1.G,V1.G)-tcrossprod(U1.R, V11.R) - tcrossprod(U11.I,V11.I),
              X21-tcrossprod(U2.G,V1.G)-tcrossprod(U2.R, V21.R)-tcrossprod(U21.I,V21.I))
    V1.C = myridge(rbind(U11.C, U21.C), Y, lambda1.C, transpose = T)
    
    Y = rbind(X12-tcrossprod(U1.G,V2.G)-tcrossprod(U1.R, V12.R)-tcrossprod(U12.I,V12.I),
              X22-tcrossprod(U2.G,V2.G)-tcrossprod(U2.R, V22.R)-tcrossprod(U22.I,V22.I))
    V2.C = myridge(rbind(U12.C, U22.C), Y, lambda2.C, transpose = T)
    
    #Update U11.C, U12.C, U21.C, U22.C
    Y = X11.t-tcrossprod(V1.G,U1.G)-tcrossprod(V11.R, U1.R)-tcrossprod(V11.I,U11.I)
    U11.C = myridge(V1.C, Y, lambda11.I, transpose = T)
    
    Y = X12.t-tcrossprod(V2.G,U1.G)-tcrossprod(V12.R, U1.R)-tcrossprod(V12.I,U12.I)
    U12.C = myridge(V2.C, Y, lambda12.I, transpose = T)
    
    Y = X21.t-tcrossprod(V1.G,U2.G)-tcrossprod(V21.R, U2.R)-tcrossprod(V21.I,U21.I)
    U21.C = myridge(V1.C,Y, lambda21.I, transpose = T)
    
    Y = X22.t-tcrossprod(V2.G,U2.G)-tcrossprod(V22.R, U2.R)-tcrossprod(V22.I,U22.I)
    U22.C = myridge(V2.C, Y, lambda22.I, transpose = T)
    
    #Update V11.R, V12.R, V21.R, V22.R
    Y = X11-tcrossprod(U1.G,V1.G)-tcrossprod(U11.C,V1.C)-tcrossprod(U11.I,V11.I)
    V11.R = myridge(U1.R, Y, lambda11.I, transpose = T)
    
    Y = X12-tcrossprod(U1.G,V2.G)-tcrossprod(U12.C,V2.C)-tcrossprod(U12.I,V12.I)
    V12.R = myridge(U1.R, Y, lambda12.I, transpose = T)
    
    Y = X21-tcrossprod(U2.G,V1.G)-tcrossprod(U21.C,V1.C)-tcrossprod(U21.I,V21.I)
    V21.R = myridge(U2.R, Y, lambda21.I, transpose = T)
    
    Y = X22-tcrossprod(U2.G,V2.G)-tcrossprod(U22.C,V2.C)-tcrossprod(U22.I,V22.I)
    V22.R = myridge(U2.R, Y, lambda22.I, transpose = T)
    
    #Update U1.R and U2.R
    Y = rbind(X11.t-tcrossprod(V1.G,U1.G)-tcrossprod(V1.C,U11.C)-tcrossprod(V11.I,U11.I),
              X12.t-tcrossprod(V2.G,U1.G)-tcrossprod(V2.C,U12.C)-tcrossprod(V12.I,U12.I))
    U1.R = myridge(rbind(V11.R, V12.R), Y, lambda1.R, transpose = T)
    
    Y = rbind(X21.t-tcrossprod(V1.G,U2.G)-tcrossprod(V1.C,U21.C)-tcrossprod(V21.I,U21.I),
              X22.t-tcrossprod(V2.G,U2.G)-tcrossprod(V2.C,U22.C)-tcrossprod(V22.I,U22.I))
    U2.R = myridge(rbind(V21.R, V22.R),Y, lambda2.R, transpose = T)
    
    #Update V11.I, V12.I, V21.I, and V22.I
    Y = X11-tcrossprod(U1.G,V1.G)-tcrossprod(U1.R, V11.R)-tcrossprod(U11.C,V1.C)
    V11.I = myridge(U11.I,Y, lambda11.I, transpose = T)
    
    Y = X12-tcrossprod(U1.G,V2.G)-tcrossprod(U1.R, V12.R)-tcrossprod(U12.C,V2.C)
    V12.I = myridge(U12.I, Y, lambda12.I, transpose = T)
    
    Y = X21-tcrossprod(U2.G,V1.G)-tcrossprod(U2.R, V21.R)-tcrossprod(U21.C,V1.C)
    V21.I = myridge(U21.I, Y, lambda21.I, transpose = T)
    
    Y = X22-tcrossprod(U2.G,V2.G)-tcrossprod(U2.R, V22.R)-tcrossprod(U22.C,V2.C)
    V22.I = myridge(U22.I, Y, lambda22.I, transpose = T)
    
    #Update U11.I, U12.I, U21.I, and U22.I
    Y = X11.t-tcrossprod(V1.G,U1.G)-tcrossprod(V11.R, U1.R)-tcrossprod(V1.C,U11.C)
    U11.I = myridge(V11.I, Y, lambda11.I, transpose = T)
    
    Y = X12.t-tcrossprod(V2.G,U1.G)-tcrossprod(V12.R, U1.R)-tcrossprod(V2.C,U12.C)
    U12.I = myridge(V12.I, Y, lambda12.I, transpose = T)
    
    Y = X21.t-tcrossprod(V1.G,U2.G)-tcrossprod(V21.R, U2.R)-tcrossprod(V1.C,U21.C)
    U21.I = myridge(V21.I,Y, lambda21.I, transpose = T)
    
    Y = X22.t-tcrossprod(V2.G,U2.G)-tcrossprod(V22.R, U2.R)-tcrossprod(V2.C,U22.C)
    U22.I = myridge(V22.I,Y, lambda22.I, transpose = T)
    
    pred11.G = tcrossprod(U1.G,V1.G); pred11.I = tcrossprod(U11.I,V11.I); pred11.R = tcrossprod(U1.R,V11.R); pred11.C = tcrossprod(U11.C,V1.C); pred11 = pred11.G+pred11.C+pred11.R+pred11.I
    pred12.G = tcrossprod(U1.G,V2.G); pred12.I = tcrossprod(U12.I,V12.I); pred12.R = tcrossprod(U1.R,V12.R); pred12.C = tcrossprod(U12.C,V2.C); pred12 = pred12.G+pred12.C+pred12.R+pred12.I
    pred21.G = tcrossprod(U2.G,V1.G); pred21.I = tcrossprod(U21.I,V21.I); pred21.R = tcrossprod(U2.R,V21.R); pred21.C = tcrossprod(U21.C,V1.C); pred21 = pred21.G+pred21.C+pred21.R+pred21.I
    pred22.G = tcrossprod(U2.G,V2.G); pred22.I = tcrossprod(U22.I,V22.I); pred22.R = tcrossprod(U2.R,V22.R); pred22.C = tcrossprod(U22.C,V2.C); pred22 = pred22.G+pred22.C+pred22.R+pred22.I
    
    crit0 = frob(X11-pred11)+frob(X12-pred12)+frob(X21-pred21)+frob(X22-pred22)+
      lambda11.I*(frob(U11.I)+frob(V11.I)+frob(U11.C)+frob(V11.R) )+
      lambda12.I*(frob(U12.I)+frob(V12.I)+frob(U12.C)+frob(V12.R) )+
      lambda21.I*(frob(U21.I)+frob(V21.I)+frob(U21.C)+frob(V21.R) )+
      lambda22.I*(frob(U22.I)+frob(V22.I)+frob(U22.C)+frob(V22.R) )+
      lambda1.C*(frob(V1.C)+frob(V1.G))+lambda2.C*(frob(V2.C)+frob(V2.G))+
      lambda1.R*(frob(U1.R)+frob(U1.G))+lambda2.R*(frob(U2.R)+frob(U2.G))
    
    crit=c(crit,crit0)
    
    if (abs(crit0.old-crit0)<eps){bool=FALSE}
    else if (count==max.iter){bool=FALSE}
    else{ count = count+1 }
  }
  
  pred11 = pred11*sigma11
  pred11.G = pred11.G*sigma11
  pred11.C = pred11.C*sigma11
  pred11.R = pred11.R*sigma11
  pred11.I = pred11.I*sigma11
  
  pred12 = pred12*sigma12
  pred12.G = pred12.G*sigma12
  pred12.C = pred12.C*sigma12
  pred12.R = pred12.R*sigma12
  pred12.I = pred12.I*sigma12
  
  pred21 = pred21*sigma21
  pred21.G = pred21.G*sigma21
  pred21.C = pred21.C*sigma21
  pred21.R = pred21.R*sigma21
  pred21.I = pred21.I*sigma21
  
  pred22 = pred22*sigma22
  pred22.G = pred22.G*sigma22
  pred22.C = pred22.C*sigma22
  pred22.R = pred22.R*sigma22
  pred22.I = pred22.I*sigma22
  
  if (out){
    U.out=list(U1.G, U1.R, U11.C, U12.C, U11.I, U12.I, U2.G, U2.R, U21.C, U22.C, U21.I, U22.I)
    V.out=list(V1.G, V11.R, V21.R, V1.C, V11.I, V21.I, V2.G, V12.R, V22.R, V2.C, V12.I, V22.I)
  }
  else{
    U.out=V.out=NULL    
  }
  
  return(list(X1.original=list(X11.original, X12.original),
              X2.original=list(X21.original, X22.original),
              X1.pred=list(pred11, pred12),
              X2.pred=list(pred21, pred22),
              X1.global=list(pred11.G, pred12.G),
              X2.global=list(pred21.G, pred22.G),
              X1.row=list(pred11.R, pred12.R),
              X2.row=list(pred21.R, pred22.R),
              X1.col=list(pred11.C, pred12.C),
              X2.col=list(pred21.C, pred22.C),
              X1.individual=list(pred11.I, pred12.I),
              X2.individual=list(pred21.I, pred22.I),
              sigma1=c(sigma11,sigma12),
              sigma2=c(sigma21,sigma22),
              U.out=U.out,
              V.out=V.out,
              crit=crit
  ))
}



impute.BIDIFAC<-function(X1,X2, noise.est="RMT", sigma1=NULL,sigma2=NULL,
                          U.start=NULL, V.start=NULL,
                          max.iter.impute=20, eps.impute=1e-3, pbar=TRUE, ...){
  # imputes missing value using jive2() iteratively
  # produces imputed matrices as well as jive decompositions
  # HT : whether or not hard-thresholding results are reported
  # m.iter : number of iterations for imputing misisng values
  X1.original=X1
  X2.original=X2
  X11<-X1[[1]]; X12<-X1[[2]]; X21<-X2[[1]]; X22=X2[[2]]
  n1<-nrow(X11); n2<-nrow(X21)
  m1<-ncol(X11); m2<-ncol(X12)
  rk<-min(n1,n2,m1,m2)
  
  impute.X11=impute.X12=impute.X21=impute.X22=0
  
  na.ind.X11<- which(is.na(X11),arr.ind = T)
  if (length(na.ind.X11)>0){
    impute.X11<-rep(NA, nrow(na.ind.X11))
    for (j in 1:nrow(na.ind.X11)){
      colmean=mean(X11[,na.ind.X11[j,2]], na.rm=T)
      rowmean=mean(X11[na.ind.X11[j,1],], na.rm=T)
      impute.X11[j]=mean(c(colmean,rowmean), na.rm=T)
    }
    X11[na.ind.X11]=impute.X11
  }
  
  na.ind.X12<- which(is.na(X12),arr.ind = T)
  if (length(na.ind.X12)>0){
    impute.X12<-rep(NA, nrow(na.ind.X12))
    for (j in 1:nrow(na.ind.X12)){
      colmean=mean(X12[,na.ind.X12[j,2]], na.rm=T)
      rowmean=mean(X12[na.ind.X12[j,1],], na.rm=T)
      impute.X12[j]=mean(c(colmean,rowmean), na.rm=T)
    }
    X12[na.ind.X12]=impute.X12
  }
  
  na.ind.X21<- which(is.na(X21),arr.ind = T)
  if (length(na.ind.X21)>0){
    impute.X21<-rep(NA, nrow(na.ind.X21))
    for (j in 1:nrow(na.ind.X21)){
      colmean=mean(X21[,na.ind.X21[j,2]], na.rm=T)
      rowmean=mean(X21[na.ind.X21[j,1],], na.rm=T)
      impute.X21[j]=mean(c(colmean,rowmean), na.rm=T)
    }
    X21[na.ind.X21]=impute.X21
  }
  
  na.ind.X22<- which(is.na(X22),arr.ind = T)
  if (length(na.ind.X22)>0){
    impute.X22<-rep(NA, nrow(na.ind.X22))
    for (j in 1:nrow(na.ind.X22)){
      colmean=mean(X22[,na.ind.X22[j,2]], na.rm=T)
      rowmean=mean(X22[na.ind.X22[j,1],], na.rm=T)
      impute.X22[j]=mean(c(colmean,rowmean), na.rm=T)
    }
    X22[na.ind.X22]=impute.X22
  }
  
  if (is.null(noise.est)){
    if (is.null(sigma1) || is.null(sigma2)){
      sigma11=sd(X11); sigma12=sd(X12); 
      sigma21=sd(X21); sigma22=sd(X22)
    }
    else{ sigma11=sigma1[1]; sigma12=sigma1[2]; sigma21=sigma2[1]; sigma22=sigma2[2] }
  }
  else if (noise.est=="RMT"){
    sigma11=findNoiseSD.RMT(X11); sigma12=findNoiseSD.RMT(X12); 
    sigma21=findNoiseSD.RMT(X21); sigma22=findNoiseSD.RMT(X22);
  } 
  
  # Iterative algorithm to impute missing values
  bool2<-TRUE
  crit<-NULL
  count2<-1
  
  if (pbar==TRUE) pb <- txtProgressBar(min = 0, max=max.iter.impute, initial=0, char="-", style = 3)
  while (bool2){
    if (pbar==TRUE){  setTxtProgressBar(pb, count2)  }
    impute.X11.old=impute.X11
    impute.X12.old=impute.X12
    impute.X21.old=impute.X21
    impute.X22.old=impute.X22
    
    fit<-BIDIFAC(X1 = list(X11,X12), X2 = list(X21, X22),noise.est = NULL, 
                  sigma1 = c(sigma11, sigma12), sigma2=c(sigma21, sigma22),
                  U.start=U.start, V.start=V.start,
                  pbar = FALSE, ...)
    U.start=fit$U.out
    V.start=fit$V.out
    
    impute.X11=fit$X1.pred[[1]][na.ind.X11]
    impute.X12=fit$X1.pred[[2]][na.ind.X12]
    impute.X21=fit$X2.pred[[1]][na.ind.X21]
    impute.X22=fit$X2.pred[[2]][na.ind.X22]
    
    crit2=(frob(impute.X11-impute.X11.old)+frob(impute.X12-impute.X12.old)+frob(impute.X21-impute.X21.old)+frob(impute.X22-impute.X22.old))/
      (nrow(na.ind.X11)+nrow(na.ind.X12)+nrow(na.ind.X21)+nrow(na.ind.X22))
    if (crit2<eps.impute){ bool2=FALSE }
    else if (count2==max.iter.impute){ bool2=FALSE }
    else{ 
      count2=count2+1 
      X11[na.ind.X11]=impute.X11
      X12[na.ind.X12]=impute.X12
      X21[na.ind.X21]=impute.X21
      X22[na.ind.X22]=impute.X22
    }
  }
  
  fit$X1.original=X1.original
  fit$X2.original=X2.original
  
  return(list(X1.imputed=list(X11, X12),X2.imputed=list(X21,X22), fit=fit))
}




rearrange <- function(fit){
  ####Get row/column orderings
  Mat_RowOrder1 = cbind(fit$X1.global[[1]]+fit$X1.row[[1]], fit$X1.global[[2]]+fit$X1.row[[2]]) 
  Mat_RowOrder2 = cbind(fit$X2.global[[1]]+fit$X2.row[[1]], fit$X2.global[[2]]+fit$X2.row[[2]]) 
  
  Mat_ColOrder1 = rbind(fit$X1.global[[1]]+fit$X1.col[[1]], fit$X2.global[[1]]+fit$X2.col[[1]]) 
  Mat_ColOrder2 = rbind(fit$X1.global[[2]]+fit$X1.col[[2]], fit$X2.global[[2]]+fit$X2.col[[2]]) 
  
  row.orders = col.orders = list()
  
  col.orders[[1]] = hclust(dist(t(Mat_ColOrder1)))$order
  col.orders[[2]] = hclust(dist(t(Mat_ColOrder2)))$order
  
  row.orders[[1]] = hclust(dist(Mat_RowOrder1))$order
  row.orders[[2]] = hclust(dist(Mat_RowOrder2))$order
  
  Image_Data1 = Image_Joint1 = Image_Row1 = Image_Col1 = Image_Individual1 =  Image_Noise1 = list()
  Image_Data2 = Image_Joint2 = Image_Row2 = Image_Col2 = Image_Individual2 =  Image_Noise2 = list()
  
  Image_Data1[[1]] = as.matrix(fit$X1.original[[1]][row.orders[[1]],col.orders[[1]]]) 
  Image_Data1[[2]] = as.matrix(fit$X1.original[[2]][row.orders[[1]],col.orders[[2]]]) 
  Image_Data2[[1]] = as.matrix(fit$X2.original[[1]][row.orders[[2]],col.orders[[1]]]) 
  Image_Data2[[2]] = as.matrix(fit$X2.original[[2]][row.orders[[2]],col.orders[[2]]]) 
  
  Image_Joint1[[1]] = as.matrix(fit$X1.global[[1]][row.orders[[1]],col.orders[[1]]]) 
  Image_Joint1[[2]] = as.matrix(fit$X1.global[[2]][row.orders[[1]],col.orders[[2]]]) 
  Image_Joint2[[1]] = as.matrix(fit$X2.global[[1]][row.orders[[2]],col.orders[[1]]]) 
  Image_Joint2[[2]] = as.matrix(fit$X2.global[[2]][row.orders[[2]],col.orders[[2]]]) 
  
  Image_Col1[[1]] = as.matrix(fit$X1.col[[1]][row.orders[[1]],col.orders[[1]]]) 
  Image_Col1[[2]] = as.matrix(fit$X1.col[[2]][row.orders[[1]],col.orders[[2]]]) 
  Image_Col2[[1]] = as.matrix(fit$X2.col[[1]][row.orders[[2]],col.orders[[1]]]) 
  Image_Col2[[2]] = as.matrix(fit$X2.col[[2]][row.orders[[2]],col.orders[[2]]]) 
  
  Image_Row1[[1]] = as.matrix(fit$X1.row[[1]][row.orders[[1]],col.orders[[1]]]) 
  Image_Row1[[2]] = as.matrix(fit$X1.row[[2]][row.orders[[1]],col.orders[[2]]]) 
  Image_Row2[[1]] = as.matrix(fit$X2.row[[1]][row.orders[[2]],col.orders[[1]]]) 
  Image_Row2[[2]] = as.matrix(fit$X2.row[[2]][row.orders[[2]],col.orders[[2]]]) 
  
  Image_Individual1[[1]] = as.matrix(fit$X1.individual[[1]][row.orders[[1]],col.orders[[1]]]) 
  Image_Individual1[[2]] = as.matrix(fit$X1.individual[[2]][row.orders[[1]],col.orders[[2]]]) 
  Image_Individual2[[1]] = as.matrix(fit$X2.individual[[1]][row.orders[[2]],col.orders[[1]]]) 
  Image_Individual2[[2]] = as.matrix(fit$X2.individual[[2]][row.orders[[2]],col.orders[[2]]]) 
  
  Joint=rbind(cbind(Image_Joint1[[1]], Image_Joint1[[2]]), cbind(Image_Joint2[[1]], Image_Joint2[[2]]))
  Col=rbind(cbind(Image_Col1[[1]], Image_Col1[[2]]), cbind(Image_Col2[[1]], Image_Col2[[2]]))
  Row=rbind(cbind(Image_Row1[[1]], Image_Row1[[2]]), cbind(Image_Row2[[1]], Image_Row2[[2]]))
  Individual=rbind(cbind(Image_Individual1[[1]], Image_Individual1[[2]]), cbind(Image_Individual2[[1]], Image_Individual2[[2]]))
  
  
  return(list(
    Global=Joint,
    Col=Col,
    Row=Row,
    Individual=Individual,
    X1.joint=Image_Joint1,
    X2.joint=Image_Joint2,
    X1.col=Image_Col1,
    X2.col=Image_Col2,
    X1.row=Image_Row1,
    X2.row=Image_Row2,
    X1.individual=Image_Individual1,
    X2.individual=Image_Individual2
  ))
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


BIDsim = function(n=100, m=100, rkG=1, rkC=c(1,1), rkR=c(1,1), rkI=c(1,1,1,1), SNR=1){
  # BIDsim simulates 2 by 2 linkced matrices with Gaussian noise.
  # Built based on the LRsim() of the denoiseR package.
  # Assumes that the sizes of the matrices are all the same (nrow=n, ncol=m).
  # the joint, column-shared, row-shared, and individual matrices are assumed to be pairwise orthogonal
  # SNR is defined by 1/sigma/sqrt(n*m), as each "signal" has frobenious norm 1.
  
  #INPUT
  # n,m  : the number of rows and columns in each matrix
  # rkG  : rank of globally shared components
  # rkR  : ranks of row-shared matrices
  #      : a vector of length 2. If constant, it will assume the same rank.
  # rkC  : ranks of column-shared matrices
  #      : a vector of length 2. If constant, it will assume the same rank.
  # rkI  : ranks of individual matrices X11,X12,X21,X22
  #      : a vector of length 4. If constant, it will assume the same rank.
  # SNR  : Singal-to noise ratio for each matrix X11,X12,X21,X22
  #      : a vector of length 4
  #      : if length is 1, BIDsim assumes the equal SNR to all matrices
  
  #OUTPUT
  #X     : generated matrices with Gaussian noise
  #      : a list of length 4
  #S     : signal components of X
  #G     : global components of X
  #R     : row-shared components of X
  #C     : column-shared components of X
  #I     : individual components of X
  #sigma : error standard deviation of each matrix
  #      : a vector of length 4
  
  if (length(SNR)==1){ SNR= rep(SNR,4) }
  if (length(rkR)==1){ rkR=rep(rkR,2) }
  if (length(rkC)==1){ rkC=rep(rkC,2) }
  if (length(rkI)==1){ rkI=rep(rkI,4) }
  
  rk=2*rkG+2*sum(rkC)+2*sum(rkR)+sum(rkI)
  rkC1=rkC[1]; rkC2=rkC[2]
  rkR1=rkR[1]; rkR2=rkR[2]
  rkI11=rkI[1]; rkI12=rkI[2]; rkI21=rkI[3]; rkI22=rkI[4]
  
  SNR11=SNR[1]; SNR12=SNR[2]; SNR21=SNR[3]; SNR22=SNR[4]
  
  rk11=rkG+rkC1+rkR1+rkI11
  rk12=rkG+rkC2+rkR1+rkI12
  rk21=rkG+rkC1+rkR2+rkI21
  rk22=rkG+rkC2+rkR2+rkI22
  
  noise.sd.11 = 1/(SNR11*sqrt(n*m)) 
  noise.sd.12 = 1/(SNR12*sqrt(n*m))
  noise.sd.21 = 1/(SNR21*sqrt(n*m))
  noise.sd.22 = 1/(SNR22*sqrt(n*m))
  
  csumC=cumsum(c(rep(rkG,2),rep(rkC,2),rkR,rkI))
  csumC1=csumC+1
  csumR=cumsum(c(rep(rkG,2),rkC,rep(rkR,2),rkI))
  csumR1=csumR+1
  
  index.u=sample(min(c(m,n)),rk)
  index.v=sample(min(c(m,n)),rk)
  
  if (rkG>0){
    index.G.u10=index.u[1:csumR[1]]
    index.G.u20=index.u[csumR1[1]:csumR[2]]
    index.G.v01=index.v[1:csumC[1]]
    index.G.v02=index.v[csumC1[1]:csumC[2]]
  }
  else{
    index.G.u10=index.G.u20=index.G.v01=index.G.v02=NULL
  }
  
  if (csumC1[2]-csumC[3]<=0 ){ index.C.u11=index.u[csumC1[2]:csumC[3]] } else{ index.C.u11=NULL }
  if (csumC1[3]-csumC[4]<=0 ){ index.C.u12=index.u[csumC1[3]:csumC[4]] } else{ index.C.u12=NULL }
  if (csumC1[4]-csumC[5]<=0 ){ index.C.u21=index.u[csumC1[4]:csumC[5]] } else{ index.C.u21=NULL }
  if (csumC1[5]-csumC[6]<=0 ){ index.C.u22=index.u[csumC1[5]:csumC[6]] } else{ index.C.u22=NULL }
  if (csumC1[6]-csumC[7]<=0 ){ index.R.u10=index.u[csumC1[6]:csumC[7]] } else{ index.R.u10=NULL }
  if (csumC1[7]-csumC[8]<=0 ){ index.R.u20=index.u[csumC1[7]:csumC[8]] } else{ index.R.u20=NULL }
  if (csumC1[8]-csumC[9]<=0 ){ index.I.u11=index.u[csumC1[8]:csumC[9]] } else{ index.I.u11=NULL }
  if (csumC1[9]-csumC[10]<=0 ){ index.I.u12=index.u[csumC1[9]:csumC[10]] } else{ index.I.u12=NULL }
  if (csumC1[10]-csumC[11]<=0 ){ index.I.u21=index.u[csumC1[10]:csumC[11]] } else{ index.I.u21=NULL }
  if (csumC1[11]-csumC[12]<=0 ){ index.I.u22=index.u[csumC1[11]:csumC[12]] } else{ index.I.u22=NULL }
  
  if (csumR1[2]-csumR[3]<=0 ){ index.C.v01=index.v[csumR1[2]:csumR[3]] } else{ index.C.v01=NULL }
  if (csumR1[3]-csumR[4]<=0 ){ index.C.v02=index.v[csumR1[2]:csumR[3]] } else{ index.C.v02=NULL }
  if (csumR1[4]-csumR[5]<=0 ){ index.R.v11=index.v[csumR1[6]:csumR[7]] } else{ index.R.v11=NULL }
  if (csumR1[5]-csumR[6]<=0 ){ index.R.v21=index.v[csumR1[7]:csumR[8]] } else{ index.R.v21=NULL }
  if (csumR1[6]-csumR[7]<=0 ){ index.R.v12=index.v[csumR1[6]:csumR[7]] } else{ index.R.v12=NULL }
  if (csumR1[7]-csumR[8]<=0 ){ index.R.v22=index.v[csumR1[7]:csumR[8]] } else{ index.R.v22=NULL }
  if (csumR1[8]-csumR[9]<=0 ){ index.I.v11=index.v[csumR1[8]:csumR[9]] } else{ index.I.v11=NULL }
  if (csumR1[9]-csumR[10]<=0 ){ index.I.v12=index.v[csumR1[9]:csumR[10]] } else{ index.I.v12=NULL }
  if (csumR1[10]-csumR[11]<=0 ){ index.I.v21=index.v[csumR1[10]:csumR[11]] } else{ index.I.v21=NULL }
  if (csumR1[11]-csumR[12]<=0 ){ index.I.v22=index.v[csumR1[11]:csumR[12]] } else{ index.I.v22=NULL }
  
  if(rk == 0){
    MU <- matrix(0, n, m)  
  } else {
    signal <- replicate(m, rnorm(n, 0, 1))
    signal <- scale(signal, scale = FALSE)
    svd.signal <- svd(signal) 
    index.u.X11=c(index.G.u10,index.C.u11,index.R.u10,index.I.u11);
    index.u.X12=c(index.G.u10,index.C.u12,index.R.u10,index.I.u12);
    index.u.X21=c(index.G.u20,index.C.u21,index.R.u20,index.I.u21);
    index.u.X22=c(index.G.u20,index.C.u22,index.R.u20,index.I.u22);
    
    index.v.X11=c(index.G.v01,index.C.v01,index.R.v11,index.I.v11);
    index.v.X12=c(index.G.v01,index.C.v02,index.R.v12,index.I.v12);
    index.v.X21=c(index.G.v02,index.C.v01,index.R.v21,index.I.v21);
    index.v.X22=c(index.G.v02,index.C.v02,index.R.v22,index.I.v22);
    
    samp.G.X11=sample(rk, rkG)
    samp.G.X12=sample(rk, rkG)
    samp.G.X21=sample(rk, rkG)
    samp.G.X22=sample(rk, rkG)
    
    samp.C.X11=sample(rk, rkC1)
    samp.C.X12=sample(rk, rkC2)
    samp.C.X21=sample(rk, rkC1)
    samp.C.X22=sample(rk, rkC2)
    
    samp.R.X11=sample(rk, rkR1)
    samp.R.X12=sample(rk, rkR1)
    samp.R.X21=sample(rk, rkR2)
    samp.R.X22=sample(rk, rkR2)
    
    samp.I.X11=sample(rk, rkI11)
    samp.I.X12=sample(rk, rkI12)
    samp.I.X21=sample(rk, rkI21)
    samp.I.X22=sample(rk, rkI22)
    
    samp.X11=c(samp.G.X11,samp.C.X11,samp.R.X11, samp.I.X11)
    samp.X12=c(samp.G.X12,samp.C.X12,samp.R.X11, samp.I.X12)
    samp.X21=c(samp.G.X21,samp.C.X21,samp.R.X11, samp.I.X21)
    samp.X22=c(samp.G.X22,samp.C.X22,samp.R.X11, samp.I.X22)
    
    G11=svd.signal$u[, index.G.u10,drop=F] %*% diag(svd.signal$d[samp.G.X11], rkG, rkG) %*% t(svd.signal$v[, index.G.v01, drop=F])/ sqrt(sum(svd.signal$d[samp.X11]^2))
    G12=svd.signal$u[, index.G.u10,drop=F] %*% diag(svd.signal$d[samp.G.X12], rkG, rkG) %*% t(svd.signal$v[, index.G.v02, drop=F])/ sqrt(sum(svd.signal$d[samp.X12]^2))
    G21=svd.signal$u[, index.G.u20,drop=F] %*% diag(svd.signal$d[samp.G.X21], rkG, rkG) %*% t(svd.signal$v[, index.G.v01, drop=F])/ sqrt(sum(svd.signal$d[samp.X21]^2))
    G22=svd.signal$u[, index.G.u20,drop=F] %*% diag(svd.signal$d[samp.G.X22], rkG, rkG) %*% t(svd.signal$v[, index.G.v02, drop=F])/ sqrt(sum(svd.signal$d[samp.X22]^2))
    
    R11=svd.signal$u[, index.R.u10,drop=F] %*% diag(svd.signal$d[samp.R.X11], rkR1, rkR1) %*% t(svd.signal$v[, index.R.v11, drop=F])/ sqrt(sum(svd.signal$d[samp.X11]^2))
    R12=svd.signal$u[, index.R.u10,drop=F] %*% diag(svd.signal$d[samp.R.X12], rkR1, rkR1) %*% t(svd.signal$v[, index.R.v12, drop=F])/ sqrt(sum(svd.signal$d[samp.X12]^2))
    R21=svd.signal$u[, index.R.u20,drop=F] %*% diag(svd.signal$d[samp.R.X21], rkR2, rkR2) %*% t(svd.signal$v[, index.R.v21, drop=F])/ sqrt(sum(svd.signal$d[samp.X21]^2))
    R22=svd.signal$u[, index.R.u20,drop=F] %*% diag(svd.signal$d[samp.R.X22], rkR2, rkR2) %*% t(svd.signal$v[, index.R.v22, drop=F])/ sqrt(sum(svd.signal$d[samp.X22]^2))
    
    C11=svd.signal$u[, index.C.u11,drop=F] %*% diag(svd.signal$d[samp.C.X11], rkC1, rkC1) %*% t(svd.signal$v[, index.C.v01, drop=F])/ sqrt(sum(svd.signal$d[samp.X11]^2))
    C12=svd.signal$u[, index.C.u12,drop=F] %*% diag(svd.signal$d[samp.C.X12], rkC2, rkC2) %*% t(svd.signal$v[, index.C.v02, drop=F])/ sqrt(sum(svd.signal$d[samp.X12]^2))
    C21=svd.signal$u[, index.C.u21,drop=F] %*% diag(svd.signal$d[samp.C.X21], rkC1, rkC1) %*% t(svd.signal$v[, index.C.v01, drop=F])/ sqrt(sum(svd.signal$d[samp.X21]^2))
    C22=svd.signal$u[, index.C.u22,drop=F] %*% diag(svd.signal$d[samp.C.X22], rkC2, rkC2) %*% t(svd.signal$v[, index.C.v02, drop=F])/ sqrt(sum(svd.signal$d[samp.X22]^2))
    
    I11=svd.signal$u[, index.I.u11,drop=F] %*% diag(svd.signal$d[samp.I.X11], rkI11, rkI11) %*% t(svd.signal$v[, index.I.v11, drop=F])/ sqrt(sum(svd.signal$d[samp.X11]^2))
    I12=svd.signal$u[, index.I.u12,drop=F] %*% diag(svd.signal$d[samp.I.X12], rkI12, rkI12) %*% t(svd.signal$v[, index.I.v12, drop=F])/ sqrt(sum(svd.signal$d[samp.X12]^2))
    I21=svd.signal$u[, index.I.u21,drop=F] %*% diag(svd.signal$d[samp.I.X21], rkI21, rkI21) %*% t(svd.signal$v[, index.I.v21, drop=F])/ sqrt(sum(svd.signal$d[samp.X21]^2))
    I22=svd.signal$u[, index.I.u22,drop=F] %*% diag(svd.signal$d[samp.I.X22], rkI22, rkI22) %*% t(svd.signal$v[, index.I.v22, drop=F])/ sqrt(sum(svd.signal$d[samp.X22]^2))
    
    S11 = G11+C11+R11+I11 
    S12 = G12+C12+R12+I12 
    S21 = G21+C21+R21+I21 
    S22 = G22+C22+R22+I22
  }
  
  X11 <- S11 + noise.sd.11*replicate(m, rnorm(n, 0, 1))
  X12 <- S12 + noise.sd.12*replicate(m, rnorm(n, 0, 1))
  X21 <- S21 + noise.sd.21*replicate(m, rnorm(n, 0, 1))
  X22 <- S22 + noise.sd.22*replicate(m, rnorm(n, 0, 1))
  return(list(X = list(X11,X12,X21,X22), S = list(S11,S12,S21,S22), 
              G = list(G11,G12,G21,G22), C= list(C11,C12,C21,C22),
              R = list(R11,R12,R21,R22), I= list(I11,I12,I21,I22),
              sigma = c(noise.sd.11,noise.sd.12,noise.sd.21,noise.sd.22)))
}


generateRank=function(rk=10, on.list=rep(1,4), minrk.list=rep(1,4)){
  # onlist = a vector of length 4, where each element is either 0 or 1.
  #          the first element=0 enforces rkg to be 0
  #          the second element=0 enforces rkR to be 0  
  #          the third element=0 enforces rkC to be 0
  #          the fourth element=0 enforces rkI to be 0
  
  # minrk.list = a vector of length 4, where each element is either 0 or 1.
  
  minrk.list1=minrk.list[on.list==1]
  x = rmultinom(n = 1, size = rk-sum(minrk.list1), prob = on.list/sum(on.list))
  x = x+minrk.list*on.list
  
  return(list(rkg=x[1], rkR=x[2], rkC=x[3], rkI=x[4]))
}


