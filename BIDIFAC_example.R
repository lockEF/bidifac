###################################################
##Simulation studies using BIDIFAC() (new)
##Author: Jun Young Park
##Contact: park1131@umn.edu
###################################################

#Example code
set.seed(1131)

#3 x 5 structure, mixed SNR
#Specifying dimensions of bidimensionally linked matrices
m.vec=c(100,150,200) 
n.vec=rep(100,5)

snr=matrix(sample(c(1/2,1, 2),15,replace = T),3,5)
simData=BIDsim(m.vec, n.vec,rkG=2, rkC=2, rkR=2, rkI=2,SNR=snr)

fit=BIDIFAC(simData$X)

# proportion of variance explained for the estimated components.
round(summary.BIDIFAC(fit),3)

#Concatenate estimate components.
Gmat=data.rearrange(fit$G)
Rmat=data.rearrange(fit$R)
Cmat=data.rearrange(fit$C)
Imat=data.rearrange(fit$I)
Smat=data.rearrange(fit$S)

#Visualization
show.image(Gmat$out)
show.image(Rmat$out)
show.image(Cmat$out)
show.image(Imat$out)
show.image(Smat$out)

#impute missing columns AND rows
simData1=simData
for (i in 1:3){
  for (j in 1:5){
    ind.row=sample(nrow(simData1$X[[i,j]]),5)
    ind.col=sample(ncol(simData1$X[[i,j]]),5)
    simData1$X[[i,j]][,ind.col]=NA
    simData1$X[[i,j]][ind.row,]=NA
  }
}

imp.fit=impute.BIDIFAC(simData1$X)
 
 
