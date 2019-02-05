###################################################
##Simulation studies using BIDIFAC() (new)
##Author: Jun Young Park
##Contact: park1131@umn.edu
###################################################

set.seed(1131)

#Generate simulated data
simData=BIDsim(m=100, n=100, rkJ = 2, rkC=2, rkR=2, rkI=3, SNR=1)

###################################################
##Simulation studies using BIDIFAC()
##Author: Jun Young Park
##Contact: park1131@umn.edu
###################################################

set.seed(1131)

###################################################
##########STEP 0: Generate simulated data##########
###################################################

#When ranks are determined randomly
# getrk=generateRank(rk=10)
# simData=BIDsim(m=100, n=100, rkJ=getrk$rkJ, rkC=getrk$rkC, rkR=getrk$rkR, rkI = getrk$rkI, SNR=2)

#When ranks are specified
simData=BIDsim(m=100, n=100, rkJ = 3, rkC=2, rkR=2, rkI=3, SNR=2)


###################################################
#########STEP 1: Run BIDIFAC for the data##########
###################################################

fit=BIDIFAC(list(simData$X[[1]], simData$X[[2]]), list(simData$X[[3]], simData$X[[4]]))

#Summarize the error terms
tbl=cbind(
  c( frob(fit$X1.joint[[1]]-simData$J[[1]])/frob(simData$J[[1]]),
     frob(fit$X1.joint[[2]]-simData$J[[2]])/frob(simData$J[[2]]),
     frob(fit$X2.joint[[1]]-simData$J[[3]])/frob(simData$J[[3]]),
     frob(fit$X2.joint[[2]]-simData$J[[4]])/frob(simData$J[[4]])),
  c( frob(fit$X1.col[[1]]-simData$C[[1]])/frob(simData$C[[1]]),
     frob(fit$X1.col[[2]]-simData$C[[2]])/frob(simData$C[[2]]),
     frob(fit$X2.col[[1]]-simData$C[[3]])/frob(simData$C[[3]]),
     frob(fit$X2.col[[2]]-simData$C[[4]])/frob(simData$C[[4]])),
  c( frob(fit$X1.row[[1]]-simData$R[[1]])/frob(simData$R[[1]]),
     frob(fit$X1.row[[2]]-simData$R[[2]])/frob(simData$R[[2]]),
     frob(fit$X2.row[[1]]-simData$R[[3]])/frob(simData$R[[3]]),
     frob(fit$X2.row[[2]]-simData$R[[4]])/frob(simData$R[[4]])),
  c( frob(fit$X1.individual[[1]]-simData$I[[1]])/frob(simData$I[[1]]),
     frob(fit$X1.individual[[2]]-simData$I[[2]])/frob(simData$I[[2]]),
     frob(fit$X2.individual[[1]]-simData$I[[3]])/frob(simData$I[[3]]),
     frob(fit$X2.individual[[2]]-simData$I[[4]])/frob(simData$I[[4]]))
)

colnames(tbl)=c("Global Error", "Column Error", "Row Error", "Individual Error")
rownames(tbl)=c("X11", "X12", "X21", "X22")

tbl
#     Global Error Column Error  Row Error Individual Error
# X11   0.08773571   0.08684621 0.07185276        0.1311311
# X12   0.08282229   0.11364179 0.07535306        0.1372327
# X21   0.07829041   0.11267822 0.12808729        0.2158613
# X22   0.15541707   0.10217324 0.08440206        0.1582459


###################################################
#########STEP 2: Visualize the BIDIFAC fit#########
###################################################

refit=rearrange(fit)
show.image(refit$Joint)
show.image(refit$Col)
show.image(refit$Row)
show.image(refit$Individual)


###################################################
##########STEP 3: Missing data imputation##########
###################################################

#Missing column case
simData$X[[1]][,c(10,20)]=NA
simData$X[[2]][,c(30,40)]=NA
simData$X[[3]][,c(50,60)]=NA
simData$X[[4]][,c(70,80)]=NA

fit<-impute.BIDIFAC(list(simData$X[[1]],simData$X[[2]]),list(simData$X[[3]],simData$X[[4]]), eps.impute=1e-10, max.iter.impute=10)
tbl.col=c(frob(fit$X1.imputed[[1]][,c(10,20)]-simData$S[[1]][,c(10,20)])/frob(simData$S[[1]][,c(10,20)]),
          frob(fit$X1.imputed[[2]][,c(30,40)]-simData$S[[2]][,c(30,40)])/frob(simData$S[[2]][,c(30,40)]),
          frob(fit$X2.imputed[[1]][,c(50,60)]-simData$S[[3]][,c(50,60)])/frob(simData$S[[3]][,c(50,60)]),
          frob(fit$X2.imputed[[2]][,c(70,80)]-simData$S[[4]][,c(70,80)])/frob(simData$S[[4]][,c(70,80)]))
names(tbl.col)=c("X11","X12","X21","X22")
tbl.col
#       X11       X12       X21       X22 
# 00.7131981 0.6048519 0.4965431 0.7093090 
 
