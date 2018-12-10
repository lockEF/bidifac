
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
simData=BIDsim(m=100, n=100, rkG = 3, rkC=2, rkR=2, rkI=3, SNR=2)


###################################################
#########STEP 1: Run BIDIFAC for the data##########
###################################################

fit=BIDIFAC(list(simData$X[[1]], simData$X[[2]]), list(simData$X[[3]], simData$X[[4]]))

#Summarize the error terms
tbl=cbind(
  c( frob(fit$X1.global[[1]]-simData$G[[1]])/frob(simData$G[[1]]),
     frob(fit$X1.global[[2]]-simData$G[[2]])/frob(simData$G[[2]]),
     frob(fit$X2.global[[1]]-simData$G[[3]])/frob(simData$G[[3]]),
     frob(fit$X2.global[[2]]-simData$G[[4]])/frob(simData$G[[4]])),
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
#    Global Error Column Error  Row Error Individual Error
#X11   0.07647018   0.08634417 0.07839756        0.1643988
#X12   0.06802652   0.08067275 0.08059401        0.1686063
#X21   0.06606632   0.08096843 0.08169410        0.2507125
#X22   0.11996727   0.08973148 0.07860975        0.1720027


###################################################
#########STEP 2: Visualize the BIDIFAC fit#########
###################################################

refit=rearrange(fit)
show.image(refit$Global)
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

fit<-impute.BIDIFAC(list(simData$X[[1]],simData$X[[2]]),list(simData$X[[3]],simData$X[[4]]), eps.impute=1e-20, max.iter.impute=10)
tbl.col=c(frob(fit$X1.imputed[[1]][,c(10,20)]-simData$S[[1]][,c(10,20)])/frob(simData$S[[1]][,c(10,20)]),
          frob(fit$X1.imputed[[2]][,c(30,40)]-simData$S[[2]][,c(30,40)])/frob(simData$S[[2]][,c(30,40)]),
          frob(fit$X2.imputed[[1]][,c(50,60)]-simData$S[[3]][,c(50,60)])/frob(simData$S[[3]][,c(50,60)]),
          frob(fit$X2.imputed[[2]][,c(70,80)]-simData$S[[4]][,c(70,80)])/frob(simData$S[[4]][,c(70,80)]))
names(tbl.col)=c("X11","X12","X21","X22")
tbl.col
#       X11       X12       X21       X22 
# 0.7043202 0.5875771 0.4798809 0.6981753 