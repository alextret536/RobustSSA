library(Rssa)
library(pcaL1)
library(matrixcalc)
set.seed(77)
M <- 10
#M<-30
L<-120

MSE.SSA <- rep(0,M)
MSE.L1svd <- rep(0,M)
MSE.loess <- rep(0,M)
MSE.med <- rep(0,M)
MSE.lowess <- rep(0,M)
MSE.or <- rep(0,M)

MAD.SSA <- rep(0,M)
MAD.L1svd <- rep(0,M)
MAD.loess <- rep(0,M)
MAD.med <- rep(0,M)
MAD.lowess <- rep(0,M)
MAD.or <- rep(0,M)


for(k in 1:M)
{
  sig <- exp((1:N)/N)+sin(2*pi*(1:N)/120+pi/6) # 1 example
  rnk<-3
  
  #sig <- exp(4*(1:N)/N)*sin(2*pi*(1:N)/30) # 2 example
  #rnk<-2
  
  #sig <- (1:N)*exp(4*(1:N)/N)*sin(2*pi*(1:N)/30) # 3 example
  #rnk<-4
  
  sig.outl<-sig
  outlier.seq<-sample(1:(N),N*0.05)
  sig.outl[outlier.seq]<-sig.outl[outlier.seq]+5*sig.outl[outlier.seq]
  #ser<-sig.outl+0.5*exp(4*(1:N)/N)*rnorm(N) #for 2 example
  ser<-sig.outl+rnorm(N)
  
  #plot(ser,type='l')
  X<-hankel(ser, L=L)
  
  s <- ssa(ser, L = L) # Basic SSA
  rec <- reconstruct(s, groups = list(c(1:rnk)))
  trend.season <- rec$F1
  
  X<-hankel(ser,L=L)
  s.L1svd<-l1pca(X,center=FALSE,projections="l1",projDim=rnk, iterations=5) #l1pca
  Pr<-s.L1svd$projPoints
  Pr.L1svd<-hankL1(Pr)
  
  Pr<-IRLS_mod(X,rnk,'loess') #IRLS modification with loess
  Pr0<-hankL2(Pr)
  
  Pr<-IRLS_mod(X,rnk,'median') #IRLS modification with median
  Pr1<-hankL2(Pr)
  
  Pr<-IRLS_mod(X,rnk,'lowess') #IRLS modification with lowess
  Pr2<-hankL2(Pr)
  
  Pr<-IRLS_orig(X,rnk) #original IRLS
  Pr3<-hankL2(Pr)
  
  
  #MSE
  MSE.SSA[k] <- mean((sig - trend.season)[1:N]^2)
  MSE.L1svd[k] <- mean((sig - Pr.L1svd)[1:N]^2)
  MSE.loess[k] <- mean((sig - Pr0)[1:N]^2)
  MSE.med[k] <- mean((sig - Pr1)[1:N]^2)
  MSE.lowess[k] <- mean((sig - Pr2)[1:N]^2)
  MSE.or[k] <- mean((sig - Pr3)[1:N]^2)
  
  #MAD
  MAD.SSA[k] <- mean(abs((sig - trend.season)[1:N]))
  MAD.L1svd[k] <- mean(abs((sig - Pr.L1svd)[1:N]))
  MAD.loess[k] <- mean(abs((sig - Pr0)[1:N]))
  MAD.med[k] <- mean(abs((sig - Pr1)[1:N]))
  MAD.lowess[k] <- mean(abs((sig - Pr2)[1:N]))
  MAD.or[k] <- mean(abs((sig - Pr3)[1:N]))
  
}

RMSE.SSA<-sqrt(mean(MSE.SSA))
RMSE.L1svd<-sqrt(mean(MSE.L1svd))
RMSE.loess<-sqrt(mean(MSE.loess))
RMSE.med<-sqrt(mean(MSE.med))
RMSE.lowess<-sqrt(mean(MSE.lowess))
RMSE.or<-sqrt(mean(MSE.or))

MAD.SSA<-mean(MAD.SSA)
MAD.L1svd<-mean(MAD.L1svd)
MAD.loess<-mean(MAD.loess)
MAD.med<-mean(MAD.med)
MAD.lowess<-mean(MAD.lowess)
MAD.or<-mean(MAD.or)

#RMSE
RMSE.SSA #error SSA
RMSE.L1svd #error L1svd
RMSE.loess #error loess
RMSE.med #error median
RMSE.lowess #error lowess
RMSE.or #error original IRLS

#MAD
MAD.SSA #error SSA
MAD.L1svd #error L1svd
MAD.loess #error loess
MAD.med #error median
MAD.lowess #error lowess
MAD.or #error original IRLS

#p-values
pval<-t.test(MSE.loess,MSE.L1svd,paired=TRUE,alternative="two.sided")$p.value
print(pval)

pval<-t.test(MSE.loess,MSE.med,paired=TRUE,alternative="two.sided")$p.value
print(pval)

pval<-t.test(MSE.SSA,MSE.loess,paired=TRUE,alternative="two.sided")$p.value
print(pval)

pval<-t.test(MSE.loess,MSE.or,paired=TRUE,alternative="two.sided")$p.value
print(pval)

pval<-t.test(MSE.loess,MSE.lowess,paired=TRUE,alternative="two.sided")$p.value
print(pval)
