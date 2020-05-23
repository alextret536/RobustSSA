library(Rssa)
library(matrixcalc)
library(pcaL1)
N<-240
Per<-120

set.seed(7)
sig <- exp(4*(1:N)/N)*sin(2*pi*(1:N)/30) 
rnk<-2
sig.outl<-sig
outlier.seq<-sample(1:(N),N*0.05)
sig.outl[outlier.seq]<-sig.outl[outlier.seq]+5*sig.outl[outlier.seq]
ser<-sig.outl+0.5*exp(4*(1:N)/N)*rnorm(N)

plot(ser,type='l')
X<-hankel(ser, L=120)

Pr<-IRLS_mod(X,rnk,'loess') #IRLS modification (trend extraction with loess)
Pr0<-hankL2(Pr)

Pr<-IRLS_mod(X,rnk,'median') #IRLS modification (trend extraction with median)
Pr1<-hankL2(Pr)

Pr<-IRLS_mod(X,rnk,'lowess') #IRLS modification (trend extraction with lowess)
Pr2<-hankL2(Pr)

Pr<-IRLS_orig(X,rnk) #IRLS original
Pr3<-hankL2(Pr)

s.L1svd<-l1pca(X,center=FALSE,projections="l1",projDim=rnk) #l1pca
Pr<-s.L1svd$projPoints
Pr.L1<-hankL1(Pr)


plot(ser,col='black',type='l')
lines( Pr0,type='l',col='violet',lw=2)
lines( Pr1,type='l',col='green',lw=2)
lines( Pr2,type='l',col='red',lw=2)
lines( Pr.L1,type='l',col='yellow',lw=2)
lines( Pr3,type='l',col='blue',lw=2)
legend('topleft', c("IRLS loess","IRLS med","IRLS lowess","IRLS orig","l1pca"),
       col=c("violet","green","red","blue","yellow"), lty=1, cex=0.8, lw=c(2,2,2,2,2))

#RMSE
sqrt(mean((sig - Pr0)^2)) #RMSE IRLS loess
sqrt(mean((sig - Pr1)^2)) #RMSE IRLS med
sqrt(mean((sig - Pr2)^2))#RMSE IRLS lowess
sqrt(mean((sig - Pr3)^2))#RMSE IRLS orig
sqrt(mean((sig - Pr.L1)^2))#RMSE IRLS l1pca
