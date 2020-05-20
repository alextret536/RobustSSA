library(Rssa)
library(matrixcalc)
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

Pr<-IRLS_mod(X,rnk,'loess')
Pr0<-hankL2(Pr)

Pr<-IRLS_mod(X,rnk,'median')
Pr1<-hankL2(Pr)

Pr<-IRLS_mod(X,rnk,'lowess')
Pr2<-hankL2(Pr)

Pr<-IRLS_orig(X,rnk)
Pr3<-hankL2(Pr)

plot(sig,type='l')
lines(ser,col='gray')
lines(Pr0,type='l',col='violet',lw=2)
lines( Pr1,type='l',col='green',lw=2)
lines( Pr2,type='l',col='red',lw=2)
lines( Pr3,type='l',col='orange',lw=2)
legend('topleft', c("signal","IRLS loess","IRLS med","IRLS lowess", "IRLS (orig.)"),
       col=c("black","violet","green","red","orange"), lty=1, cex=0.8, lw=c(1,2,2,2,2))

#RMSE
sqrt(mean((sig - Pr0)^2)) #RMSE IRLS loess
sqrt(mean((sig - Pr1)^2)) #RMSE IRLS med
sqrt(mean((sig - Pr2)^2))#RMSE IRLS lowess
sqrt(mean((sig - Pr3)^2))#RMSE IRLS orig
