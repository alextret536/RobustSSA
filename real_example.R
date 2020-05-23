library(Rssa)
library(matrixcalc)
library(pcaL1)
df <- read.csv(file = 'XTEXVA01CNM667S.csv')
head(df)
series<-df$XTEXVA01CNM667S
plot(series,type='l')
length(series)
L<-60
s <- ssa(series, L = L)

plot(wcor(s))
X<-hankel(series, L=L)
rnk<-8

Pr<-IRLS_mod(X,rnk,'loess')
Pr0<-hankL2(Pr)

Pr<-IRLS_mod(X,rnk,'median')
Pr1<-hankL2(Pr)

Pr<-IRLS_mod(X,rnk,'lowess')
Pr2<-hankL2(Pr)

Pr<-IRLS_orig(X,rnk)
Pr3<-hankL2(Pr)

s.L1svd<-l1pca(X,center=FALSE,projections="l1",projDim=rnk)
Pr<-s.L1svd$projPoints
Pr.L1<-hankL1(Pr)

rec <- reconstruct(s, groups = list(c(1:rnk)))
trend.season <- rec$F1


plot(series,col='black',type='l')
lines(trend.season,type='l',col='red',lw=2)
lines(Pr0,type='l',col='violet',lw=2)
lines( Pr1,type='l',col='green',lw=2)
lines( Pr2,type='l',col='orange',lw=2)
lines( Pr.L1,type='l',col='yellow',lw=2)
lines( Pr3,type='l',col='blue',lw=2)
legend('topleft', c("SSA","IRLS loess","IRLS med","IRLS lowess","IRLS orig","l1pca"),
       col=c("red","violet","green","orange","blue","yellow"), lty=1, cex=0.8, lw=c(2,2,2,2,2,2,2))

