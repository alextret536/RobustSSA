library(Rssa)
library(matrixcalc)
library(pcaL1)
df <- read.csv(file = 'IMP5130.csv')
head(df)
series<-df$IMP5130[0:240]
series[95]<-series[95]*3
series[234]<-series[234]*2
plot(series,type='l')
length(series)
L<-60
s <- ssa(series, L = L)

plot(wcor(s))
X<-hankel(series, L=L)
rnk<-5

Pr<-IRLS_mod(X,rnk,'loess') #IRLS modification (trend extraction with loess)
Pr0<-hankL2(Pr)

#Pr<-IRLS_mod(X,rnk,'median') #IRLS modification (trend extraction with median)
#Pr1<-hankL2(Pr)

Pr<-IRLS_mod(X,rnk,'lowess') #IRLS modification (trend extraction with lowess)
Pr2<-hankL2(Pr)

Pr<-IRLS_orig(X,rnk) #IRLS original
Pr3<-hankL2(Pr)

#s.L1svd<-l1pca(X,center=FALSE,projections="l1",projDim=rnk) #l1pca
#Pr<-s.L1svd$projPoints
#Pr.L1<-hankL1(Pr)

rec <- reconstruct(s, groups = list(c(1:rnk))) 
trend.season <- rec$F1


plot(series,col='black',type='l')
lines(trend.season,type='l',col='red',lw=2)
lines( Pr2,type='l',col='orange',lw=2)
lines( Pr3,type='l',col='green',lw=2)
lines(Pr0,type='l',col='blue',lw=2)
legend('topleft', c("SSA","IRLS loess","IRLS lowess","IRLS orig"),
       col=c("red","blue","orange","green"), lty=1, cex=0.8, lw=c(2,2,2,2))
