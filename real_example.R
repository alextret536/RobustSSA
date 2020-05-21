library(Rssa)
library(matrixcalc)
df <- read.csv(file = 'FRAPROINDMISMEI.csv')
head(df)
series<-df$FRAPROINDMISMEI[0:528]
plot(series,type='l')
length(series)
L=48
s <- ssa(series, L = L)

plot(wcor(s))
ser<-series
X<-hankel(ser, L=L)
rnk<-8

Pr<-IRLS_mod(X,rnk,'loess')
Pr0<-hankL2(Pr)

Pr<-IRLS_mod(X,rnk,'median')
Pr1<-hankL2(Pr)

Pr<-IRLS_mod(X,rnk,'lowess')
Pr2<-hankL2(Pr)


rec <- reconstruct(s, groups = list(c(1:rnk)))
trend.season <- rec$F1


plot(ser,col='black',type='l')
lines(Pr0,type='l',col='violet',lw=2)
lines(trend.season,type='l',col='red',lw=2)
lines( Pr1,type='l',col='green',lw=2)
lines( Pr2,type='l',col='blue',lw=2)
legend('topleft', c("SSA","IRLS loess","IRLS med","IRLS lowess"),
       col=c("red","violet","green","blue"), lty=1, cex=0.8, lw=c(1,2,2,2))
