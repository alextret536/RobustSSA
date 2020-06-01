library(Rssa)
library(matrixcalc)
library(pcaL1)
df <- read.csv(file = 'IMP5130.csv')
head(df)
cut <- 12
series<-df$IMP5130[cut:240]
series[95-cut]<-series[95-cut]*3
series[234-cut]<-series[234-cut]*2
plot(series,type='l')
title(main='Series plot')
length(series)
L<-60
s <- ssa(series, L = L)

plot(wcor(s))
plot(s, type = "series", groups = 1:8)
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

gseries<-series
gseries[95-cut]<-NA
gseries[234-cut]<-NA
s.temp <- ssa(gseries, L = L, force.decompose = FALSE)
#ig <- igapfill(s.temp, groups = list(1:rnk))
ig <- gapfill(s.temp, groups = list(1:rnk))
sg <- ssa(ig, L = L)
rec.g <- reconstruct(sg, groups = list(c(1:rnk))) 
trend.season.g <- rec.g$F1

plot(series,col='black',type='l')
lines(trend.season,type='l',col='red',lw=2)
lines(trend.season.g,type='l',col='magenta',lw=2)
lines( Pr2,type='l',col='orange',lw=2)
lines( Pr3,type='l',col='green',lw=2)
lines(Pr0,type='l',col='blue',lw=2)
legend('topleft', c("SSA", "SSA no gaps", "IRLS loess","IRLS lowess","IRLS orig"),
       col=c("red","magenta","blue","orange","green"), lty=1, cex=1.3, lw=c(2,2,2,2))
title(main='Reconstructed series')
