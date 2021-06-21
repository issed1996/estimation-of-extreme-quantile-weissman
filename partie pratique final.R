
set.seed(789)
library(evmix)
loc=0
scale=1
shape=2
#n=100000
n=1000
data=rfrechet(n,shape,loc,scale)


library(ReIns)
hist(data,100)

#estimator of EVI:gamma
hill_estimators=hill(data) #Bn=110
H <- Hill(data,plot=TRUE)
abline(h=1/shape,col="blue")
legend("topleft",
       c("hill estimators","exact gamma"),
       fill=c("black","blue")
)


plot(H$k,H$gamma,type="l")
abline(h=alpha,col="blue")

x_max=400
h=20
l=n-h-1
v=c()
for (k0 in c(1:x_max)) {
  k1=k0+h
  m=var(H$gamma[c(k0:k1)])
  v=append(v,m) 
  #print(k0)
}

library(NNTbiomarker)
min(v)
kn=argmin(v)
kn=178

gamma_estimator=H$gamma[kn]
H$gamma[kn]
#estimator of quantile of order 1-pn while pn tends to 0

p <- 10^(-6)

res=Quant(data, gamma=H$gamma, p=p, plot=FALSE)
#the weissman estimator of quantile
weissman_estimator=res$Q[kn]
weissman_estimator
#classic estimator on quantile
classic_estimator=sort(data,decreasing = TRUE)[trunc(p*n)+1]
classic_estimator
#exact quantile
exact_quantile=log(1/(1-p))^(-1/shape)
exact_quantile

x <- 10^-seq(1.1,3.2,0.1)

y_weissman=c()
y_classic=c()
y_exact=c()
for (p  in x) {
  res=Quant(data, gamma=H$gamma, p=p, plot=FALSE)
  weissman_estimator=res$Q[kn]
  y_weissman=append(y_weissman,weissman_estimator)
  
  i=trunc(p*n)
  if (i==0) {
    i=1
  }
  classic_estimator=sort(data,decreasing = TRUE)[i]
  y_classic=append(y_classic,classic_estimator)
  
  exact_quantile=log(1/(1-p))^(-1/shape)
  y_exact=append(y_exact,exact_quantile)
  
}

y_weissman

y_classic

plot(log10(x),y_weissman,type="l",col="red",main="comparaison des estimateurs de quantiles ",ylab = 'quantile',xlab = 'log10(pn)')
lines(log10(x),y_exact, col="blue")
lines(log10(x),y_classic, col="green")
legend("topleft",
       c("y_weissman","y_exact","y_classic"),
       fill=c("red","blue","green")
)


#plot(log(x),y_classic,type="l",col="blue")
#plot(log(x),y_exact,type="l",col="green")




library(ReIns)
?soa

# Look at last 500 observations of SOA data
SOAdata <- sort(soa$size)[length(soa$size)-(0:499)]

hist(SOAdata,100)
# Hill estimator
H <- Hill(SOAdata,plot=T)
# Bias-reduced estimator (QV)
H_QV <- Hill.2oQV(SOAdata)

# Exceedance probability
p <- 10^(-5)
# Weissman estimator
Quant(SOAdata, gamma=H$gamma, p=p, plot=TRUE)



x_max=400
h=20
l=n-h-1
v=c()
for (k0 in c(1:x_max)) {
  k1=k0+h
  m=var(H$gamma[c(k0:k1)])
  v=append(v,m) 
  #print(k0)
}

library(NNTbiomarker)
min(v)
kn=argmin(v)
kn

gamma_estimator=H$gamma[kn]
H$gamma[kn]

x <- 10^-seq(1.1,3,0.1)

y_weissman=c()
y_classic=c()
y_exact=c()
for (p  in x) {
  res=Quant(SOAdata, gamma=H$gamma, p=p, plot=FALSE)
  weissman_estimator=res$Q[kn]
  y_weissman=append(y_weissman,weissman_estimator)
  
  classic_estimator=sort(SOAdata,decreasing = TRUE)[trunc(p*n)+1]
  y_classic=append(y_classic,classic_estimator)
  
  #exact_quantile=log(1/(1-p))^(-1/shape)
  #y_exact=append(y_exact,exact_quantile)
  
}
plot(log10(x),y_weissman,type="l",col="red",main="comparaison des estimateurs de quantiles ",ylab = 'quantile',xlab = 'log10(pn)')
#lines(log10(x),y_exact, col="blue")
lines(log10(x),y_classic, col="green")
legend("topleft",
       c("y_weissman","y_classic"),
       fill=c("red","green")
)

plot(log(x),y_classic,type="l",col="blue")


n=300
X=runif(n)
X
gamma_x=(1/2)*(1/10+sin(3.14*X))*(11/10-(1/2)*exp(-64*(X-1/2)^2))
gamma_x








set.seed(789)

loc=1
shape=10
data=rpareto(n,shape,loc)
hist(data,100)

H <- Hill(data,plot=TRUE)
abline(h=1/shape,col="blue")
legend("topleft",
       c("hill estimators","exact gamma"),
       fill=c("black","blue")
)


plot(H$k,H$gamma,type="l")
abline(h=alpha,col="blue")

x_max=400
h=20
l=n-h-1
v=c()
for (k0 in c(1:x_max)) {
  k1=k0+h
  m=var(H$gamma[c(k0:k1)])
  v=append(v,m) 
  #print(k0)
}

library(NNTbiomarker)
min(v)
kn=argmin(v)
kn
gamma_estimator=H$gamma[kn]
#estimateur de gamma=1/shape
H$gamma[kn]


x <- 10^-seq(1.1,4,0.1)

y_weissman=c()
y_classic=c()
y_exact=c()
for (p  in x) {
  res=Quant(data, gamma=H$gamma, p=p, plot=FALSE)
  weissman_estimator=res$Q[kn]
  y_weissman=append(y_weissman,weissman_estimator)
  
  i=trunc(p*n)
  if (i==0) {
    i=1
  }
  classic_estimator=sort(data,decreasing = TRUE)[i]
  y_classic=append(y_classic,classic_estimator)
  
  exact_quantile=(p^(-1/shape))*loc
  y_exact=append(y_exact,exact_quantile)
  
}

y_weissman

y_classic

plot(log10(x),y_weissman,type="l",col="red",main="comparaison des estimateurs de quantiles ",ylab = 'quantile',xlab = 'log10(pn)')
lines(log10(x),y_exact, col="blue")
lines(log10(x),y_classic, col="green")
legend("topleft",
       c("y_weissman","y_exact","y_classic"),
       fill=c("red","blue","green")
)
