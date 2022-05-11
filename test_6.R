library(ggplot2)
library(lava)
library(MASS)
library(dplyr)
library(patchwork)#加载包
mvg<-function(n=100,p){
  mean<-rep(0,p)
  sigma<-diag(p)
  mvrnorm(n, mean, sigma) #这里是用了正态性假设
}

#传统卡方近似
chi_method<-function(n,p,x){
  -(n-1-(2*p+5)/6)*log(det(cor(x)))
}
draw_chi<-function(n,p){
  f<-(1/2)*(p-1)*p
  x<-replicate(10000,
               chi_method(n,p,mvg(n,p)))
  hist(x,
       breaks=40,
       freq=F,
       xlab=NULL,
       ylim=c(0,0.01),
       main = paste("n=",n,",","p=",p))
  curve(dchisq(x, f), 
        add=TRUE, 
        col="darkblue", 
        lwd=2) 
}

#中心极限方法
cen_method<-function(n,p,x){
  u_n<-(p-n+3/2)*log(1-p/(n-1))-(n-2)*p/(n-1)
  sigma2<--2*(p/(n-1)+log(1-p/(n-1)))
  (log(det(cor(x)))-u_n)/sqrt(sigma2)
}
draw_cen<-function(n,p){
  x<-replicate(10000,
               cen_method(n,p,mvg(n,p)))
  hist(x,
       breaks=40,
       freq=F,
       xlab=NULL,
       ylim=c(0,0.01),
       main = paste("n=",n,",","p=",p))
  curve(dnorm(x),
        add=TRUE, 
        col="darkblue", 
        lwd=2) 
}