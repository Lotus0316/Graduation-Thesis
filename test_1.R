library(MASS)
library(Matrix)
mvg<-function(n=100,p){
  mean<-rep(0,p)
  sigma<-diag(p)
  mvrnorm(n, mean, sigma) #这里是用了正态性假设
}

#传统卡方近似
chi_method<-function(n,p,x){
  S<-((n-1)/n)*cov(x)
  rho<-1-(2*p^2+p+2)/((6*(n-1)*p))
  -(n-1)*rho*(log(det(S))+(-p)*(log(tr(S))-log(p)))
}
draw_chi<-function(n,p){
  f<-(1/2)*(p-1)*(p+2)
  x<-replicate(10000,
               chi_method(n,p,mvg(n,p)))
  hist(x,
       breaks=40,
       freq=F,
       xlab=NULL,
       main = paste("n=",n,",","p=",p))
  curve(dchisq(x, f), 
        add=TRUE, 
        col="darkblue", 
        lwd=2) 
}

#中心极限方法
cen_method<-function(n,p,x){
  S<-((n-1)/n)*cov(x)
  log_V<-log(det(S))+(-p)*(log(tr(S))-log(p))
  u_n<--p-(n-p-3/2)*log(1-p/(n-1))
  sigma2<--2*(p/(n-1)+log(1-p/(n-1)))
  (log_V-u_n)/sqrt(sigma2)
}
draw_cen<-function(n,p){
  x<-replicate(10000,
               cen_method(n,p,mvg(n,p)))
  hist(x,
       breaks=40,
       freq=F,
       xlab=NULL,
       main = paste("n=",n,",","p=",p))
  curve(dnorm(x),
        add=TRUE, 
        col="darkblue", 
        lwd=2) 
}