library(MASS)
library(Matrix)
mvg<-function(n,p){
  mean<-rep(0,sum(p))
  sigma<-diag(sum(p))
  mvrnorm(n, mean, sigma) #这里是用了正态性假设

}
chi_method<-function(n,p,x){
  A<-(n-1)*cov(x)
  log_lam<-n*p/2-n*p*log(n)/2+(n/2)*log(det(A))-tr(A)/2-n*colMeans(x)%*%colMeans(x)/2
  rho<-1-(2*p^2+9*p+11)/(6*n*(p+3))
  -2*rho*log_lam
}
draw_chi<-function(n,p){
  f<-(1/2)*p*(p+3)
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
cen_method<-function(n,p,x){
  A<-(n-1)*cov(x)
  log_lam<-n*p/2-n*p*log(n)/2+(n/2)*log(det(A))-tr(A)/2-n*colMeans(x)%*%colMeans(x)/2
  u_n<--(1/4)*(n*(2*n-2*p-3)*log(1-p/(n-1))+2*(n+1)*p)
  sigma2<--(1/2)*(p/(n-1)+log(1-p/(n-1)))
  (log_lam-u_n)/(n*sqrt(sigma2))
}
draw_cen<-function(n,p){
  x<-replicate(10000,
               cen_method(n,p,mvg(n,p)))
  hist(x,
       breaks=40,
       freq=F,
       xlab=NULL,
       main = paste("n=",paste(n,collapse=","),",","p=",p))
  curve(dnorm(x),
        add=TRUE, 
        col="darkblue", 
        lwd=2) 
}