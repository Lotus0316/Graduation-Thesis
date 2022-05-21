library(MASS)
library(Matrix)
mvg<-function(n,p){
  mean<-rep(0,sum(p))
  sigma<-diag(sum(p))
  y<-array(0,dim=c(max(n),p,length(n)))
  for (i in 1:length(n)) {
    y[,,i]<-mvrnorm(n[i], mean, sigma)
  }
  y
}
chi_method<-function(n,p,x){
  k<-length(n)
  y<-matrix(0,k,p)
  N<-sum(n)
  for (i in 1:k) {
    y[i,]<-colMeans(x[,,i])
  }
  #这样就构造了一个y矩阵，其每一行分别对应着第i个总体的样本均值
  A<-matrix(0,p,p)
  for (i in 1:k) {
    A<-A+n[i]*(y[i,]-colMeans(y))%o%(y[i,]-colMeans(y))
  }
  #以上为A矩阵的计算
  b<-rep(0,k)
  B<-matrix(0,p,p)
  for (i in 1:k) {
    Bi<-diag(0,p)
    for (j in 1:n[i]) {
      Bi<-Bi+(x[j,,i]-y[i,])%o%(x[j,,i]-y[i,])
    }
    b[i]<-det(Bi)
    B<-B+Bi
  }
  #以上为B的计算，b为一个行向量，其每个元素分别保存了Bi的行列式
  #接下来计算log Lambda_n
  log_Lam<-sum((n/2)*log(b))-(N/2)*log(det(A+B))+(N*p/2)*log(N)-sum((p*n/2)*log(n))
  rho<-1-(sum(N/n)-1)*(2*p^2+9*p+11)/(6*(k-1)*(p+3)*N)
  -2*rho*log_Lam
}
draw_chi<-function(n,p){
  N<-sum(n)
  k<-length(n)
  f<-(1/2)*p*(k-1)*(p+3)
  x<-replicate(10000,
               chi_method(n,p,mvg(n,p)))
  hist(x,
       breaks=40,
       freq=F,
       xlab=NULL,
       main = paste("n=",paste(n,collapse=","),",","p=",p))
  curve(dchisq(x, f), 
        add=TRUE, 
        col="darkblue", 
        lwd=2) 
}

cen_method<-function(n,p,x){
  k<-length(n)
  y<-matrix(0,k,p)
  N<-sum(n)
  for (i in 1:k) {
    y[i,]<-colMeans(x[,,i])
  }
  #这样就构造了一个y矩阵，其每一行分别对应着第i个总体的样本均值
  A<-matrix(0,p,p)
  for (i in 1:k) {
    A<-A+n[i]*(y[i,]-colMeans(y))%o%(y[i,]-colMeans(y))
  }
  #以上为A矩阵的计算
  b<-rep(0,k)
  B<-matrix(0,p,p)
  for (i in 1:k) {
    Bi<-diag(0,p)
    for (j in 1:n[i]) {
      Bi<-Bi+(x[j,,i]-y[i,])%o%(x[j,,i]-y[i,])
    }
    b[i]<-det(Bi)
    B<-B+Bi
  }
  #以上为B的计算，b为一个行向量，其每个元素分别保存了Bi的行列式
  #接下来计算log Lambda_n
  log_Lam<-sum((n/2)*log(b))-(N/2)*log(det(A+B))+(N*p/2)*log(N)-sum((p*n/2)*log(n))
  #接下来计算渐进均值和方差
  u_n<-(1/4)*(-2*k*p-sum(p/n)+N*(-log(1-p/N))*(2*p-2*N+3)-sum(n*(-log(1-p/(n-1)))*(2*p-2*n+3)))
  sigma2<-(1/2)*(sum((-log(1-p/(n-1)))*(n^2)/(N^2))-(-log(1-p/N)))
  (log_Lam-u_n)/(N*sqrt(sigma2))
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