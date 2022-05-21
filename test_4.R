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
#上述样本生成器的生成样本格式为将每个总体的样本依次拼接，为行向量形式，即x2j依次接在x1j下方，这样一来
#第i个总体的样本就位于sum(n[1:i-1])+1到sum(n[1:i])行之间
chi_method<-function(n,p,x){
  k<-length(n)
  y<-matrix(0,k,p)
  N<-sum(n)
  a<-rep(0,k)
  A<-matrix(0,p,p)
  for (i in 1:k) {
    Ai<-(n[i]-1)*cov(x[,,i])
    a[i]<-det(Ai)
    A<-A+Ai
  }
  #以上为A的计算，a为一个行向量，其每个元素分别保存了Ai的行列式
  log_Lam<-sum(((n-1)/2)*log(a))-((N-k)/2)*log(det(A))+((N-k)*p/2)*log(N-k)-sum((p*(n-1)/2)*log(n-1))
  rho<-1-(sum((N-k)/(n-1))-1)*(2*p^2+3*p-1)/(6*(k-1)*(p+1)*(N-k))
  -2*rho*log_Lam
}
draw_chi<-function(n,p){
  N<-sum(n)
  k<-length(n)
  f<-(1/2)*p*(k-1)*(p+1)
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
  a<-rep(0,k)
  A<-matrix(0,p,p)
  for (i in 1:k) {
    Ai<-(n[i]-1)*cov(x[,,i])
    a[i]<-det(Ai)
    A<-A+Ai
  }
  log_Lam<-sum(((n-1)/2)*log(a))-((N-k)/2)*log(det(A))+((N-k)*p/2)*log(N-k)-sum((p*(n-1)/2)*log(n-1))
  u_n<-(1/4)*((N-k)*log(1-p/(N-k))*(2*N-2*p-2*k-1)-sum((n-1)*log(1-p/(n-1))*(2*n-2*p-3)))
  sigma2<-(1/2)*(log(1-p/(N-k))-sum(log(1-p/(n-1))*((n-1)/(N-k))^2))
  (log_Lam-u_n)/((N-k)*sqrt(sigma2))
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
