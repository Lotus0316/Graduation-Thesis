library(MASS)
library(Matrix)
mvg<-function(n=100,p){
  mean<-rep(0,sum(p))
  sigma<-diag(sum(p))
  mvrnorm(n, mean, sigma) #这里是用了正态性假设
}

#传统卡方检验
chi_method<-function(n,p=rep(1,3),x){
  C<-cov(x)*(n-1)
  y<-1
  for (i in 2:length(p)) {
    if(length((sum(p[1:i-1])+1):sum(p[1:i]))==1){
      y<-y*C[c((sum(p[1:i-1])+1):sum(p[1:i]))
             ,c((sum(p[1:i-1])+1):sum(p[1:i]))]
    }
    else{
      y<-y*det(C[c((sum(p[1:i-1])+1):sum(p[1:i]))
                 ,c((sum(p[1:i-1])+1):sum(p[1:i]))])
    }
  
  }
  denomi<-det(C[c(1:p[1]),c(1:p[1])])*y
  Lambda_n<-(det(C)/denomi)
  rho<-1-(2*(sum(p)^3-sum(p^3))+9*(p^2-sum(p)^2))/(6*n*(sum(p)^2-sum(p^2)))
  -2*rho*log(Lambda_n)*(n/2)
  #原本Lambda_n的定义是W_n的n/2次方，但为了避免乘方过大导致Lambda_n溢出直接
  #log处补上系数n/2
}

#这里的p应该是个向量，其元素为p1-pk
draw_chi<-function(n,p){
  f<-(1/2)*(sum(p)^2-sum(p^2))
  x<-replicate(10000,
               chi_method(n,p,mvg(n,p)))
  hist(x,
       breaks=40,
       freq=F,
       xlab=NULL,
       ylim=c(0,0.01),
       main = paste("n=",n,",","p=",paste(p,collapse=",")))
  curve(dchisq(x, f), 
        add=TRUE, 
        col="darkblue", 
        lwd=2) 
}

cen_method<-function(n,p=rep(1,3),x){
  C<-cov(x)*(n-1)
  y<-1
  for (i in 2:length(p)) {
    if(length((sum(p[1:i-1])+1):sum(p[1:i]))==1){
      y<-y*C[c((sum(p[1:i-1])+1):sum(p[1:i]))
             ,c((sum(p[1:i-1])+1):sum(p[1:i]))]
    }
    else{
      y<-y*det(C[c((sum(p[1:i-1])+1):sum(p[1:i]))
                 ,c((sum(p[1:i-1])+1):sum(p[1:i]))])
    }
    
  }
  denomi<-det(C[c(1:p[1]),c(1:p[1])])*y
  Lambda_n<-(det(C)/denomi)
  u_n<--(-log(1-sum(p)/(n-1)))*(sum(p)-(n-1)+1/2)+sum((-log(1-p/(n-1)))*(p-(n-1)+1/2))
  sigma2<-2*(-log(1-sum(p)/(n-1)))-2*sum((-log(1-p/(n-1))))
  (log(Lambda_n)-u_n)/sqrt(sigma2)
  #由于文章中定义的n=N-1，而代码中直接用n记作样本量，所以公式中的n在代码中均为n-1，文中的n+1均为n
}

draw_cen<-function(n,p){
  x<-replicate(10000,
               cen_method(n,p,mvg(n,p)))
  hist(x,
       breaks=40,
       freq=F,
       xlab=NULL,
       main = paste("n=",n,",","p=",paste(p,collapse=",")))
  curve(dnorm(x),
        add=TRUE, 
        col="darkblue", 
        lwd=2) 
}