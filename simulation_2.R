mvg<-function(n=100,p,sigma){
  mean<-rep(0,sum(p))
  mvrnorm(n, mean, sigma) #这里是用了正态性假设
}

size_cen<-function(n,p,a,sigma){
  k<-0
  for (i in 1:10000) {
    cen<-cen_method(n,c(2*p/5,2*p/5,p/5),mvg(n,p,sigma))
    if(abs(cen)>=qnorm(1-a/2)){
      k<-k+1
    }
  }
  k/10000
}

size_chi<-function(n,p,a,sigma){
  k<-0
  f<-(1/2)*(sum(c(2*p/5,2*p/5,p/5))^2-sum(c(2*p/5,2*p/5,p/5)^2))#要看情况改动
  for (i in 1:10000) {
    chi<-chi_method(n,c(2*p/5,2*p/5,p/5),mvg(n,p,sigma))
    if(chi>=qchisq(1-a,f)){
      k<-k+1
    }
  }
  k/10000
}

power_chi<-function(n,p,a,sigma){
  k<-0
  f<-(1/2)*(sum(c(2*p/5,2*p/5,p/5))^2-sum(c(2*p/5,2*p/5,p/5)^2))#要看情况改动
  for (i in 1:10000) {
    chi<-chi_method(n,c(2*p/5,2*p/5,p/5),mvg(n,p,sigma))
    if(chi>=qchisq(1-a,f)){
      k<-k+1
    }
  }
  k/10000
}

power_cen<-function(n,p,a,sigma){
  k<-0
  for (i in 1:10000) {
    cen<-cen_method(n,c(2*p/5,2*p/5,p/5),mvg(n,p,sigma))
    if(abs(cen)>=qnorm(1-a/2)){
      k<-k+1
    }
  }
  k/10000
}


sim<-function(n,p,a){
  simu<-matrix(0,4,4)
  for (i in 1:4) {
    Sigma_1<-diag(rep(1,p[i]))
    Sigma_2<-0.15*matrix(1,p[i],p[i])+0.85*diag(p[i])
    simu[i,1]<-size_cen(n,p[i],a,Sigma_1)
    simu[i,2]<-size_chi(n,p[i],a,Sigma_1)
    simu[i,3]<-power_cen(n,p[i],a,Sigma_2)
    simu[i,4]<-power_chi(n,p[i],a,Sigma_2)
  }
  simu
}