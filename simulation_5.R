mvg_0<-function(n=100,p){
  mean<-rep(0,p)
  sigma<-diag(p)
  mvrnorm(n, mean, sigma) #这里是用了正态性假设
}
mvg_1<-function(n=100,p){
  mean<-c(rep(0.1,floor(p/2)),rep(0,p-floor(p/2)))
  sigma<-diag(p)
  for (i in 1:p) {
    for (j in 1:p) {
      if(0<abs(i-j) && abs(i-j)<=3){
        sigma[i,j]<-0.1
      }
    }
  }
  mvrnorm(n, mean, sigma) #这里是用了正态性假设
}

size_cen<-function(n,p,a){
  k<-0
  for (i in 1:10000) {
    cen<-cen_method(n,p,mvg_0(n,p))
    if(abs(cen)>=qnorm(1-a/2)){
      k<-k+1
    }
  }
  k/10000
}

size_chi<-function(n,p,a){
  k<-0
  q<-length(n)
  f<-(1/2)*p*(p+3)
  for (i in 1:10000) {
    chi<-chi_method(n,p,mvg_0(n,p))
    if(chi>=qchisq(1-a,f)){
      k<-k+1
    }
  }
  k/10000
}

power_chi<-function(n,p,a){
  k<-0
  q<-length(n)
  f<-(1/2)*p*(p+3)
  for (i in 1:10000) {
    chi<-chi_method(n,p,mvg_1(n,p))
    if(chi>=qchisq(1-a,f)){
      k<-k+1
    }
  }
  k/10000
}

power_cen<-function(n,p,a){
  k<-0
  for (i in 1:10000) {
    cen<-cen_method(n,p,mvg_1(n,p))
    if(abs(cen)>=qnorm(1-a/2)){
      k<-k+1
    }
  }
  k/10000
}


sim<-function(n,p,a){
  simu<-matrix(0,4,4)
  for (i in 1:4) {
    simu[i,1]<-size_cen(n,p[i],a)
    simu[i,2]<-size_chi(n,p[i],a)
    simu[i,3]<-power_cen(n,p[i],a)
    simu[i,4]<-power_chi(n,p[i],a)
  }
  simu
}