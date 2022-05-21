mvg_0<-function(n,p){
  mean<-rep(0,sum(p))
  sigma<-diag(sum(p))
  y<-array(0,dim=c(max(n),p,length(n)))
  for (i in 1:length(n)) {
    y[,,i]<-mvrnorm(n[i], mean, sigma)
  }
  y
}

mvg_1<-function(n,p){
  mean<-array(0,dim=c(p,1,length(n)))
  sigma<-array(0,dim=c(p,p,length(n)))
  mean[,,1]<-rep(0,p)
  mean[,,2]<-rep(0.1,p)
  mean[,,3]<-rep(0,p)
  sigma[,,1]<-0.5*matrix(1,p,p)+0.5*diag(p)
  sigma[,,2]<-0.6*matrix(1,p,p)+0.4*diag(p)
  sigma[,,3]<-0.5*matrix(1,p,p)+0.31*diag(p)
  y<-array(0,dim=c(max(n),p,length(n)))
  for (i in 1:length(n)) {
    y[,,i]<-mvrnorm(n[i], mean[,,i], sigma[,,i])
  }
  y
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
  f<-(1/2)*p*(q-1)*(p+3)
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
  f<-(1/2)*p*(q-1)*(p+3)
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