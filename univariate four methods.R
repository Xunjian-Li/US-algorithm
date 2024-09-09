root.UHE2 = function(a3,a2,a1,a0,a.max=10,m,x0=0){
  xt.save = c(); yt.save = c()
  xt = x0; y1t = a3*xt^m+a2*xt^2+a1*xt+a0
  xt.save[1] = xt; yt.save[1] = y1t
  k = 1; esp = 1
  while(esp>1e-8){### Applying SLUF method
    if(y1t>0){ b2 = a3*m*(m-1)*a.max^(m-2) }else{ b2 = 0 }
    A2 = b2/2+a2
    A1 = a3*m*xt^(m-1) - b2*xt + a1
    A0 = a3*(1-m)*xt^m + b2/2*xt^2 + a0
    Delta.det = A1^2-4*A2*A0
    if(Delta.det>=0){
      xt.save[k+1] = -(A1+sqrt(Delta.det))/(2*A2)
      yt.save[k+1] = a3*xt.save[k+1]^m+a2*xt.save[k+1]^2+a1*xt.save[k+1] + a0
      esp = abs(yt.save[k+1])
      xt = xt.save[k+1]
      y1t = yt.save[k+1]
      k = k+1
    }else{
      stop(paste('There exists no root in (0,',a.max,")"))
    }
  }
  return(c(xt,k-1))
  # return(rbind(xt.save,yt.save))
}


root.UHE3 = function(a3,a2,a1,a0,a.max=10,m,x0=0){
  xt.save = c(); yt.save = c()
  xt = x0; y1t = a3*xt^m+a2*xt^2+a1*xt+a0
  xt.save[1] = xt; yt.save[1] = y1t
  k = 1; esp = 1
  while(esp>1e-8){### Applying TLF method
    A2 = a3*m*(m-1)/2*xt^(m-2)+a2
    A1 = a3*m*(2-m)*xt^(m-1) + a1
    A0 = a3*(m-1)*(m/2-1)*xt^m + a0
    Delta.det = A1^2-4*A2*A0
    if(Delta.det>=0){
      xt.save[k+1] = -(A1+sqrt(Delta.det))/(2*A2)
      yt.save[k+1] = a3*xt.save[k+1]^m+a2*xt.save[k+1]^2+a1*xt.save[k+1] + a0
      esp = abs(yt.save[k+1])
      xt = xt.save[k+1]
      y1t = yt.save[k+1]
      k = k+1
    }else{
      stop(paste('There exists no root in (0,',a.max,")"))
    }
  }
  return(c(xt,k-1))
  # return(rbind(xt.save,yt.save))
}


root.NR = function(a3,a2,a1,a0,m){
  xt.save.2 = c(); yt.save.2 = c(); xt = x0
  y1t = a3*xt^m+a2*xt^2+a1*xt+a0
  xt.save.2[1] = xt; yt.save.2[1] = y1t
  k = 1; esp = 1
  while(esp>1e-12){
    g.prime = m*a3*xt^(m-1)+2*a2*xt+a1
    xt.save.2[k+1] = xt - y1t/g.prime
    yt.save.2[k+1] = a3*xt.save.2[k+1]^m+a2*xt.save.2[k+1]^2+a1*xt.save.2[k+1]+a0
    esp = abs(yt.save.2[k+1])
    xt = xt.save.2[k+1]
    y1t = yt.save.2[k+1]
    k = k+1
    if(is.infinite(xt)){
      # return(rbind(xt.save.2,yt.save.2))
      return(c(xt,k-1))
      stop("Please choose a proper initial value")
    }
  }
  return(c(xt,k-1))
  # return(rbind(xt.save.2,yt.save.2))
}


root.BS = function(a3,a2,a1,a0,m){
  a.bound = c(0,2)
  xt.save.2 = c(); yt.save.2 = c(); xt = sum(a.bound)/2
  y1t = a3*xt^m+a2*xt^2+a1*xt+a0
  xt.save.2[1] = xt; yt.save.2[1] = y1t
  k = 1; esp = 1
  while(esp>1e-12){
    a.bound[2-sum(yt.save.2[k]>0)] = sum(a.bound)/2
    xt.save.2[k+1] = sum(a.bound)/2
    yt.save.2[k+1] = a3*xt.save.2[k+1]^m+a2*xt.save.2[k+1]^2+a1*xt.save.2[k+1] + a0
    esp = abs(yt.save.2[k+1])
    k = k+1
  }
  return(c(xt.save.2[k],k-1))
  # return(rbind(xt.save.2,yt.save.2))
}
# est.res5 = root.BS(a3,a2,a1,a0,m)
# est.res = est.res5