#### case 1: Root of a high-order polynomial equation

# calculating unique higher-order equation with SLUB method
root.UHE = function(a3, a2, a1, a0, a.max=10, m, x0 = 0)
{ 
  xt.save = c(); yt.save = c(); 
  xt = x0; y1t = a3*xt^m+a2*xt^2+a1*xt+a0
  xt.save[1] = xt; yt.save[1] = y1t
  k = 1; esp = 1
  while(esp>1e-8){
    if(a3<0){                ### When a3<0, we apply SLUB method
      if(y1t>0){ b2 = a3*m*(m-1)*a.max^(m-2) }else{b2 = 0}
      A2 = b2/2+a2
      A1 = a3*m*xt^(m-1)-b2*xt + a1
      A0 = a3*(1-m)*xt^m+b2/2*xt^2 + a0 
    }else{                   ### When a3>0, we apply TLF method
      A2 = a3*m*(m-1)/2*xt^(m-2)+a2
      A1 = a3*m*(2-m)*xt^(m-1)+a1
      A0 = a3*(m-1)*(m/2-1)*xt^m+a0
    }
    Delta.det = A1^2-4*A2*A0
    if(Delta.det>=0){
      xt.save[k+1] = -(A1+sqrt(Delta.det))/(2*A2)
      yt.save[k+1] = a3*xt.save[k+1]^m+a2*xt.save[k+1]^2+
        a1*xt.save[k+1]+a0
      esp = abs(yt.save[k+1])
      xt = xt.save[k+1]
      y1t = yt.save[k+1]
      k = k+1
    }else{
      stop(paste('There exists no root in (0,',a.max,")"))
    }
  }
  return(c(xt,k-1))
}

# calculating unique higher-order equation with TLB method
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

# calculating unique higher-order equation with NR method
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

# calculating unique higher-order equation with Bisection method
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

## comparisons of NR Bisec and US
a3 = -1; a2 = 1; a1 = -1; a0 = 1; m = 3; n.rep = 100000
start.time = Sys.time()
num.inva0 = 0
num.iter0 = 0
for(i in c(1:n.rep)){
  x0 = runif(1,0,2)
  est.res3 = root.UHE3(a3,a2,a1,a0,a.max=2,m,x0=0)
  if(is.infinite(est.res3[1])){
    num.inva0 = num.inva0 + 1
  }else{
    num.iter0 = num.iter0 + est.res3[2]
  }
}
end.time = Sys.time()
time.us0 = end.time - start.time

a3 = 1; a2 = -3; a1 = -1; a0 = 1; m = 3; n.rep = 100000
start.time = Sys.time()
num.inva1 = 0
num.iter1 = 0
for(i in c(1:n.rep)){
  x0 = runif(1,0,2)
  est.res3 = root.UHE3(a3,a2,a1,a0,a.max=2,m,x0=0)
  if(is.infinite(est.res3[1])){
    num.inva1 = num.inva1 + 1
  }else{
    num.iter1 = num.iter1 + est.res3[2]
  }
}
end.time = Sys.time()
time.us = end.time - start.time

start.time = Sys.time()
num.inva2 = 0
num.iter2 = 0
for(i in c(1:n.rep)){
  x0 = runif(1,0,2)
  est.res3 = root.NR(a3,a2,a1,a0,m)
  if(is.infinite(est.res3[1])){ 
    num.inva2 = num.inva2 + 1
  }else{
    num.iter2 = num.iter2 + est.res3[2]
  }
}
end.time = Sys.time()
time.nr = end.time - start.time

start.time = Sys.time()
num.inva3 = 0
num.iter3 = 0
for(i in c(1:n.rep)){
  x0 = runif(1,0,2)
  est.res3 = root.BS(a3,a2,a1,a0,m)
  if(is.infinite(est.res3[1])){ 
    num.inva3 = num.inva3 + 1
  }else{
    num.iter3 = num.iter3 + est.res3[2]
  }
}
end.time = Sys.time()
time.bs = end.time - start.time

time.tra = rbind(c(time.us0, num.inva0,num.iter0/n.rep/(1-num.inva0)),
                 # c(time.us, num.inva1,num.iter1/n.rep/(1-num.inva1)), 
                 c(time.nr, num.inva2,num.iter2/(n.rep-num.inva2)), 
                 c(time.bs, num.inva3,num.iter3/n.rep/(1-num.inva3)))
time.tra
