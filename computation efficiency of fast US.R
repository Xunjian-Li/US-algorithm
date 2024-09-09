est_theta_func_NR = function(theta0 = 1){
  esp = 1
  thetat = theta0
  k = 1
  while(esp>1e-8){
    C1 = 0
    C2 = 0
    for(i in c(1:n)){
      x.len = min(x[i],1000)
      C1 = C1 - sum(1/(c(0:(x.len-1))+thetat+1))
      C2 = C2 + sum(1/(c(0:(x.len-1))+thetat+1)^2)
    }
    gt = n/thetat + C1
    g.prime = -n/thetat^2 + C2
    theta.t1 = thetat - gt/g.prime
    esp = abs(thetat-theta.t1)
    thetat = theta.t1
    k= k+1
    
    if(is.na(esp)){thetat=0; break}
  }
  return(c(k-1,thetat))
}

est_theta_func_fast_US = function(theta0 = 1){
  C= n
  esp = 1
  thetat = theta0
  k = 1
  while(esp>1e-8){
    C1 = 0
    C2 = 0
    for(i in c(1:n)){
      x.len = min(x[i],1000)
      C1 = C1 - sum(1/(c(0:(x.len-1))+thetat+1))
      C2 = C2 + sum(1/(c(0:(x.len-1))+thetat+1)^2)
    }
    gt = n/thetat + C1
    g.prime = -n/thetat^2 + C2
    b_theta = -n/thetat^2 + C/(thetat+1)^2
    A1 = gt - n/thetat + C/(thetat+1)
    A2 = A1 + n - C
    theta.t1 = -(A2+sqrt(A2^2-4*A1*n))/(2*A1)
    s.ratio = min(b_theta/g.prime*(g.prime<0) + (g.prime>0),2)
    theta.t1 = thetat + s.ratio*(theta.t1 - thetat)
    if(theta.t1<0){theta.t1 = 1e-8}
    esp = abs(thetat-theta.t1)
    thetat = theta.t1
    k= k+1
  }
  return(c(k-1,thetat))
}
  
est_theta_func_US = function(theta0 = 1){
  C= n
  esp = 1
  thetat = theta0
  k = 1
  while(esp>1e-8){
    C1 = 0
    for(i in c(1:n)){
      x.len = min(x[i],1000)
      C1 = C1 - sum(1/(c(0:(x.len-1))+thetat+1))
    }
    gt = n/thetat + C1
    b_theta = -n/thetat^2 + C/(thetat+1)^2
    A1 = gt - n/thetat + C/(thetat+1)
    A2 = A1 + n - C
    theta.t1 = -(A2+sqrt(A2^2-4*A1*n))/(2*A1)
    esp = abs(thetat-theta.t1)
    thetat = theta.t1
    k= k+1
  }
  return(c(k-1,thetat))
}


est_theta_func_fast_FP = function(theta0 = 1){
  esp = 1
  thetat = theta0
  k = 1
  while(esp>1e-8){
    C1 = 0
    C2 = 0
    for(i in c(1:n)){
      x.len = min(x[i],1000)
      C1 = C1 - sum(1/(c(0:(x.len-1))+thetat+1))
      C2 = C2 + sum(1/(c(0:(x.len-1))+thetat+1)^2)
    }
    gt = n/thetat + C1
    g.prime = -n/thetat^2 + C2
    b_theta = -n/thetat^2
    A1 = gt - n/thetat
    A2 = A1 + n
    theta.t1 = -(A2+sqrt(A2^2-4*A1*n))/(2*A1)
    s.ratio = min(b_theta/g.prime*(g.prime<0) + (g.prime>0),2)
    theta.t1 = thetat + s.ratio*(theta.t1 - thetat)
    if(theta.t1<0){theta.t1 = 1e-8}
    esp = abs(thetat-theta.t1)
    thetat = theta.t1
    k= k+1
  }
  return(c(k-1,thetat))
}

est_theta_func_FP = function(theta0 = 1){
  thetat = theta0
  k = 1
  esp = 1
  while(esp>1e-8){
    C1 = 0
    for(i in c(1:n)){
      x.len = min(x[i],1000)
      C1 = C1 - sum(1/(c(0:(x.len-1))+thetat+1))
    }
    gt = n/thetat + C1
    b_theta = -n/thetat^2
    A1 = gt - n/thetat
    A2 = A1 + n
    theta.t1 = -(A2+sqrt(A2^2-4*A1*n))/(2*A1)
    esp = abs(thetat-theta.t1)
    thetat = theta.t1
    k= k+1
  }
  return(c(k-1,thetat))
}


est_theta_func_SLF = function(theta0 = 1){
  esp = 1
  theta.all = c()
  thetat = theta0
  theta.all[1] = thetat
  n = length(x)
  k = 1
  while(esp>1e-8){
    C1 = 0
    G = 0
    a3.tra = 0
    a4.tra = 0
    for(i in c(1:n)){
      x.len = min(x[i],1000)
      r.i = 1/(c(0:(x.len-1))+thetat+1)
      C1 = C1 - sum(r.i)
      G = G + sum(r.i^2)
      a3.tra = a3.tra + sum(r.i^3)
      a4.tra = a4.tra + sum(c(1:x.len)*r.i^3)
    }
    a3 = G^3/a3.tra^2
    a4 = a4.tra/a3.tra
    gt = n/thetat + C1
    A1 = gt - n/thetat + a3/(thetat+a4)
    A2 = A1*a4 + n-a3
    A3 = n*a4
    theta.t1 = -(A2+sqrt(A2^2-4*A1*A3))/(2*A1)
    esp = abs(thetat-theta.t1)
    thetat = theta.t1
    k= k+1
    theta.all[k] = theta.t1
  }
  return(c(k,theta.all))
}

library(VGAM)
library(lubridate)
rest.all = matrix(0,24,4)
k.num = 1
times.rep = 100

for(kk in c(1:1)){
  for(tt in c(3:3)){
    theta.ori = c(0.5,1,5,10)[kk]
    n = c(100,200,400)[tt]
    for(zz in c(3:3)){
      time_start = Sys.time()
      times.all = 1
      k.all = c()
      est.all = c()
      while(times.all <= times.rep){
        x = ryules(n, theta.ori)
        theta0 = runif(1,1,5)
        if(zz == 1){
          est1 = est_theta_func_US(theta0)
        }else if(zz == 2){
          est1 = est_theta_func_fast_US(theta0)
        }else if(zz == 3){
          est1 = est_theta_func_NR(theta0)
        }else if(zz == 4){
          est1 = est_theta_func_FP(theta0)
        }else if(zz == 5){
          est1 = est_theta_func_fast_FP(theta0)
        }else{
          est1 = est_theta_func_SLF(theta0)
        }
        k.all[times.all] = est1[1]
        est.all[times.all] = est1[2]
        times.all = times.all +1
      }
      time_end = Sys.time()
      rest.all[k.num,] = c(k.num, as.duration(time_end - time_start),mean(k.all[est.all>0]),sum(est.all<0)/length(est.all))
      k.num = k.num+1
    }
  }
}

rest.all


sum(rest.all[,2])














time_start = Sys.time()
times.all = 1
k.all = c()
est.all = c()
theta.ori = 0.5
theta0 = runif(1,1,5)
while(times.all <= 100){
  x = ryules(n, theta.ori)
  theta0 = runif(1,1,5)
  # est1 = est_theta_func_fast_FP(theta0)
  # est1 = est_theta_func_FP(theta0)
  est1 = est_theta_func_fast_US(theta0)
  # est1 = est_theta_func_US(theta0)
  # est1 = est_theta_func_NR(theta0)
  k.all[times.all] = est1[1]
  est.all[times.all] = est1[2]
  times.all = times.all +1
}
time_end = Sys.time()
time_end - time_start
mean(k.all)
sum(est.all<=0)/length(est.all)