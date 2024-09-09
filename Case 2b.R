#####case 2: 0.5-th quantile of normal distribution
############comparing efficiency of four methods###############

mu = 1; sig = 1; p = 0.05; x0 = 0
x0 = -1.5
###################SLUB algorithm#################
estimation_fun2.2 = function(mu,sig,p,x0=0){
  xt.save.2 = c(); yt.save.2 = c(); xt = x0
  y1t = p - pnorm(xt,mu,sig)
  xt.save.2[1] = xt; yt.save.2[1] = y1t
  k = 1; esp = 1
  while(esp>1e-12){
    g.prime = -dnorm(xt,mu,sig)
    if(y1t>0){
      b2 = -1/(sig^3*sqrt(2*pi)*exp(1/2))
    }else{
      b2 = 1/(sig^3*sqrt(2*pi)*exp(1/2))
    }
    A1 = b2/2; A2 = -dnorm(xt,mu,sig); A3 = p-pnorm(xt,mu,sig)
    # A1*(xt-xt)^2+A2*(xt-xt)+A3
    xt.save.2[k+1] = xt -(A2+sqrt(A2^2-4*A1*A3))/(2*A1)
    yt.save.2[k+1] = p-pnorm(xt.save.2[k+1],mu,sig)
    esp = abs(yt.save.2[k+1])
    xt = xt.save.2[k+1]
    y1t = yt.save.2[k+1]
    k = k+1
  }
  return(c(xt,k-1))
}

###################TLB algorithm##################
estimation_fun2.3 = function(mu,sig,p,x0=0){
  xt.save.2 = c(); yt.save.2 = c(); xt = x0
  y1t = p - pnorm(xt,mu,sig)
  xt.save.2[1] = xt; yt.save.2[1] = y1t
  k = 1; esp = 1
  while(esp>1e-12){
    g.prime = -dnorm(xt,mu,sig)
    g.primeprime = (xt-mu)/sig^2*dnorm(xt,mu,sig)
    a3 = -1/(3*sig^3*sqrt(2*pi)*exp(3/2))
    a2 = g.primeprime/2
    a1 = g.prime
    a0 = y1t
    # res11 = root.UHE.revise(a3,a2,a1,a0,a.max=10,m=3,x0=0)
    # xt.save.2[k+1] = xt + res11[1,dim(res11)[2]]
    res11 = root_ploynomial_3(a3,a2,a1,a0)
    if(length(res11)==1){
      xt.save.2[k+1] = xt + res11[1]
    }else{
      if(y1t>0){
        xt.save.2[k+1] = xt + min(res11[res11>=0])
      }else{
        xt.save.2[k+1] = xt + max(res11[res11<=0])
      }
    }
    yt.save.2[k+1] = p-pnorm(xt.save.2[k+1],mu,sig)
    esp = abs(yt.save.2[k+1])
    xt = xt.save.2[k+1]
    y1t = yt.save.2[k+1]
    k = k+1
  }
  return(c(xt,k-1))
}

estimation_fun2.NR = function(mu,sig,p,x0=0){
  # p = 0.00000001
  xt.save.2 = c(); yt.save.2 = c(); xt = x0
  y1t = p - pnorm(xt,mu,sig)
  xt.save.2[1] = xt; yt.save.2[1] = y1t
  k = 1; esp = 1
  while(esp>1e-12){
    g.prime = -dnorm(xt,mu,sig)
    xt.save.2[k+1] = xt - y1t/g.prime
    yt.save.2[k+1] = p-pnorm(xt.save.2[k+1],mu,sig)
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
}
estimation_fun2.BS = function(mu,sig,p,x0=0){
  # p = 0.00000001
  a.bound = c(-20,20)
  xt.save.2 = c(); yt.save.2 = c(); xt = sum(a.bound)/2
  y1t = p - pnorm(xt,mu,sig)
  xt.save.2[1] = xt; yt.save.2[1] = y1t
  k = 1; esp = 1
  while(esp>1e-12){
    a.bound[2-sum(yt.save.2[k]>0)] = sum(a.bound)/2
    xt.save.2[k+1] = sum(a.bound)/2
    yt.save.2[k+1] = p-pnorm(xt.save.2[k+1],mu,sig)
    esp = abs(yt.save.2[k+1])
    k = k+1
  }
  return(c(xt,k-1))
}

mu = -2; sig = 1; p = 0.01;  ## x.star = -4.326348
mu = 2;  sig = 1; p = 0.01;  ## x.star = -0.3263479
mu = -2; sig = 1; p = 0.9;   ## x.star = -0.7184484
mu = 2;  sig = 1; p = 0.9;   ## x.star = 3.281552

n.rep = 100000

start.time = Sys.time()
num.inva0 = 0
num.iter0 = 0
for(i in c(1:n.rep)){
  x0 = runif(1,-4,4)
  est.res3 = estimation_fun2.2(mu,sig,p,x0)
  if(is.infinite(est.res3[1])){ 
    num.inva0 = num.inva0 + 1
  }else{
    num.iter0 = num.iter0 + est.res3[2]
  }
}
end.time = Sys.time()
time.us0 = end.time - start.time

start.time = Sys.time()
num.inva1 = 0
num.iter1 = 0
for(i in c(1:n.rep)){
  x0 = runif(1,-4,4)
  est.res3 = estimation_fun2.3(mu,sig,p,x0)
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
  x0 = runif(1,-4,4)
  est.res3 = estimation_fun2.NR(mu,sig,p,x0)
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
  x0 = runif(1,-4,4)
  est.res3 = estimation_fun2.BS(mu,sig,p,x0)
  if(is.infinite(est.res3[1])){ 
    num.inva3 = num.inva3 + 1
  }else{
    num.iter3 = num.iter3 + est.res3[2]
    }
}
end.time = Sys.time()
time.bs = end.time - start.time

time.tra = rbind(c(time.us0, num.inva0,num.iter0/n.rep/(1-num.inva0)),
                 c(time.us, num.inva1,num.iter1/n.rep/(1-num.inva1)), 
                 c(time.nr, (n.rep -num.inva2)/n.rep,num.iter2/(n.rep-num.inva2)), 
                 c(time.bs, num.inva3,num.iter3/n.rep/(1-num.inva3)))
time.tra
