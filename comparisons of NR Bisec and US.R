a3 = -1; a2 = 1; a1 = -1; a0 = 1; m = 3; n.rep = 100000
start.time = Sys.time()
num.inva0 = 0
num.iter0 = 0
for(i in c(1:n.rep)){
  x0 = runif(1,0,2)
  est.res3 = root.UHE2(a3,a2,a1,a0,a.max=2,m,x0=0)
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