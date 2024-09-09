estimation_fun1.2 = function(a1,a2,a3,a4,a.max=10,m,x0=0){
  xt.save = c()
  yt.save = c()
  flag = 0
  xt = x0
  y1t = a1*xt^m+a2*xt^2+a3*xt + a4
  # x = seq(-2,2,0.1)
  # plot(-(a1*(-1)^m*x^m+a2*(x)^2-a3*(x) + a4), type='l',ylim = c(-0.6,0.6))
  # plot((a1*(x)^m+a2*(x)^2+a3*(x) + a4), type='l',col='red',ylim = c(-0.6,0.6))
  # lines((a1*(x)^m+a2*(x)^2+a3*(x) + a4), type='l',col='red',ylim = c(-0.6,0.6))
  if(y1t<0){
    a1 = -a1*(-1)^m
    a2 = -a2
    a4 = -a4
    y1t = a1*xt^m+a2*xt^2+a3*xt + a4
    flag = 1
  }
  xt.save[1] = xt
  yt.save[1] = y1t
  k = 1
  esp = 1
  
  while(esp>1e-8){
    b2 = a1*m*(m-1)*a.max^(m-2)
    A1 = b2+a2
    A2 = a1*m*xt^(m-1) -2*b2*xt + a3
    A3 = a1*(1-m)*xt^m + b2*xt^2 + a4
    if(A2^2-4*A1*A3 >= 0){
      xt.save[k+1] = -(A2+sqrt(A2^2-4*A1*A3))/(2*A1)
    }else{
      b2 = 0
      A1 = b2+a2
      A2 = a1*m*xt^(m-1) -2*b2*xt + a3
      A3 = a1*(1-m)*xt^m + b2*xt^2 + a4
      xt.save[k+1] = -(A2+sqrt(A2^2-4*A1*A3))/(2*A1)
    }
    # A1*xt^2+A2*xt+A3
    yt.save[k+1] = a1*xt.save[k+1]^m+a2*xt.save[k+1]^2+a3*xt.save[k+1] + a4
    esp = abs(yt.save[k+1])
    xt = xt.save[k+1]
    y1t = yt.save[k+1]
    k = k+1
  }
  if(flag){
    xt.save = -xt.save
    yt.save = -yt.save
  }
  return(rbind(xt.save,yt.save))
}


x = seq(-3,3,0.01)
y1.1 = a1*(x-xt)^m+a2*(x-xt)^2+a3*(x-xt) + a4
plot(x,y1,type='l', ylim = c(-0.5,1),
     xlab = expression(theta),ylab = '',cex.lab = 1.5,lwd = 3,col=1, lty=1)
lines(c(-3,10),c(0,0),lty=2 ,lwd = 1)


a1 = -2
a4 = 100
est.res0 = estimation_fun1.2(a1,a2,a3,a4,a.max=4,m,x0=4)
est.res = est.res0
# est.res1 = estimation_fun1(a1,a2,a3,a4,m,x0=0.6)
# est.res = est.res1
# est.res1.2 = estimation_fun1(a1,a2,a3,a4,m,x0=0)
# est.res1 = cbind(est.res.Newton0,est.res.Newton1,est.res0,est.res1)
x = seq(0,5,0.01)
y1 = a1*x^m+a2*x^2+a3*x + a4
plot(x,y1,type='l', ylim = c(-200,200),
     xlab = expression(theta),ylab = '',cex.lab = 1.5,lwd = 3,col=1, lty=1)
lines(c(-3,10),c(0,0),lty=2 ,lwd = 1)
for(i in c(1:3)){
  xt = est.res[1,i]
  y1t = est.res[2,i]
  y2 = c()
  b2 = 0
  A1 = b2+a2
  A2 = a1*m*xt^(m-1) -2*b2*xt + a3
  A3 = a1*(1-m)*xt^m + b2*xt^2 + a4
  y2 = A1*x[x<=xt]^2 + A2*x[x<=xt] + A3
  b2 = a1*m*(m-1)*a.max^(m-2)
  A1 = b2+a2
  A2 = a1*m*xt^(m-1) -2*b2*xt + a3
  A3 = a1*(1-m)*xt^m + b2*xt^2 + a4
  
  y2 = c(y2, A1*x[x>xt]^2 + A2*x[x>xt] + A3)
  adjus1 = c(0.012,0.14)
  plot_points(x,y2,xt,y1t,i,thea=1,adjus1)
}
xt = est.res[1,i+1]
y1t = est.res[2,i+1]
plot_points(x,y2,xt,y1t, i+1,thea=1,adjus1)
xt = est.res[1,dim(est.res)[2]]
y1t = est.res[2,dim(est.res)[2]]
adjus1 = c(-0.05, 0)
plot_points(x,y2,xt,y1t,i+1,adjus1,T,thea=1, poi="bottomleft")



############solving a3*xt^3+a2*xt^2+a1*xt + a0 = 0 with closed-form solution###
root_ploynomial_3 = function(a3,a2,a1,a0){
  nroot <- function(x,n){
    abs(x)^(1/n)*sign(x)
  }
  y1t = a3*xt^3+a2*xt^2+a1*xt + a0
  p = (3*a3*a1-a2^2)/(3*a3^2)
  q = (2*a2^3-9*a3*a2*a1+27*a3^2*a0)/(27*a3^3)
  disterm = (q/2)^2+(p/3)^3
  if(disterm>=0){
    u1 = nroot(-q/2+sqrt((q/2)^2+(p/3)^3),3)
    u2 = nroot(-q/2-sqrt((q/2)^2+(p/3)^3),3)
    t1 = u1+u2
    omega = complex(real=- 1,imaginary=sqrt(3))/2
    t2 = omega*u1 + omega^2*u2
    t3 = omega^2*u1 + omega*u2
    
    x1 = t1 - a2/(3*a3)
    return(x1)
  }else{
    r = sqrt(-(p/3)^3)
    theta1 = 1/3*acos(-q/(2*r))
    t1 = 2*nroot(r,3)*cos(theta1)
    t2 = 2*nroot(r,3)*cos(theta1+2/3*pi)
    t3 = 2*nroot(r,3)*cos(theta1+4/3*pi)
    
    x1 = t1 - a2/(3*a3)
    x2 = t2 - a2/(3*a3)
    x3 = t3 - a2/(3*a3)
    return(c(x1,x2,x3))
  }
}

root_ploynomial_3(a3,a2,a1,a0)


