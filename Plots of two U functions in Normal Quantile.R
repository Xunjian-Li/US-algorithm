#####0.05-th quantile of normal distribution
############plots of two U functions###############
library(latex2exp)

plot_points = function(x,y2,xt,y1t,i,adjus1,star1=F,thea=1, poi="topright"){
  if(!star1){
    lines(x,y2,type='l', xlab = expression(theta),ylab = 'b',lwd = 3,col=2, lty=3)
    lines(c(xt,xt),c(y1t,0),type = 'l',lty=2)
    points(xt,0, lwd = 3, col = "blue", pch=19)
    if(thea==1){
      text((xt+adjus1[1]),adjus1[2],as.expression(
        substitute(theta^(alpha),list(alpha = i-1))),cex=1.2)
    }else{
      text((xt+adjus1[1]),adjus1[2],as.expression(
        substitute(x^(alpha),list(alpha = i-1))),cex=1.2)
    }
  }
  if(star1){
    points(xt,0, lwd = 3, col = "green" , pch=4)
    legend((xt+adjus1[1]),adjus1[2], bty='n',legend = c(TeX('$\\textit{theta}^*$'),TeX('$\\textit{x}^*$'))[thea],cex=1.2)
    lenge1 = c(TeX('$\\textit{g(\\theta)}$'), TeX('$\\textit{U(\\theta|\\theta^{(t)})}$'))
    lenge2 = c(TeX('$\\textit{g(x)}$'), TeX('$\\textit{U(x|x^{(t)})}$'))
    legend(poi, inset=.05, legend = c(lenge1,lenge2)[c(1,2)+(thea-1)*2],
           lty=c(1,3),col=c(1,2),
           bty='n',lwd = c(2,2),cex=1.2)
  }
}


############solving a3*xt^3+a2*xt^2+a1*xt + a0 = 0 with closed-form solution###
root_ploynomial_3 = function(a3, a2, a1, a0)
{ # Function name:(a3, a2, a1, a0)################################
  # ---------------------------Aim--------------------------------
  # Calculating the equation: a3*x^3+a2*x^2+a1*x+a0, -inf <x < inf
  # ---------------------------Input------------------------------
  #    a3: the coefficient of x^3
  #    a2: the coefficient of x^2
  #    a1: the coefficient of x
  #    a0: the constant, that is a positive real number
  # ---------------------------Output-----------------------------
  # The real roots 
  ################################################################
  nroot <- function(x,n)
  { # Function for calculating x^(1/n), n = 1,3,5,...     
    abs(x)^(1/n)*sign(x)
  }
  p = (3*a3*a1-a2^2)/(3*a3^2)
  q = (2*a2^3-9*a3*a2*a1+27*a3^2*a0)/(27*a3^3)
  disterm = (q/2)^2+(p/3)^3
  if(disterm>=0){    ## When disterm>0, there exists one real root
    u1 = nroot(-q/2+sqrt((q/2)^2+(p/3)^3),3)
    u2 = nroot(-q/2-sqrt((q/2)^2+(p/3)^3),3)
    t1 = u1+u2
    omega = complex(real=- 1,imaginary=sqrt(3))/2
    t2 = omega*u1 + omega^2*u2
    t3 = omega^2*u1 + omega*u2
    x1 = t1 - a2/(3*a3)
    return(x1)
  }else{           ## When disterm>0, there exist three real roots
    r = sqrt(-(p/3)^3)
    theta = 1/3*acos(-q/(2*r))
    t1 = 2*nroot(r,3)*cos(theta)
    t2 = 2*nroot(r,3)*cos(theta+2/3*pi)
    t3 = 2*nroot(r,3)*cos(theta+4/3*pi)
    x1 = t1 - a2/(3*a3)
    x2 = t2 - a2/(3*a3)
    x3 = t3 - a2/(3*a3)
    return(c(x1,x2,x3))
  }
}

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
#   return(c(xt,k-1))
  return(rbind(xt.save.2,yt.save.2))
}
est.res2 = estimation_fun2.2(mu,sig,p,x0)
est.res = est.res2
x = seq(-2.5,-0.5,0.01)
y1 = p - pnorm(x,mu,sig)
plot(x,y1,type='l', ylim=c(-0.03,0.07),
     xlab = c(expression(theta),'x')[2],ylab = '',cex.lab = 1.5,lwd = 3,col=1, lty=1)
lines(c(-100,100),c(0,0),lty=2 ,lwd = 1)

for(i in c(1:1)){
  xt = est.res[1,i]
  y1t = est.res[2,i]
  b2 = ((x<=xt)-(x>xt))/(sig^3*sqrt(2*pi)*exp(1/2))
  y2 = y1t - dnorm(xt,mu,sig)*(x-xt) + b2/2*(x-xt)^2
  adjus1 = c(0.0,-0.04)
  plot_points(x,y2,xt,y1t,i,adjus1,thea=2)
}
xt = est.res[1,i+1];              y1t = est.res[2,i+1]
plot_points(x,y2,xt,y1t,i+1,adjus1,thea=2)
xt = est.res[1,dim(est.res)[2]];  y1t = est.res[2,dim(est.res)[2]]
adjus1 = c(-0.15,0.1)
plot_points(x,y2,xt,y1t,i+1,adjus1,star1 = T,thea=2)

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
#   return(c(xt,k-1))
  return(rbind(xt.save.2,yt.save.2))
}
est.res3 = estimation_fun2.3(mu,sig,p,x0)
est.res = est.res3
x = seq(-2.5,-0.5,0.01)
y1 = p - pnorm(x,mu,sig)
plot(x,y1,type='l', ylim=c(-0.03,0.07),
     xlab = c(expression(theta),'x')[2],ylab = '',cex.lab = 1.5,lwd = 3,col=1, lty=1)
lines(c(-100,100),c(0,0),lty=2 ,lwd = 1)

for(i in c(1:1)){
  xt = est.res[1,i]
  y1t = est.res[2,i]
  y2 = y1t - dnorm(xt,mu,sig)*(x-xt) + (xt-mu)/sig^2*dnorm(xt,mu,sig)/2*(x-xt)^2 +
    -1/(3*sig^3*sqrt(2*pi)*exp(3/2))*(x-xt)^3
  adjus1 = c(-0.04,-0.005)
  plot_points(x,y2,xt,y1t,i,adjus1,thea=2)
}
xt = est.res[1,i+1];              y1t = est.res[2,i+1]
plot_points(x,y2,xt,y1t,i+1,adjus1,thea=2)
xt = est.res[1,dim(est.res)[2]];  y1t = est.res[2,dim(est.res)[2]]
adjus1 = c(-0.04,0.01)
plot_points(x,y2,xt,y1t,i+1,adjus1,star1 = T,thea=2)

