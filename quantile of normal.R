#####case 2: 0.5-th quantile of normal distribution
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

mu = 1; sig = 1; p = 0.05; x0 = 0
x0 = -1.5
###################FLB algorithm##################
estimation_fun2 = function(mu,sig,p,x0=0){
  xt.save.2 = c(); yt.save.2 = c(); xt = x0
  y1t = p - pnorm(xt,mu,sig)
  xt.save.2[1] = xt; yt.save.2[1] = y1t
  k = 1; esp = 1
  while(esp>1e-12){
    xt.save.2[k+1] = xt + y1t/dnorm(mu,mu,sig)
    yt.save.2[k+1] = p - pnorm(xt.save.2[k+1],mu,sig)
    esp = abs(yt.save.2[k+1])
    xt = xt.save.2[k+1]
    y1t = yt.save.2[k+1]
    k = k+1
  }
  return(rbind(xt.save.2,yt.save.2))
}
est.res1 = estimation_fun2(mu,sig,p,x0)
est.res = est.res1
x = seq(-1.5,0,0.01)
y1 = p - pnorm(x,mu,sig)
plot(x,y1,type='l', ylim=c(-0.15,0.15),
     xlab = c(expression(theta),'x')[2],ylab = '',cex.lab = 1.5,lwd = 3,col=1, lty=1)
lines(c(-100,100),c(0,0),lty=2 ,lwd = 1)

for(i in c(1:2)){
  xt = est.res[1,i]
  y1t = est.res[2,i]
  y2 = y1t - dnorm(mu,mu,sig)*(x-xt)
  adjus1 = c(0.0,-0.04)
  plot_points(x,y2,xt,y1t,i,adjus1,thea=2)
}
xt = est.res[1,i+1];              y1t = est.res[2,i+1]
plot_points(x,y2,xt,y1t,i+1,adjus1,thea=2)
xt = est.res[1,dim(est.res)[2]];  y1t = est.res[2,dim(est.res)[2]]
adjus1 = c(-0.15,0.1)
plot_points(x,y2,xt,y1t,i+1,adjus1,star1 = T,thea=2)
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
  # return(c(xt,k-1))
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
  return(c(xt,k-1))
  # return(rbind(xt.save.2,yt.save.2))
}
est.res3 = estimation_fun2.3(mu,sig,p,x0)
est.res = est.res3
x = seq(-1.5,0,0.01)
y1 = p - pnorm(x,mu,sig)
plot(x,y1,type='l', ylim=c(-0.05,0.05),
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
  # return(rbind(xt.save.2,yt.save.2))
}
est.res4 = estimation_fun2.NR(mu,sig,p,x0)
est.res = est.res4
x = seq(-1.5,0,0.01)
y1 = p - pnorm(x, mu, sig)
plot(x,y1,type='l', ylim=c(-0.15,0.15),
     xlab = c(expression(theta),'x')[2],ylab = '',cex.lab = 1.5,lwd = 3,col=1, lty=1)
lines(c(-100,100),c(0,0),lty=2 ,lwd = 1)

for(i in c(1:2)){
  xt = est.res[1,i]
  y1t = est.res[2,i]
  y2 = y1t - dnorm(xt,mu,sig)*(x-xt)
  adjus1 = c(0.0,-0.04)
  plot_points(x,y2,xt,y1t,i,adjus1,thea=2)
}
xt = est.res[1,i+1];              y1t = est.res[2,i+1]
plot_points(x,y2,xt,y1t,i+1,adjus1,thea=2)
xt = est.res[1,dim(est.res)[2]];  y1t = est.res[2,dim(est.res)[2]]
adjus1 = c(-0.15,0.1)
plot_points(x,y2,xt,y1t,i+1,adjus1,star1 = T,thea=2)


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
  # return(rbind(xt.save.2,yt.save.2))
}
est.res5 = estimation_fun2.BS(mu,sig,p,x0)
est.res = est.res5

est.res1
est.res2
est.res3
est.res4
est.res5
