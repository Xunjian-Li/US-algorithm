######### plots of UC functions
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

#######case 3: trigonometric function
x = seq(-2,1,0.1); y1 = cos(pi/2*x) - x; xt = -1; 
y1t = cos(pi/2*xt) - xt; y2 = y1t - (pi/2+1)*(x-xt)

estimation_fun3 = function(x0){
  xt.save = c(); yt.save = c(); xt = x0
  y1t = cos(pi/2*xt) - xt
  xt.save[1] = xt; yt.save[1] = y1t
  k = 1; esp = 1
  while(esp>1e-8){
    xt.save[k+1] = y1t/(pi/2+1) + xt
    yt.save[k+1] = cos(pi/2*xt.save[k+1]) - xt.save[k+1]
    esp = abs(yt.save[k+1])
    xt = xt.save[k+1]
    y1t = yt.save[k+1]
    k = k+1
  }
  return(rbind(xt.save,yt.save))
}

est.res3 = estimation_fun3(-1)
est.res = est.res3

x = seq(-2,1,0.1)
y1 = cos(pi/2*x) - x
plot(x,y1,type='l', ylim=c(-1,2),
     xlab = 'x',ylab = '',cex.lab = 1.5,lwd = 3,col=1, lty=1)
lines(c(-100,100),c(0,0),lty=1 ,lwd = 0.5)
for(i in c(1:2)){
  xt = est.res[1,i]
  y1t = est.res[2,i]
  y2 = y1t - (pi/2+1)*(x-xt)
  adjus1 = c(-0.05,-0.15)
  plot_points(x,y2,xt,y1t,i,thea=2,adjus1)
}

xt = est.res[1,i+1]
y1t = est.res[2,i+1]
plot_points(x,y2,xt,y1t,i+1,thea=2,adjus1)

xt = est.res[1,dim(est.res)[2]]
y1t = est.res[2,dim(est.res)[2]]
adjus1 = c(-0.1, 0.3)
plot_points(x,y2,xt,y1t,i,thea=2,adjus1,T)

est.res3.1 = estimation_fun3(3)
est.res = est.res3.1

x = seq(-0,4,0.1)
y1 = cos(pi/2*x) - x
plot(x,y1,type='l', ylim=c(-3.5,1.5),
     xlab = 'x',ylab = '',cex.lab = 1.5,lwd = 3,col=1, lty=1)
lines(c(-100,100),c(0,0),lty=1 ,lwd = 0.5)
for(i in c(1:2)){
  xt = est.res[1,i]
  y1t = est.res[2,i]
  y2 = y1t - (pi/2+1)*(x-xt)
  adjus1 = c(0.08,0.15)
  plot_points(x,y2,xt,y1t,i,thea=2,adjus1)
}

xt = est.res[1,i+1]
y1t = est.res[2,i+1]
plot_points(x,y2,xt,y1t,i+1,thea=2,adjus1)

xt = est.res[1,dim(est.res)[2]]
y1t = est.res[2,dim(est.res)[2]]
adjus1 = c(-0.2, 0.0)
plot_points(x,y2,xt,y1t,i,thea=2,adjus1,T)

