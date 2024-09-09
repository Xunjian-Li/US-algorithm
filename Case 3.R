#####case 3: Solving an equation with multiple roots
############plot of multiroot###############
library(latex2exp)

a2 = -0.5
a1 = -2
a0 = 1
root.sin = function(a2,a1,a0,x0=0){
  xt.save = c(); yt.save = c()
  xt = x0; y1t = a2*xt+a1*sin(xt)+a0
  xt.save[1] = xt; yt.save[1] = y1t
  k = 1; esp = 1
  while(esp>1e-8){
    if(y1t>0){b1 = 0.4}else{
      a2 = - a2
      a1 = - a1
      a0 = -a0
      b1 = 2/3
      y1t = a2*xt+a1*sin(xt)+a0
      }
    xt.save[k+1] = b1*y1t + xt
    yt.save[k+1] = a2*xt.save[k+1]+a1*sin(xt.save[k+1])+a0
    esp = abs(yt.save[k+1])
    xt = xt.save[k+1]
    y1t = yt.save[k+1]
    k = k+1
  }
  return(rbind(xt.save,yt.save))
}
res1 = root.sin(a2,a1,a0,x0=0)
res2 = root.sin(a2,a1,a0,x0=4.090497e-01+0.00001)
res3 = root.sin(a2,a1,a0,x0=3.535612e+00+0.00001)

num1 = length(res1[1,])
num2 = length(res2[1,])
num3 = length(res3[1,])

par(mfrow = c(1,3))
plot(c(res1[1,]),type='b',ylab = TeX('$\\textit{X}_^{(t)}$'),xlab = 't',main = '(b1)')
points(num1,res1[1,num1], lwd = 3, col = "green" , pch=4)
legend(num1-1.5,res1[1,num1], bty='n',legend = TeX('$\\textit{X}_1^*$'),cex=1.2)

plot(c(res2[1,]),type='b',ylab = TeX('$\\textit{X}_^{(t)}$'),xlab = 't',main = '(b2)')
points(num2,res2[1,num2], lwd = 3, col = "green" , pch=4)
legend(num2-8,res2[1,num2], bty='n',legend = TeX('$\\textit{X}_2^*$'),cex=1.2)

plot(c(res3[1,]),type='b',ylab = TeX('$\\textit{X}_^{(t)}$'),xlab = 't',main = '(b3)')
points(num3,res3[1,num3], lwd = 3, col = "green" , pch=4)
legend(num3-8,res3[1,num3], bty='n',legend = TeX('$\\textit{X}_3^*$'),cex=1.2)


plot(c(res1[1,],res2[1,],res3[1,]),type='b',ylab = TeX('$\\textit{X}_^{(t)}$'),xlab = 't')
points(num1,res1[1,num1], lwd = 3, col = "green" , pch=4)
legend(num1-4,res1[1,num1], bty='n',legend = TeX('$\\textit{X}_1^*$'),cex=1.2)

points(num2+num1,res2[1,num2], lwd = 3, col = "green" , pch=4)
legend(num2+num1-4,res2[1,num2]+0.6, bty='n',legend = TeX('$\\textit{X}_2^*$'),cex=1.2)

points(num3+num2+num1,res3[1,num3], lwd = 3, col = "green" , pch=4)
legend(num3+num2+num1-5,res3[1,num3], bty='n',legend = TeX('$\\textit{X}_3^*$'),cex=1.2)





plot(c(res1[1,],res2[1,],res3[1,]),type='b',ylab = TeX('$\\textit{X}_^{(t)}$'),xlab = 't')
points(num1,res1[1,num1], lwd = 3, col = "green" , pch=4)
legend(num1-4,res1[1,num1]+0.6, bty='n',legend = TeX('$\\textit{X}_1^*$'),cex=1.2)

points(num2+num1,res2[1,num2], lwd = 3, col = "green" , pch=4)
legend(num2+num1-4,res2[1,num2]+0.6, bty='n',legend = TeX('$\\textit{X}_2^*$'),cex=1.2)

points(num3+num2+num1,res3[1,num3], lwd = 3, col = "green" , pch=4)
legend(num3+num2+num1-5,res3[1,num3], bty='n',legend = TeX('$\\textit{X}_3^*$'),cex=1.2)
