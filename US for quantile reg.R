####### US for quantile regression#######
x = seq(-0.3,0.3,0.01)
y = abs(x)/2

######################################################
xt = -0.1
q = 1/2
y1 = (x^2/abs(xt) + (4*q-2)*x+abs(xt))/4
y2 = (2*x-xt)^2/8/abs(xt) + 3*abs(xt)/8

plot(x,y,type = 'l')
lines(x,y1,type = 'l',col='red')
lines(x,y2,type = 'l',col='green')


######################################################
xt = 0.1
q = 1/3
y = c(-(1-q)*x[x<=0],q*x[x>0])
y1 = (x^2/abs(xt) + (4*q-2)*x+abs(xt))/4
y2 = (2*x-xt)^2/8/abs(xt) - (1/2-q)*x + 3*abs(xt)/8

plot(x,y,type = 'l')
lines(x,y1,type = 'l',col='red')
lines(x,y2,type = 'l',col='green')