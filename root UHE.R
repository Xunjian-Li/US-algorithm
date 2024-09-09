root.UHE = function(a3, a2, a1, a0, a.max=10, m, x0 = 0)
{ # Function name:(a3, a2, a1, a0, a.max=10, m, x0 = 0)###########
  # ---------------------------Aim--------------------------------
  # Calculating the equation: a3*x^m+a2*x^2+a1*x+a0, 0 < x < a.max
  # ---------------------------Input------------------------------
  #    a3: the coefficient of x^m
  #    a2: the coefficient of x^2
  #    a1: the coefficient of x
  #    a0: the constant, that is a positive real number
  # a.max: the maximum of x
  #     m: the highest order of the equation
  #    x0: an initial value and the default value is 0
  # ---------------------------Output-----------------------------
  # The unique root in (0,a.max)
  ################################################################
  xt.save = c(); yt.save = c(); 
  xt = x0; y1t = a3*xt^m+a2*xt^2+a1*xt+a0
  xt.save[1] = xt; yt.save[1] = y1t
  k = 1; esp = 1
  while(esp>1e-8){
    if(a3<0){                ### When a3<0, we apply SLUF method
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
  # return(c(xt,k-1))
  return(rbind(xt.save,yt.save))
}
# --------Example: x^3-3x^2-x+1 = 0----------------------------
# -------------------------------------------------------------
a3 = 1; a2 = -3; a1 = -1; a0 = 1; m = 3
root.UHE(a3,a2,a1,a0,a.max=10,m=3,x0=0)
x = seq(0,2,0.01)
plot(x, a3*x^m+a2*x^2+a1*x + a0,type='l',ylab = "U function")
lines(c(-10,10),c(0,0),lty=2)

root.UHE.revise = function(a3,a2,a1,a0,a.max=10,m,x0=0){
  flag = 0
  xt = x0
  y1t = a3*xt^m+a2*xt^2+a1*xt + a0
  if(y1t<0){
    a3 = -a3*(-1)^m; a2 = -a2; a0 = -a0; y1t = a3*xt^m+a2*xt^2+a1*xt+a0
    flag = 1
  }
  if(flag){
    -root.UHE(a3,a2,a1,a0,a.max=10,m=3,x0=0)
  }else{
    root.UHE(a3,a2,a1,a0,a.max=10,m=3,x0=0)
  }
}

a3 = -1; a2 = 1; a1 = -1; a0 = -1; m = 3
est.res0 = root.UHE.revise(a3,a2,a1,a0,a.max=4,m,x0=4)
est.res = est.res0

root.UHE.revise(a3,a2,a1,a0,a.max,m=3,x0=0)

# xt = est.res[1,i+1]
# y1t = est.res[2,i+1]
# plot_points(x,y2,xt,y1t, i+1,thea=1,adjus1)
# xt = est.res[1,dim(est.res)[2]]
# y1t = est.res[2,dim(est.res)[2]]
# adjus1 = c(-0.05, 0)
# plot_points(x,y2,xt,y1t,i+1,adjus1,T,thea=1, poi="bottomleft")




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

root_ploynomial_3(1,-3,1,1)
