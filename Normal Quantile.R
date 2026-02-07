############################################################
#### case 2: Normal quantile equation
#### Solve g(x) = p - Phi(x; mu, sig) = 0
############################################################

## ---------------------------------------------------------
## Base function and derivatives
## g(x) = p - pnorm(x; mu, sig)
## ---------------------------------------------------------
norm_g_int <- function(x, mu, sig, p) {
  z <- (x - mu) / sig
  p * x - (x - mu) * pnorm(z) - sig * dnorm(z)  # + C
}
norm_g <- function(x, mu, sig, p) {
  p - pnorm(x, mean = mu, sd = sig)
}
norm_gp <- function(x, mu, sig, p) {
  -dnorm(x, mean = mu, sd = sig)            # g'(x)
}
norm_gpp <- function(x, mu, sig, p) {
  # g''(x) = (x-mu)/sig^2 * dnorm(x;mu,sig)
  (x - mu) / (sig^2) * dnorm(x, mean = mu, sd = sig)
}
norm_gppp <- function(x, mu, sig, p) {
  # Let z=(x-mu)/sig, phi=dnorm(x;mu,sig)= (1/sig)*phi_std(z)
  # derivative: d/dx [ (x-mu)/sig^2 * phi ] = (1/sig^2)*phi + (x-mu)/sig^2 * phi'
  # phi' = -(x-mu)/sig^2 * phi
  # => g''' = (1/sig^2)*phi - (x-mu)^2/sig^4 * phi = phi*(1/sig^2 - (x-mu)^2/sig^4)
  phi <- dnorm(x, mean = mu, sd = sig)
  phi * (1/(sig^2) - (x - mu)^2/(sig^4))
}

## ---------------------------------------------------------
## Newton-Raphson (no damping)
## ---------------------------------------------------------
root.NR.norm <- function(mu, sig, p, x0,
                         tol = 1e-12, maxit = 200) {
  
  xt.save <- c(); yt.save <- c()
  xt <- x0
  g  <- norm_g(xt, mu, sig, p)
  xt.save[1] <- xt; yt.save[1] <- g
  
  k <- 1
  esp <- abs(g)
  
  while (esp > tol && k < maxit) {
    
    g1 <- norm_gp(xt, mu, sig, p)
    if (!is.finite(g1) || abs(g1) < .Machine$double.eps) return(c(Inf, k-1))
    
    xnew <- xt - g/g1
    
    xt.save[k+1] <- xnew
    gnew <- norm_g(xnew, mu, sig, p)
    yt.save[k+1] <- gnew
    
    esp <- abs(gnew)
    xt <- xnew
    g  <- gnew
    k <- k + 1
    
    if (!is.finite(xt)) return(c(Inf, k-1))
  }
  
  if (esp > tol) xt <- Inf
  
  c(xt, k-1)
}

## ---------------------------------------------------------
## Damped Newton-Raphson with Armijo on |Psi|
##   solve g(x)=0,  g(x)=p - Phi(x;mu,sig)
##   Psi(x)=norm_g_int(x) is an antiderivative of g(x)
##   Armijo test: Psitry >= Psi0 + c * t * grad_Psi_dir
## ---------------------------------------------------------
root.newton.damped.norm <- function(mu, sig, p, x0, 
                                    tol = 1e-12, maxit = 200, 
                                    eta = 0.5, c = 1e-4, maxls = 50) {
  
  xt.save <- c(); yt.save <- c()
  xt <- x0
  g  <- norm_g(xt, mu, sig, p)
  xt.save[1] <- xt; yt.save[1] <- g
  
  k <- 1
  esp <- abs(g)
  
  while (esp > tol && k < maxit) {
    
    g1 <- norm_gp(xt, mu, sig, p)
    if (!is.finite(g1) || abs(g1) < .Machine$double.eps) return(c(Inf, k-1))
    
    # Newton direction
    step <- -g / g1
    
    Psi0 <- norm_g_int(xt, mu, sig, p)
    grad_Psi_dir <- g * step
    
    t <- 1
    ls <- 0
    accepted <- FALSE
    
    while (ls < maxls) {
      xnew <- xt + t * step
      Psitry <- norm_g_int(xnew, mu, sig, p)
      
      if (is.finite(Psitry) && Psitry >= Psi0 + c * t * grad_Psi_dir) {
        accepted <- TRUE
        break
      }
      t <- eta * t
      ls <- ls + 1
    }
    
    if (!accepted) return(c(Inf, k-1))
    
    xt <- xnew
    g  <- norm_g(xt, mu, sig, p)
    
    xt.save[k+1] <- xt
    yt.save[k+1] <- g
    
    esp <- abs(g)
    k <- k + 1
    
    if (!is.finite(xt)) return(c(Inf, k-1))
  }
  
  if (esp > tol) xt <- Inf
  
  c(xt, k-1)
}


## ---------------------------------------------------------
## Bisection method (bracketed root)
## (x0 unused; provide a bracket [L,U])
## ---------------------------------------------------------
root.BS.norm <- function(mu, sig, p,
                         L = -20, U = 20,
                         tol = 1e-12, maxit = 500) {
  
  gL <- norm_g(L, mu, sig, p)
  gU <- norm_g(U, mu, sig, p)
  if (!is.finite(gL) || !is.finite(gU) || gL*gU > 0) return(c(Inf, 0))
  
  xt.save <- c(); yt.save <- c()
  xt <- (L+U)/2
  g  <- norm_g(xt, mu, sig, p)
  xt.save[1] <- xt; yt.save[1] <- g
  
  k <- 1
  esp <- abs(g)
  
  while (esp > tol && k < maxit) {
    
    if (gL*g <= 0) {
      U <- xt; gU <- g
    } else {
      L <- xt; gL <- g
    }
    
    xt <- (L+U)/2
    g  <- norm_g(xt, mu, sig, p)
    
    xt.save[k+1] <- xt
    yt.save[k+1] <- g
    
    esp <- abs(g)
    k <- k + 1
  }
  
  c(xt, k-1)
}

## ---------------------------------------------------------
## Halley's method
## x_{k+1} = x_k - 2 g g' / (2 g'^2 - g g'')
## ---------------------------------------------------------
root.halley.norm <- function(mu, sig, p, x0,
                             tol = 1e-12, maxit = 200) {
  
  xt.save <- c(); yt.save <- c()
  xt <- x0
  g  <- norm_g(xt, mu, sig, p)
  xt.save[1] <- xt; yt.save[1] <- g
  
  k <- 1
  esp <- abs(g)
  
  while (esp > tol && k < maxit) {
    
    g1 <- norm_gp(xt, mu, sig, p)
    g2 <- norm_gpp(xt, mu, sig, p)
    
    denom <- 2*g1^2 - g*g2
    if (!is.finite(denom) || abs(denom) < .Machine$double.eps) return(c(Inf, k-1))
    
    xnew <- xt - (2*g*g1)/denom
    
    xt.save[k+1] <- xnew
    gnew <- norm_g(xnew, mu, sig, p)
    yt.save[k+1] <- gnew
    
    esp <- abs(gnew)
    xt <- xnew
    g  <- gnew
    k <- k + 1
    
    if (!is.finite(xt)) return(c(Inf, k-1))
  }
  
  c(xt, k-1)
}

## ---------------------------------------------------------
## Householder (order-3)
## x_{k+1} = x_k - (6 g g'^2 - 3 g^2 g'') / (6 g'^3 - 6 g g' g'' + g^2 g''')
## ---------------------------------------------------------
root.householder4.norm <- function(mu, sig, p, x0,
                                   tol = 1e-12, maxit = 200) {
  
  xt.save <- c(); yt.save <- c()
  xt <- x0
  g  <- norm_g(xt, mu, sig, p)
  xt.save[1] <- xt; yt.save[1] <- g
  
  k <- 1
  esp <- abs(g)
  
  while (esp > tol && k < maxit) {
    
    g1 <- norm_gp(xt, mu, sig, p)
    g2 <- norm_gpp(xt, mu, sig, p)
    g3 <- norm_gppp(xt, mu, sig, p)
    
    denom <- 6*g1^3 - 6*g*g1*g2 + g^2*g3
    numer <- 6*g*g1^2 - 3*g^2*g2
    if (!is.finite(denom) || abs(denom) < .Machine$double.eps) return(c(Inf, k-1))
    
    xnew <- xt - numer/denom
    
    xt.save[k+1] <- xnew
    gnew <- norm_g(xnew, mu, sig, p)
    yt.save[k+1] <- gnew
    
    esp <- abs(gnew)
    xt <- xnew
    g  <- gnew
    k <- k + 1
    
    if (!is.finite(xt)) return(c(Inf, k-1))
  }
  
  c(xt, k-1)
}


## ---------------------------------------------------------
## Hirano–Murota-style globalized Newton for quantiles
## Solve g(x) = F(x) - p = 0, with g' = f(x)
## Using inner line-search:
##   x_new = x + step, step = -μ * g / g'
##   accept if |g_new| <= (1 - (1-β) μ) |g|
##   else shrink μ := μ/(1+δ)
## ---------------------------------------------------------
root.hirano.norm <- function(mu, sig, p, x0, 
                             tol       = 1e-12, 
                             maxit     = 200, 
                             beta      = 0.5,
                             delta     = 0.5,
                             inner_max = 50) {
  
  if (p <= 0 || p >= 1) stop("Probability p must be in (0, 1)")
  if (sig <= 0) stop("Sigma must be positive")
  
  n_degree <- 2
  z <- (x0 - mu) / sig
  
  for (iter in 1:maxit) {
    z_eval <- Re(z) 
    pdf_val <- dnorm(z_eval)
    
    if (pdf_val < .Machine$double.eps) {
      return(c(root = mu + z_eval * sig, iter = maxit))
    }
    
    cdf_val <- pnorm(z_eval) - p
    
    a <- complex(length = n_degree + 1)
    a[1] <- cdf_val             
    a[2] <- pdf_val
    a[3] <- (-z_eval * pdf_val) / 2
    
    if (abs(Mod(a[1])) < tol) {
      return(c(root = mu + Re(z) * sig, iter = iter - 1))
    }
    
    mu_param <- 1
    success <- FALSE
    
    for (j in 1:inner_max) {
      candidates <- complex(length = n_degree)
      
      for (k in 1:n_degree) {
        denom <- a[k+1]
        if (Mod(denom) < 1e-100) { 
          candidates[k] <- complex(real=Inf, imaginary=Inf)
        } else {
          term <- -mu_param * a[1] / denom
          candidates[k] <- term^(1/k)
        }
      }
      
      moduli <- Mod(candidates)
      if (all(is.infinite(moduli))) break 
      
      m <- which.min(moduli)
      zeta_m <- candidates[m]
      
      poly_val <- complex(real = 0, imaginary = 0)
      for (k in 0:n_degree) {
        poly_val <- poly_val + a[k+1] * (zeta_m^k)
      }
      
      lhs <- Mod(poly_val)
      rhs <- (1 - (1 - beta) * mu_param) * Mod(a[1])
      
      if (is.finite(lhs) && is.finite(rhs) && lhs <= rhs) {
        z <- z + zeta_m
        success <- TRUE
        break
      } else {
        mu_param <- mu_param / (1 + delta)
      }
    }
    
    if (!success) {
      return(c(root = mu + Re(z) * sig, iter = maxit)) 
    }
  }
  
  return(c(root = mu + Re(z) * sig, iter = maxit))
}


## ---------------------------------------------------------
## Your US variants (SLUB-like and TLB-like)
## Keep your formulas, but harden:
##  - maxit
##  - denom checks
##  - return Inf on failure
## ---------------------------------------------------------

## US-2 (your estimation_fun2.2: SLUB-like quadratic surrogate step)
root.US2.norm <- function(mu, sig, p, x0,
                          tol = 1e-12, maxit = 200) {
  
  xt.save <- c(); yt.save <- c()
  xt <- x0
  g  <- norm_g(xt, mu, sig, p)
  xt.save[1] <- xt; yt.save[1] <- g
  
  k <- 1
  esp <- abs(g)
  
  while (esp > tol && k < maxit) {
    
    A2 <- -dnorm(xt, mu, sig)           # your A2
    A3 <- g                             # your A3
    if (!is.finite(A2) || abs(A2) < .Machine$double.eps) return(c(Inf, k-1))
    
    # your b2 choice (sign depends on g)
    if (g > 0) {
      b2 <- -1/(sig^3*sqrt(2*pi)*exp(1/2))
    } else {
      b2 <-  1/(sig^3*sqrt(2*pi)*exp(1/2))
    }
    A1 <- b2/2
    
    disc <- A2^2 - 4*A1*A3
    if (!is.finite(disc) || disc < 0 || abs(A1) < .Machine$double.eps) return(c(Inf, k-1))
    
    xnew <- xt - (A2 + sqrt(disc))/(2*A1)
    
    xt.save[k+1] <- xnew
    gnew <- norm_g(xnew, mu, sig, p)
    yt.save[k+1] <- gnew
    
    esp <- abs(gnew)
    xt <- xnew
    g  <- gnew
    k <- k + 1
    
    if (!is.finite(xt)) return(c(Inf, k-1))
  }
  
  c(xt, k-1)
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

## US-3 (your estimation_fun2.3: cubic surrogate step on increment z = x_{k+1}-x_k)
root.US3.norm <- function(mu, sig, p, x0,
                          tol = 1e-12, maxit = 200) {
  
  xt.save <- c(); yt.save <- c()
  xt <- x0
  g  <- norm_g(xt, mu, sig, p)
  xt.save[1] <- xt; yt.save[1] <- g
  
  k <- 1
  esp <- abs(g)
  
  while (esp > tol && k < maxit) {
    
    g1 <- norm_gp(xt, mu, sig, p)
    g2 <- norm_gpp(xt, mu, sig, p)
    
    # your coefficients
    a3 <- -1/(3*sig^3*sqrt(2*pi)*exp(3/2))
    a2 <- g2/2
    a1 <- g1
    a0 <- g
    
    rr = root_ploynomial_3(a3,a2,a1,a0)
    
    # your selection rule based on sign of g
    if (length(rr) == 1) {
      z <- rr[1]
    } else {
      if (g > 0) {
        cand <- rr[rr >= 0]
        if (length(cand) == 0) return(c(Inf, k-1))
        z <- min(cand)
      } else {
        cand <- rr[rr <= 0]
        if (length(cand) == 0) return(c(Inf, k-1))
        z <- max(cand)
      }
    }
    
    xnew <- xt + z
    
    xt.save[k+1] <- xnew
    gnew <- norm_g(xnew, mu, sig, p)
    yt.save[k+1] <- gnew
    
    esp <- abs(gnew)
    xt <- xnew
    g  <- gnew
    k <- k + 1
    
    if (!is.finite(xt)) return(c(Inf, k-1))
  }
  
  c(xt, k-1)
}

############################################################
#### Monte Carlo comparisons (shared x0, warmup, time + iter)
############################################################
compare_norm_quantile <- function(mu, sig, p,
                                  n.rep = 100000,
                                  xlow = -4, xhigh = 4,
                                  seed = NULL,
                                  warmup = 5000) {
  
  if (!is.null(seed)) set.seed(seed)
  
  ## shared initial values
  x0_vec <- runif(n.rep, xlow, xhigh)
  
  methods <- list(
    Newton       = function(x0) root.NR.norm(mu,sig,p,x0=x0),
    Bisection    = function(x0) root.BS.norm(mu,sig,p, L=-20, U=20),
    Halley       = function(x0) root.halley.norm(mu,sig,p,x0=x0),
    Householder4 = function(x0) root.householder4.norm(mu,sig,p,x0=x0),
    DampedNewton = function(x0) root.newton.damped.norm(mu,sig,p,x0=x0),
    Hirano       = function(x0) root.hirano.norm(mu,sig,p,x0=x0),
    US2          = function(x0) root.US2.norm(mu,sig,p,x0=x0),
    US3          = function(x0) root.US3.norm(mu,sig,p,x0=x0)
  )
  
  ## -------- warmup --------
  n.warm <- min(warmup, n.rep)
  for (j in seq_along(methods)) {
    for (i in seq_len(n.warm)) {
      invisible(methods[[j]](x0_vec[i]))
    }
  }
  
  ## time_success_sec = total time over successful runs only
  ## avg_time_success = average time per successful run
  out <- matrix(NA_real_, nrow = length(methods), ncol = 5)
  rownames(out) <- names(methods)
  colnames(out) <- c("time_success_sec", "avg_time_success_us",
                     "num_invalid", "num_success", "avg_iter_success")
  
  for (j in seq_along(methods)) {
    
    num.inva <- 0L
    num.iter <- 0.0
    num.succ <- 0L
    time.ok  <- 0.0  # seconds, successful runs only
    
    for (i in seq_len(n.rep)) {
      t0 <- proc.time()[3]
      est <- methods[[j]](x0_vec[i])
      t1 <- proc.time()[3]
      
      xhat <- est[1]
      it   <- est[2]
      
      ok <- is.finite(xhat) && !is.nan(xhat)
      if (!ok) {
        num.inva <- num.inva + 1L
        next
      }
      
      ## success: accumulate time and iterations
      num.succ <- num.succ + 1L
      num.iter <- num.iter + it
      time.ok  <- time.ok + (t1 - t0)
    }
    
    avg_iter <- if (num.succ > 0) num.iter/num.succ else NA_real_
    avg_time_us <- if (num.succ > 0) (time.ok/num.succ)*1e6 else NA_real_
    
    out[j, ] <- c(time.ok, avg_time_us, num.inva, num.succ, avg_iter)
  }
  
  out
}

############################################################
#### Example run (your scenarios)
############################################################
mu <- 2; sig <- 1; p <- 0.9     # x* ~ 3.281552
n.rep <- 10000

time.tra1 <- compare_norm_quantile(mu, sig, p,
                                  n.rep = n.rep,
                                  xlow = -4, xhigh = 4,
                                  seed = 1)

mu <- 2; sig <- 1; p <- 0.1     # x* ~ 3.281552
n.rep <- 10000

time.tra2 <- compare_norm_quantile(mu, sig, p,
                                  n.rep = n.rep,
                                  xlow = -4, xhigh = 4,
                                  seed = 1)

results1 = cbind(time.tra1[,c(2,4,5)], time.tra2[,c(2,4,5)])

library(xtable)

tab <- as.data.frame(results1)
tab[, c(2,5)] = tab[, c(2,5)]/100
tab[, c(2,5)] <- lapply(tab[, c(2,5)], function(x) sprintf("%.2f\\%%", x ))

xt <- xtable(tab,
             caption = "Performance comparison of solvers",
             label   = "tab:perf")

print(xt,
      include.rownames = TRUE,
      comment = FALSE,
      sanitize.text.function = identity)

