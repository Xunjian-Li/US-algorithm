############################################################
#### case 1: Root of a high-order polynomial equation
############################################################

## ---------------------------------------------------------
## Base polynomial and derivatives
## g(x) = a3*x^m + a2*x^2 + a1*x + a0
## ---------------------------------------------------------
poly_F <- function(x, a3, a2, a1, a0, m){
  a3/(m+1)*x^(m+1) + a2/3*x^3 + a1/2*x^2 + a0*x
}
poly_f <- function(x, a3, a2, a1, a0, m){
  a3*x^m + a2*x^2 + a1*x + a0
}
poly_fp <- function(x, a3, a2, a1, a0, m){
  m*a3*x^(m-1) + 2*a2*x + a1
}
poly_fpp <- function(x, a3, a2, a1, a0, m){
  m*(m-1)*a3*x^(m-2) + 2*a2
}
poly_fppp <- function(x, a3, a2, a1, a0, m){
  m*(m-1)*(m-2)*a3*x^(m-3)
}

## ---------------------------------------------------------
## Halley's method
## ---------------------------------------------------------
root.halley = function(a3, a2, a1, a0, m, x0,
                       tol = 1e-12, maxit = 200) {
  
  xt.save.2 = c(); yt.save.2 = c()
  xt = x0
  g  = poly_f(xt, a3, a2, a1, a0, m)
  xt.save.2[1] = xt; yt.save.2[1] = g
  
  k = 1
  esp = abs(g)
  
  while (esp > tol && k < maxit) {
    
    g1 = poly_fp(xt, a3, a2, a1, a0, m)
    g2 = poly_fpp(xt, a3, a2, a1, a0, m)
    
    denom = 2*g1^2 - g*g2
    if (!is.finite(denom) || abs(denom) < .Machine$double.eps) {
      return(c(Inf, k-1))
    }
    
    xnew = xt - (2*g*g1)/denom
    
    xt.save.2[k+1] = xnew
    gnew = poly_f(xnew, a3, a2, a1, a0, m)
    yt.save.2[k+1] = gnew
    
    esp = abs(gnew)
    xt = xnew
    g  = gnew
    k = k + 1
    
    if (!is.finite(xt)) {
      return(c(Inf, k-1))
    }
  }
  
  return(c(xt, k-1))
}

## ---------------------------------------------------------
## Householder (order-3) method
## x_{k+1} = x_k - (6 g g'^2 - 3 g^2 g'') / (6 g'^3 - 6 g g' g'' + g^2 g''')
## ---------------------------------------------------------
root.householder4 = function(a3, a2, a1, a0, m, x0,
                             tol = 1e-12, maxit = 200) {
  
  xt.save.2 = c(); yt.save.2 = c()
  xt = x0
  g  = poly_f(xt, a3, a2, a1, a0, m)
  xt.save.2[1] = xt; yt.save.2[1] = g
  
  k = 1
  esp = abs(g)
  
  while (esp > tol && k < maxit) {
    
    g1 = poly_fp(xt, a3, a2, a1, a0, m)
    g2 = poly_fpp(xt, a3, a2, a1, a0, m)
    g3 = poly_fppp(xt, a3, a2, a1, a0, m)
    
    denom = 6*g1^3 - 6*g*g1*g2 + g^2*g3
    numer = 6*g*g1^2 - 3*g^2*g2
    
    if (!is.finite(denom) || abs(denom) < .Machine$double.eps) {
      return(c(Inf, k-1))
    }
    
    xnew = xt - numer/denom
    
    xt.save.2[k+1] = xnew
    gnew = poly_f(xnew, a3, a2, a1, a0, m)
    yt.save.2[k+1] = gnew
    
    esp = abs(gnew)
    xt = xnew
    g  = gnew
    k = k + 1
    
    if (!is.finite(xt)) {
      return(c(Inf, k-1))
    }
  }
  
  return(c(xt, k-1))
}

## ---------------------------------------------------------
## Damped Newton-Raphson with backtracking on |g|
## ---------------------------------------------------------
root.newton.damped <- function(a3, a2, a1, a0, m, x0,
                                 tol = 1e-12, maxit = 200,
                                 eta = 0.5, c = 1e-4, maxls = 50){
  
  xt <- x0
  g  <- poly_f(xt, a3, a2, a1, a0, m)
  Fx <- poly_F(xt, a3, a2, a1, a0, m)
  
  k <- 0
  while (abs(g) > tol && k < maxit){
    
    g1 <- poly_fp(xt, a3, a2, a1, a0, m)
    if (!is.finite(g1) || abs(g1) < .Machine$double.eps){
      return(c(Inf, k))
    }
    
    ## Newton direction for solving g(x)=0
    p <- -g / g1
    
    ## Ensure descent direction for minimizing F: g*p < 0
    if (!is.finite(p) || g*p >= 0){
      p <- -g   # steepest descent for F since F'(x)=g(x)
    }
    
    t <- 1
    for (ls in 1:maxls){
      xtry <- xt + t*p
      Ftry <- poly_F(xtry, a3, a2, a1, a0, m)
      
      ## Armijo sufficient decrease on F:
      ## F(x + t p) <= F(x) + c t g p
      if (is.finite(Ftry) && Ftry <= Fx + c*t*g*p) break
      t <- eta*t
    }
    
    xt <- xt + t*p
    if (!is.finite(xt)) return(c(Inf, k))
    
    g  <- poly_f(xt, a3, a2, a1, a0, m)
    Fx <- poly_F(xt, a3, a2, a1, a0, m)
    
    k <- k + 1
  }
  
  if (esp > tol) xt <- Inf
  
  return(c(xt, k))
}


## ---------------------------------------------------------
## Hirano–Murota-style globalized iteration for polynomial root
## Solve g(x) = a3*x^m + a2*x^2 + a1*x + a0 = 0
## Using inner line-search:
##   x_new = x + ζ_m(μ),  ζ_k(μ) = (-μ a0 / a_k)^(1/k), k = 1,2
##   accept if |∑_{k=0}^2 a_k ζ_m^k| <= (1 - (1-β) μ) |a0|
##   else shrink μ := μ/(1+δ)
## ---------------------------------------------------------
root.hirano <- function(a3, a2, a1, a0, m, x0,
                             tol       = 1e-12,
                             maxit     = 200,
                             beta      = 0.5,   # 0 < beta < 1
                             delta     = 0.5,   # delta > 0
                             inner_max = 50) {  # max inner line-search steps
  
  if (m < 2) stop("m must be >= 2 for second derivative to make sense")
  
  n_degree <- 2L
  x <- x0
  g0 <- poly_f(x, a3, a2, a1, a0, m)
  
  x <- as.complex(x)
  g0 <- as.complex(g0)
  
  for (iter in 1:maxit) {
    a <- complex(length = n_degree + 1L)
    a[1] <- g0
    
    x_eval <- Re(x) 
    g1 <- poly_fp (x_eval, a3, a2, a1, a0, m)
    g2 <- poly_fpp(x_eval, a3, a2, a1, a0, m)
    
    a[2] <- g1
    a[3] <- g2 / 2
    
    if (Mod(a[1]) < tol) {
      return(c(root = Re(x), iter = iter - 1L))
    }
    
    ## μ 初始值
    mu_param <- 1
    success  <- FALSE
    
    ## ---------- inner Hirano–Murota line search ----------
    for (j in 1:inner_max) {
      candidates <- complex(length = n_degree)
      eps <- 1e-100
      
      for (k in 1:n_degree) {
        denom <- a[k + 1L]
        if (!is.finite(Re(denom)) || Mod(denom) < eps) {
          candidates[k] <- complex(real = Inf, imaginary = Inf)
        } else {
          term <- -mu_param * a[1] / denom
          candidates[k] <- term^(1 / k)
        }
      }
      
      moduli <- Mod(candidates)
      if (all(is.infinite(moduli))) break
      m_idx  <- which.min(moduli)
      zeta_m <- candidates[m_idx]
      
      poly_val <- complex(real = 0, imaginary = 0)
      for (k in 0:n_degree) {
        poly_val <- poly_val + a[k + 1L] * (zeta_m^k)
      }
      
      lhs <- Mod(poly_val)
      rhs <- (1 - (1 - beta) * mu_param) * Mod(a[1])
      
      if (is.finite(lhs) && is.finite(rhs) && lhs <= rhs) {
        x <- x + zeta_m
        
        g0 <- poly_f(Re(x), a3, a2, a1, a0, m)
        g0 <- as.complex(g0)
        
        success <- TRUE
        break
      } else {
        mu_param <- mu_param / (1 + delta)
      }
    }
    
    if (!success) {
      return(c(root = Re(x), iter = maxit))
    }
  }
  
  return(c(root = Re(x), iter = maxit))
}

## ---------------------------------------------------------
## Your original SLUB / TLF method for unique higher-order equation
## (modified: Delta.det < 0 -> return Inf rather than stop)
## ---------------------------------------------------------
root.UHE = function(a3, a2, a1, a0, a.max = 10, m, x0,
                    tol = 1e-12) {
  
  xt.save = c(); yt.save = c()
  xt = x0
  y1t = poly_f(xt, a3, a2, a1, a0, m)
  xt.save[1] = xt; yt.save[1] = y1t
  
  k = 1
  esp = 1
  
  while (esp > tol) {
    
    if (a3 < 0) {  # SLUB
      if (y1t > 0) { b2 = a3*m*(m-1)*a.max^(m-2) } else { b2 = 0 }
      A2 = b2/2 + a2
      A1 = a3*m*xt^(m-1) - b2*xt + a1
      A0 = a3*(1-m)*xt^m + b2/2*xt^2 + a0
    } else {       # TLF
      A2 = a3*m*(m-1)/2*xt^(m-2) + a2
      A1 = a3*m*(2-m)*xt^(m-1) + a1
      A0 = a3*(m-1)*(m/2-1)*xt^m + a0
    }
    
    Delta.det = A1^2 - 4*A2*A0
    
    if (Delta.det >= 0 && is.finite(Delta.det)) {
      xt.save[k+1] = -(A1 + sqrt(Delta.det))/(2*A2)
      yt.save[k+1] = poly_f(xt.save[k+1], a3, a2, a1, a0, m)
      esp = abs(yt.save[k+1])
      xt = xt.save[k+1]
      y1t = yt.save[k+1]
      k = k + 1
      if (!is.finite(xt)) return(c(Inf, k-1))
    } else {
      return(c(Inf, k-1))
    }
  }
  
  return(c(xt, k-1))
}

## ---------------------------------------------------------
## TLB / TLF method (your root.UHE3)
## (modified: Delta.det < 0 -> return Inf rather than stop)
## ---------------------------------------------------------
root.UHE3 = function(a3, a2, a1, a0, a.max = 10, m, x0,
                     tol = 1e-12) {
  
  xt.save = c(); yt.save = c()
  xt = x0
  y1t = poly_f(xt, a3, a2, a1, a0, m)
  xt.save[1] = xt; yt.save[1] = y1t
  
  k = 1
  esp = 1
  
  while (esp > tol) {
    
    A2 = a3*m*(m-1)/2*xt^(m-2) + a2
    A1 = a3*m*(2-m)*xt^(m-1) + a1
    A0 = a3*(m-1)*(m/2-1)*xt^m + a0
    
    Delta.det = A1^2 - 4*A2*A0
    
    if (Delta.det >= 0 && is.finite(Delta.det)) {
      xt.save[k+1] = -(A1 + sqrt(Delta.det))/(2*A2)
      yt.save[k+1] = poly_f(xt.save[k+1], a3, a2, a1, a0, m)
      esp = abs(yt.save[k+1])
      xt = xt.save[k+1]
      y1t = yt.save[k+1]
      k = k + 1
      if (!is.finite(xt)) return(c(Inf, k-1))
    } else {
      return(c(Inf, k-1))
    }
  }
  
  return(c(xt, k-1))
}

## ---------------------------------------------------------
## Your original Newton-Raphson (no damping)
## (slightly hardened: derivative ~0 -> Inf)
## ---------------------------------------------------------
root.NR = function(a3, a2, a1, a0, m, x0,
                   tol = 1e-12) {
  
  xt.save.2 = c(); yt.save.2 = c()
  xt = x0
  y1t = poly_f(xt, a3, a2, a1, a0, m)
  xt.save.2[1] = xt; yt.save.2[1] = y1t
  
  k = 1
  esp = 1
  
  while (esp > tol) {
    
    g.prime = poly_fp(xt, a3, a2, a1, a0, m)
    if (!is.finite(g.prime) || abs(g.prime) < .Machine$double.eps) {
      return(c(Inf, k-1))
    }
    
    xt.save.2[k+1] = xt - y1t/g.prime
    yt.save.2[k+1] = poly_f(xt.save.2[k+1], a3, a2, a1, a0, m)
    
    esp = abs(yt.save.2[k+1])
    xt = xt.save.2[k+1]
    y1t = yt.save.2[k+1]
    k = k + 1
    
    if (!is.finite(xt)) {
      return(c(Inf, k-1))
    }
  }
  
  if (esp > tol) xt <- Inf
  return(c(xt, k-1))
}

## ---------------------------------------------------------
## Your bisection method
## (kept as-is; note x0 is unused in your experiment)
## ---------------------------------------------------------
root.BS = function(a3, a2, a1, a0, m,
                   tol = 1e-12){
  
  a.bound = c(0, 2)
  xt.save.2 = c(); yt.save.2 = c()
  xt = sum(a.bound)/2
  y1t = poly_f(xt, a3, a2, a1, a0, m)
  
  xt.save.2[1] = xt; yt.save.2[1] = y1t
  k = 1; esp = 1
  
  while (esp > tol) {
    a.bound[2 - sum(yt.save.2[k] > 0)] = sum(a.bound)/2
    xt.save.2[k+1] = sum(a.bound)/2
    yt.save.2[k+1] = poly_f(xt.save.2[k+1], a3, a2, a1, a0, m)
    esp = abs(yt.save.2[k+1])
    k = k + 1
  }
  
  return(c(xt.save.2[k], k-1))
}

############################################################
#### Monte Carlo comparisons
############################################################
compare_roots <- function(a3, a2, a1, a0, m,
                          n.rep = 100000,
                          xlow = 0, xhigh = 3,
                          a.max = 2,
                          seed = NULL,
                          roots_ref = NULL,
                          tol_root = 1e-8,
                          warmup = 5000) {
  
  ## -------- helper: classify root (local) --------
  classify_root <- function(x, roots_ref, tol = 1e-8) {
    if (!is.finite(x) || length(x) != 1) return(NA_integer_)
    d <- abs(x - roots_ref)
    j <- which.min(d)
    if (is.finite(d[j]) && d[j] <= tol) j else 0L
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  ## shared initial values
  x0_vec <- runif(n.rep, xlow, xhigh)
  
  methods <- list(
    NR            = function(x0) root.NR(a3,a2,a1,a0,m=m,x0=x0),
    Bisection     = function(x0) root.BS(a3,a2,a1,a0,m=m),
    Halley        = function(x0) root.halley(a3,a2,a1,a0,m=m,x0=x0),
    Householder4  = function(x0) root.householder4(a3,a2,a1,a0,m=m,x0=x0),
    DampedNewton  = function(x0) root.newton.damped(a3,a2,a1,a0,m=m,x0=x0),
    Hirano        = function(x0) root.hirano(a3,a2,a1,a0,m=m,x0=x0),
    UHE           = function(x0) root.UHE(a3,a2,a1,a0,a.max=a.max,m=m,x0=x0),
    UHE3          = function(x0) root.UHE3(a3,a2,a1,a0,a.max=a.max,m=m,x0=x0)
  )
  
  ## reference roots (default from NR with 3 starts)
  if (is.null(roots_ref)) {
    NR1 <- root.NR(a3,a2,a1,a0,m=m,x0= 10)[1]
    NR2 <- root.NR(a3,a2,a1,a0,m=m,x0=  1)[1]
    NR3 <- root.NR(a3,a2,a1,a0,m=m,x0=-10)[1]
    roots_ref <- c(NR1, NR2, NR3)
  }
  roots_ref <- as.numeric(roots_ref)
  
  K <- length(roots_ref)
  root_names  <- paste0("root", seq_len(K))
  count_names <- c(root_names, "other")
  
  ## -------- warmup --------
  n.warm <- min(warmup, n.rep)
  for (j in seq_along(methods)) {
    for (i in seq_len(n.warm)) {
      invisible(methods[[j]](x0_vec[i]))
    }
  }
  
  out <- matrix(NA_real_, nrow = length(methods), ncol = (K+1))
  rownames(out) <- names(methods)
  colnames(out) <- paste0("freq_", count_names)
  
  time_mat <- matrix(NA_real_, nrow = length(methods), ncol = 3)
  rownames(time_mat) <- names(methods)
  colnames(time_mat) <- c("time_success_sec",
                          "avg_time_success_us",
                          "avg_iter_success")
  
  for (j in seq_along(methods)) {
    
    num.inva    <- 0L
    num.iter    <- 0
    counts      <- integer(K + 1)
    time_valid  <- 0
    n_success   <- 0L
    
    for (i in seq_len(n.rep)) {
      
      x0 <- x0_vec[i]
      
      t_call0 <- proc.time()[3]
      est     <- methods[[j]](x0)
      t_call1 <- proc.time()[3] 
      dt      <- t_call1 - t_call0
      
      xhat <- est[1]
      it   <- est[2]
      
      if (!is.finite(xhat) || is.nan(xhat)) {
        num.inva <- num.inva + 1L
        next
      }
      
      cls <- classify_root(xhat, roots_ref, tol = tol_root)
      if (is.na(cls)) {
        num.inva <- num.inva + 1L
        next
      }
      
      num.iter  <- num.iter + it
      n_success <- n_success + 1L
      time_valid <- time_valid + dt
      
      if (cls == 0L) {
        counts[K + 1L] <- counts[K + 1L] + 1L
      } else {
        counts[cls] <- counts[cls] + 1L
      }
    }
    
    denom_success <- n_success
    avg_iter <- if (denom_success > 0) num.iter / denom_success else NA_real_
    freq <- if (denom_success > 0) counts / denom_success else rep(NA_real_, K+1)
    
    out[j, ] <- freq
    
    time_success_sec <- time_valid
    avg_time_success_us <- if (denom_success > 0)
      (time_valid / denom_success) * 1e6 else NA_real_
    
    time_mat[j, ] <- c(time_success_sec,
                       avg_time_success_us,
                       avg_iter)
  }
  
  return(list(summary   = out,
              time      = time_mat,
              roots_ref = roots_ref,
              tol_root  = tol_root))
}

## ---------------------------------------------------------
## Run the experiment (your example)
## ---------------------------------------------------------
a3 = 1; a2 = -3; a1 = -1; a0 = 1; m = 3
n.rep = 100000

time.tra = compare_roots(a3,a2,a1,a0,m,
                         n.rep = n.rep,
                         xlow = 0, xhigh = 3,
                         a.max = 2,
                         seed = 1)


results1 = cbind(time.tra$summary[,c(1,2,3)], time.tra$time[,c(2,3)])

library(xtable)

tab <- as.data.frame(results1)
tab[, 1:3] <- lapply(tab[, 1:3], function(x) sprintf("%.2f\\%%", x * 100))

xt <- xtable(tab,
             caption = "Performance comparison of solvers",
             label   = "tab:perf")

print(xt,
      include.rownames = TRUE,
      comment = FALSE,
      sanitize.text.function = identity)

