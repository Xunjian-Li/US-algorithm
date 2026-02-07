## ---------------------------------------------------------
## Base score and derivatives: g, g', g'', g'''
## ---------------------------------------------------------
yule_loglik <- function(theta, x, max_trunc = 1000) {
  if (theta <= 0 || !is.finite(theta)) return(-Inf)
  n  <- length(x)
  C0 <- 0
  for (i in seq_len(n)) {
    x.len <- min(x[i], max_trunc)
    if (x.len <= 0) next
    idx   <- 0:(x.len - 1)
    denom <- theta + 1 + idx
    C0    <- C0 + sum(log(denom))
  }
  n * log(theta) - C0
}

yule_g <- function(theta, x, max_trunc = 1000) {
  n <- length(x)
  C1 <- 0
  for (i in seq_len(n)) {
    x.len <- min(x[i], max_trunc)
    if (x.len <= 0) next
    idx <- 0:(x.len - 1)
    denom <- theta + 1 + idx
    r1 <- 1/denom
    C1 <- C1 - sum(r1)
  }
  n/theta + C1
}

yule_gp <- function(theta, x, max_trunc = 1000) {
  n <- length(x)
  C2 <- 0
  for (i in seq_len(n)) {
    x.len <- min(x[i], max_trunc)
    if (x.len <= 0) next
    idx <- 0:(x.len - 1)
    denom <- theta + 1 + idx
    r1 <- 1/denom
    r2 <- r1^2
    C2 <- C2 + sum(r2)
  }
  
  -n/theta^2 + C2
}

yule_gpp <- function(theta, x, max_trunc = 1000) {
  n <- length(x)
  C3 <- 0
  
  for (i in seq_len(n)) {
    x.len <- min(x[i], max_trunc)
    if (x.len <= 0) next
    idx <- 0:(x.len - 1)
    denom <- theta + 1 + idx
    r1 <- 1/denom
    r3 <- r1^3
    C3 <- C3 + sum(r3)
  }
  2*n/theta^3 - 2*C3
}

yule_gppp <- function(theta, x, max_trunc = 1000) {
  n <- length(x)
  C4 <- 0
  for (i in seq_len(n)) {
    x.len <- min(x[i], max_trunc)
    if (x.len <= 0) next
    idx <- 0:(x.len - 1)
    denom <- theta + 1 + idx
    r1 <- 1/denom
    r4 <- r1^4
    C4 <- C4 + sum(r4)
  }
  
  -6*n/theta^4 + 6*C4
}

root.NR.yule <- function(x, theta0,
                         tol = 1e-8, maxit = 200,
                         max_trunc = 1000) {
  theta <- theta0
  g  <- yule_g(theta, x, max_trunc)
  k  <- 1L
  esp <- abs(g)
  
  while (esp > tol && k < maxit) {
    g1 <- yule_gp(theta, x, max_trunc)
    if (!is.finite(g1) || abs(g1) < .Machine$double.eps) return(c(Inf, k - 1L))
    
    ## Newton on g(θ) = 0:
    theta_new <- theta - g / g1
    
    if (!is.finite(theta_new)) return(c(Inf, k - 1L))
    if (theta_new <= 0) theta_new <- 1e-8
    
    g_new <- yule_g(theta_new, x, max_trunc)
    esp   <- abs(g_new)
    
    theta <- theta_new
    g     <- g_new
    k     <- k + 1L
  }
  
  c(theta, k - 1L)
}

root.BS.yule <- function(x,
                         L = 1e-6, U = 50,
                         tol = 1e-8, maxit = 500,
                         max_trunc = 1000) {
  gL <- yule_g(L, x, max_trunc)
  gU <- yule_g(U, x, max_trunc)
  if (!is.finite(gL) || !is.finite(gU) || gL*gU > 0) return(c(Inf, 0L))
  
  theta <- (L + U) / 2
  g  <- yule_g(theta, x, max_trunc)
  k  <- 1L
  esp <- abs(g)
  
  while (esp > tol && k < maxit) {
    if (gL*g <= 0) {
      U <- theta; gU <- g
    } else {
      L <- theta; gL <- g
    }
    theta <- (L + U)/2
    g     <- yule_g(theta, x, max_trunc)
    esp   <- abs(g)
    k     <- k + 1L
  }
  
  c(theta, k - 1L)
}

root.halley.yule <- function(x, theta0,
                             tol = 1e-8, maxit = 200,
                             max_trunc = 1000) {
  theta <- theta0
  g  <- yule_g(theta, x, max_trunc)
  k  <- 1L
  esp <- abs(g)
  
  while (esp > tol && k < maxit) {
    g1 <- yule_gp(theta, x, max_trunc)
    g2 <- yule_gpp(theta, x, max_trunc)
    
    denom <- 2*g1^2 - g*g2
    if (!is.finite(denom) || abs(denom) < .Machine$double.eps) return(c(Inf, k - 1L))
    
    theta_new <- theta - (2*g*g1)/denom
    if (!is.finite(theta_new)) return(c(Inf, k - 1L))
    if (theta_new <= 0) theta_new <- 1e-8
    
    g_new <- yule_g(theta_new, x, max_trunc)
    esp   <- abs(g_new)
    theta <- theta_new
    g     <- g_new
    k     <- k + 1L
  }
  
  c(theta, k - 1L)
}


## ---------------------------------------------------------
## Householder (order-3)
## x_{k+1} = x_k - (6 g g'^2 - 3 g^2 g'') / (6 g'^3 - 6 g g' g'' + g^2 g''')
## ---------------------------------------------------------
root.householder4.yule <- function(x, theta0,
                                   tol = 1e-8, maxit = 200,
                                   max_trunc = 1000) {
  theta <- theta0
  g  <- yule_g(theta, x, max_trunc)
  k  <- 1L
  esp <- abs(g)
  
  while (esp > tol && k < maxit) {
    g1 <- yule_gp(theta, x, max_trunc)
    g2 <- yule_gpp(theta, x, max_trunc)
    g3 <- yule_gppp(theta, x, max_trunc)
    
    denom <- 6*g1^3 - 6*g*g1*g2 + g^2*g3
    numer <- 6*g*g1^2 - 3*g^2*g2
    if (!is.finite(denom) || abs(denom) < .Machine$double.eps) return(c(Inf, k - 1L))
    
    theta_new <- theta - numer/denom
    if (!is.finite(theta_new)) return(c(Inf, k - 1L))
    if (theta_new <= 0) theta_new <- 1e-8
    
    g_new <- yule_g(theta_new, x, max_trunc)
    esp   <- abs(g_new)
    theta <- theta_new
    g     <- g_new
    k     <- k + 1L
  }
  
  c(theta, k - 1L)
}

root.newton.damped.yule <- function(x, theta0,
                                    tol        = 1e-8,
                                    maxit      = 200,
                                    lambda0    = 1.0,
                                    lambda_min = 1e-8,
                                    lambda_max = 1e8,
                                    max_trunc  = 1000) {
  
  theta  <- max(theta0, 1e-8)
  ell    <- yule_loglik(theta, x, max_trunc)
  lambda <- lambda0
  iter_used <- 0L
  
  for (k in 1:maxit) {
    g <- yule_g(theta, x, max_trunc)
    if (!is.finite(g)) return(c(Inf, k - 1L))
    if (abs(g) < tol)
      return(c(theta, k - 1L))
    
    H <- yule_gp(theta, x, max_trunc)
    if (!is.finite(H)) H <- 0
    
    accepted <- FALSE
    
    for (inner in 1:20) {
      H_eff <- H - lambda
      if (H_eff >= -1e-8) {
        H_eff <- -(abs(H) + lambda)
      }
      if (!is.finite(H_eff) || abs(H_eff) < .Machine$double.eps) {
        if (lambda >= lambda_max) {
          break
        } else {
          lambda <- min(lambda * 10, lambda_max)
          next
        }
      }
      
      step <- -g / H_eff
      if (!is.finite(step) || step == 0) {
        if (lambda >= lambda_max) {
          break
        } else {
          lambda <- min(lambda * 10, lambda_max)
          next
        }
      }
      
      theta_try <- theta + step
      theta_try <- max(theta_try, 1e-8)
      
      ell_try <- yule_loglik(theta_try, x, max_trunc)
      if (!is.finite(ell_try)) {
        if (lambda >= lambda_max) {
          break
        } else {
          lambda <- min(lambda * 10, lambda_max)
          next
        }
      }
      
      # armijo criterion
      pred_gain <- g * step + 0.5 * H * step^2
      if (!is.finite(pred_gain) || pred_gain <= 0) {
        pred_gain <- g * step
      }
      
      act_gain <- ell_try - ell
      rho <- if (pred_gain != 0) act_gain / pred_gain else -Inf
      
      if (!is.finite(rho)) {
        if (lambda >= lambda_max) {
          break
        } else {
          lambda <- min(lambda * 10, lambda_max)
          next
        }
      }
      
      if (rho > 0) {
        theta <- theta_try
        ell   <- ell_try
        
        if (rho > 0.75) {
          lambda <- max(lambda * 0.3, lambda_min)
        } else if (rho < 0.25) {
          lambda <- min(lambda * 2.0, lambda_max)
        }
        accepted  <- TRUE
        iter_used <- k
        break
      } else {
        if (lambda >= lambda_max) {
          break
        } else {
          lambda <- min(lambda * 2.0, lambda_max)
        }
      }
    }  # end inner loop
    if (!accepted) break
  }
  
  c(theta, iter_used)
}



root.hirano.yule <- function(x, theta0,
                             tol       = 1e-8,
                             maxit     = 200,
                             beta      = 0.5,
                             delta     = 0.5,
                             inner_max = 50,
                             max_trunc = 1000) {
  theta <- theta0
  
  for (iter in 1:maxit) {
    g0 <- yule_g(theta, x, max_trunc)
    if (abs(g0) < tol) return(c(theta, iter - 1L))
    
    g1 <- yule_gp(theta, x, max_trunc)
    g2 <- yule_gpp(theta, x, max_trunc)
    
    a0 <- g0
    a1 <- g1
    a2 <- g2/2
    
    mu_param <- 1
    success  <- FALSE
    
    for (j in 1:inner_max) {
      candidates <- complex(length = 2L)
      eps <- 1e-100
      
      ## k = 1
      if (is.finite(a1) && abs(a1) > eps) {
        candidates[1L] <- -mu_param * a0 / a1
      } else {
        candidates[1L] <- complex(real = Inf, imaginary = Inf)
      }
      
      ## k = 2
      if (is.finite(a2) && abs(a2) > eps) {
        u2 <- -mu_param * a0 / a2
        if (u2 >= 0) {
          candidates[2L] <- sqrt(u2)
        } else {
          candidates[2L] <- complex(real = Inf, imaginary = Inf)
        }
      } else {
        candidates[2L] <- complex(real = Inf, imaginary = Inf)
      }
      
      moduli <- Mod(candidates)
      if (all(is.infinite(moduli))) break
      
      m     <- which.min(moduli)
      zeta  <- candidates[m]
      
      poly_val <- a0 + a1*zeta + a2*zeta^2
      lhs <- Mod(poly_val)
      rhs <- (1 - (1 - beta) * mu_param) * abs(a0)
      
      if (is.finite(lhs) && is.finite(rhs) && lhs <= rhs) {
        theta_new <- Re(theta + zeta)
        if (theta_new <= 0) theta_new <- 1e-8
        theta <- theta_new
        success <- TRUE
        break
      } else {
        mu_param <- mu_param / (1 + delta)
      }
    }
    
    if (!success) return(c(theta, maxit))
  }
  
  c(theta, maxit)
}

## ---- US ----
root.US.yule <- function(x, theta0 = 1,
                         tol = 1e-8, maxit = 100, max_trunc = 1000) {
  theta <- theta0
  k     <- 0L
  esp   <- Inf
  n     <- length(x)
  
  while (esp > tol && k < maxit) {
    C1 <- 0
    for (i in 1:n) {
      x.len <- min(x[i], max_trunc)
      if (x.len <= 0) next
      C1 <- C1 - sum(1/(0:(x.len-1) + theta + 1))
    }
    
    gt      <- n/theta + C1
    b_theta <- -n/theta^2 + n/(theta+1)^2
    A1      <- gt - n/theta + n/(theta+1)
    disc    <- A1^2 - 4*A1*n
    if (disc < 0) disc <- 0
    
    theta_new <- -(A1 + sqrt(disc))/(2*A1)
    if (theta_new <= 0) theta_new <- 1e-8
    
    esp   <- abs(theta_new - theta)
    theta <- theta_new
    k     <- k + 1L
  }
  c(theta, k)
}

## ---- Fixed Point ----
root.FP.yule <- function(x, theta0 = 1,
                         tol = 1e-8, maxit = 200, max_trunc = 1000) {
  
  n <- length(x)
  theta <- max(theta0, 1e-8)
  k <- 0L
  esp <- Inf
  
  while (esp > tol && k < maxit) {
    
    C1 <- 0
    for (i in seq_len(n)) {
      x.len <- min(x[i], max_trunc)
      if (x.len <= 0) next
      C1 <- C1 - sum(1 / (0:(x.len - 1) + theta + 1))
    }
    
    A1 <- C1
    A2 <- A1 + n
    
    if (!is.finite(A1) || abs(A1) < 1e-14) {
      return(c(Inf, k))
    }
    
    disc <- A2^2 - 4*A1*n
    if (!is.finite(disc)) return(c(Inf, k))
    if (disc < 0) disc <- 0
    
    theta_new <- -(A2 + sqrt(disc)) / (2*A1)
    if (!is.finite(theta_new)) return(c(Inf, k))
    theta_new <- max(theta_new, 1e-8)
    
    esp <- abs(theta_new - theta)
    theta <- theta_new
    k <- k + 1L
  }
  
  if (k >= maxit && esp > tol) return(c(Inf, k))
  c(theta, k)
}

## ---- fast FP ----
root.fast_FP.yule <- function(x, theta0 = 1,
                              tol = 1e-8, maxit = 100, max_trunc = 1000) {
  theta <- theta0
  k     <- 0L
  esp   <- Inf
  n     <- length(x)
  
  while (esp > tol && k < maxit) {
    C1 <- 0
    C2 <- 0
    for (i in 1:n) {
      x.len <- min(x[i], max_trunc)
      if (x.len <= 0) next
      idx <- 0:(x.len - 1)
      denom <- theta + 1 + idx
      r1 <- 1/denom
      C1 <- C1 - sum(r1)
      C2 <- C2 + sum(r1^2)
    }
    
    A1 <- C1
    A2 <- A1 + n
    
    gt     <- n/theta + C1                    # g(θ)
    gprime <- -n/theta^2 + C2                 # g'(θ)
    b_theta <- -n/theta^2
    
    if (!is.finite(A1) || abs(A1) < 1e-14) {
      return(c(Inf, k))
    }
    
    disc <- A2^2 - 4*A1*n
    if (!is.finite(disc)) return(c(Inf, k))
    if (disc < 0) disc <- 0
    
    theta_quad <- -(A2 + sqrt(disc)) / (2*A1)
    if (!is.finite(theta_quad)) return(c(Inf, k))
    theta_quad <- max(theta_quad, 1e-8)
    
    s.ratio <- min(b_theta/gprime*(gprime < 0) + (gprime > 0), 2)
    theta_new <- theta + s.ratio*(theta_quad - theta)
    if (theta_new <= 0) theta_new <- 1e-8
    
    esp   <- abs(theta_new - theta)
    theta <- theta_new
    k     <- k + 1L
  }
  c(theta, k)
}

## ---- fast US ----
root.fast_US.yule <- function(x, theta0 = 1,
                              tol = 1e-8, maxit = 100, max_trunc = 1000) {
  theta <- theta0
  k     <- 0L
  esp   <- Inf
  n     <- length(x)
  
  while (esp > tol && k < maxit) {
    C1 <- 0
    C2 <- 0
    for (i in 1:n) {
      x.len <- min(x[i], max_trunc)
      if (x.len <= 0) next
      idx <- 0:(x.len - 1)
      denom <- theta + 1 + idx
      r1 <- 1/denom
      C1 <- C1 - sum(r1)
      C2 <- C2 + sum(r1^2)
    }
    
    gt     <- n/theta + C1                    # g(θ)
    gprime <- -n/theta^2 + C2                 # g'(θ)
    b_theta <- -n/theta^2 + n/(theta+1)^2
    A1 <- gt - n/theta + n/(theta+1)
    disc <- A1^2 - 4*A1*n
    if (disc < 0) disc <- 0
    theta_quad <- -(A1 + sqrt(disc))/(2*A1)
    
    s.ratio <- min(b_theta/gprime*(gprime < 0) + (gprime > 0), 2)
    theta_new <- theta + s.ratio*(theta_quad - theta)
    if (theta_new <= 0) theta_new <- 1e-8
    
    esp   <- abs(theta_new - theta)
    theta <- theta_new
    k     <- k + 1L
  }
  c(theta, k)
}


library(VGAM)

compare_yule_mle <- function(theta_true,
                             n,
                             n.rep      = 1000,
                             theta0_min = 0,
                             theta0_max = 5,
                             seed       = NULL,
                             max_trunc  = 1000) {
  if (!is.null(seed)) set.seed(seed)
  
  methods <- list(
    Newton        = function(x, theta0) root.NR.yule(x, theta0, max_trunc = max_trunc),
    Bisection     = function(x, theta0) root.BS.yule(x, max_trunc = max_trunc),
    Halley        = function(x, theta0) root.halley.yule(x, theta0, max_trunc = max_trunc),
    Householder   = function(x, theta0) root.householder4.yule(x, theta0, max_trunc = max_trunc),
    DampedNewton  = function(x, theta0) root.newton.damped.yule(x, theta0, max_trunc = max_trunc),
    Hirano        = function(x, theta0) root.hirano.yule(x, theta0, max_trunc = max_trunc),
    FP            = function(x, theta0) root.FP.yule(x, theta0, max_trunc = max_trunc),
    Fast_FP       = function(x, theta0) root.fast_FP.yule(x, theta0, max_trunc = max_trunc),
    US            = function(x, theta0) root.US.yule(x, theta0, max_trunc = max_trunc),
    Fast_US       = function(x, theta0) root.fast_US.yule(x, theta0, max_trunc = max_trunc)
  )
  
  out <- matrix(NA_real_, nrow = length(methods), ncol = 5)
  rownames(out) <- names(methods)
  colnames(out) <- c("time_success_sec", "avg_time_success_us",
                     "num_invalid", "num_success", "avg_iter_success")
  
  for (j in seq_along(methods)) {
    num.inva <- 0L
    num.iter <- 0.0
    num.succ <- 0L
    time.ok  <- 0.0
    
    for (i in seq_len(n.rep)) {
      x       <- ryules(n, theta_true)
      theta0  <- runif(1, theta0_min, theta0_max)
      
      t0 <- proc.time()[3]
      est <- methods[[j]](x, theta0)
      t1 <- proc.time()[3]
      
      theta_hat <- est[1]
      it        <- est[2]
      
      ok <- is.finite(theta_hat) && !is.nan(theta_hat) && theta_hat > 0
      if (!ok) {
        num.inva <- num.inva + 1L
        next
      }
      
      num.succ <- num.succ + 1L
      num.iter <- num.iter + it
      time.ok  <- time.ok + (t1 - t0)
    }
    
    avg_iter <- if (num.succ > 0) num.iter/num.succ else NA_real_
    avg_time_ms <- if (num.succ > 0) (time.ok/num.succ) else NA_real_
    
    out[j, ] <- c(time.ok*1e3, avg_time_ms*1e3, num.inva, num.succ, avg_iter)
  }
  
  out
}


theta_true <- 0.5
n          <- 400
n.rep      <- 100

res_yule1 <- compare_yule_mle(theta_true, n,
                             n.rep = n.rep,
                             theta0_min = 0,
                             theta0_max = 5,
                             seed = 123)

theta_true <- 5
n          <- 400
n.rep      <- 100

res_yule2 <- compare_yule_mle(theta_true, n,
                             n.rep = n.rep,
                             theta0_min = 0,
                             theta0_max = 5,
                             seed = 123)

results1 = cbind(res_yule1, res_yule2)

library(xtable)

tab <- as.data.frame(results1[,c(2,4,5,7,9,10)])
tab[, c(2,5)] <- lapply(tab[, c(2,5)], function(x) sprintf("%.2f\\%%", x ))
tab[, c(1,3)] <- lapply(tab[, c(1,3)], function(x) sprintf("%.2f", x ))

xt <- xtable(tab,
             caption = "Performance comparison of solvers",
             label   = "tab:perf")

print(xt,
      include.rownames = TRUE,
      comment = FALSE,
      sanitize.text.function = identity)



print(
  xtable(results1,
         digits = 4,
         caption = "Performance comparison of solvers",
         label = "tab:perf"),
  include.rownames = TRUE,
  comment = FALSE
)

