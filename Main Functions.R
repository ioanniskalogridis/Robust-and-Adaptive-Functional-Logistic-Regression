###################################################################################################################################
############## Robust and Adaptive Functional Linear Regression (Kalogridis (2023)) ###############################################
require(fda)
require(MASS)
require(maxLik)

dpd.f <- function(x, y, m = 2, nbasis = NULL,  norder = 4, toler = 1e-08, maxiter = 1000, nsteps = 20, tuning = NULL,
                  alpha.cand = c(seq(1e-04, 2, len = 20))){
  
  # Main function
  # x is a matrix contained the values of the discretized predictors
  # y is the response variable
  # nbasis is the number of B-spline basis functions, by default nbasis = [min(n/4, 30)]
  # norder is the order of the spline basis, by default a cubic spline basis
  # m is the order of the penalty, by default the penalty is placed on the integrated squared second derivative
  # toler is the tolerance level of the Fisher-scoring algorithm
  # maxiter is the naximum number of allowed iterations of the algorithm
  # tuning refers to the way of selecting the tuning parameter. The default "tuning = NULL" corresponds to selection in
  # a robust and adaptive way, as outlined in the paper. For a fixed tuning parameter, e.g., equal to $1$, set "tuning = 1",
  # alpha.cand is a vector of candidate values for the tuning parameter (kappa in the notation of the paper)
  # a pilot estimator corresponding to kappa = 2.5 is used, as discussed in the paper
  # nsteps denotes the number of allowed iterations for the selection of the tuning parameter
  
  x <- as.matrix(x)
  n <- length(y)
  y <- as.vector(y)
  if(is.null(nbasis)){
    nbasis = floor(min(n/4, 30))
  } else{
    nbasis = nbasis
  }
  b.sp <- create.bspline.basis(c(0, 1), nbasis = nbasis, norder = norder)
  b.sp.e <- eval.basis(seq(1/dim(x)[2], 1-1/dim(x)[2], len = dim(x)[2]), b.sp)
  x.p.ni <- x%*%b.sp.e/dim(x)[2]
  x.p <- cbind(rep(1, n), x.p.ni)
  
  p.m <- bsplinepen(b.sp, Lfdobj = m )
  p.m <- rbind( rep(0, nbasis+1 ), cbind(rep(0, nbasis), p.m) )
  
  inv.logit <- function(x) 1/(1+exp(-x))
  beta.in <- rep(0, dim(x.p)[2])
  
  obj.f <- function(beta, pen, tun){
    lambda = pen
    alpha = tun
    t1 <- (1-inv.logit(x.p%*%beta))^{1+alpha} + inv.logit(x.p%*%beta)^{1+alpha}
    t2 <- (1+1/alpha)*dbinom(y, size = 1, prob = inv.logit(x.p%*%beta))^{alpha}
    contr <- t1 - t2
    obj <- as.numeric(mean(contr) + lambda*t(beta)%*%(p.m%*%beta))
    return(-obj)
  }
  
  gradient <- function(beta, pen, tun){
    lambda = pen
    alpha = tun
    gr1 <- (0-inv.logit(x.p%*%beta))*(1-inv.logit(x.p%*%beta))^{1+alpha}+ (1-inv.logit(x.p%*%beta))*(inv.logit(x.p%*%beta))^{1+alpha}
    gr2 <-  (y-inv.logit(x.p%*%beta))*dbinom(y, size = 1 , prob = inv.logit( x.p%*%beta ))^{alpha}
    gr.t <- gr1 - gr2
    gr <- (1+alpha)*diag(c(gr.t))%*%x.p
    gradient <- as.vector(colSums(gr)/n + 2*lambda*(p.m%*%beta))
    return(-gradient)
  }
  
  hessian <-  function(beta, pen, tun){
    lambda = pen
    alpha = tun
    wt <- inv.logit(x.p%*%beta)^2*(1-inv.logit(x.p%*%beta))^{1+alpha} + (1-inv.logit(x.p%*%beta ) )^2*(inv.logit(x.p%*%beta))^{1+alpha}
    wt <- (1+alpha)*wt
    hs <- scale(t(x.p), center = FALSE, scale = c(1/wt))%*%x.p/n + 2*lambda*p.m
    return(-hs)
  }
  
  newraph <- function(lambda, beta.in, maxit = maxiter, tol = toler, alpha){
    mnr <- maxNR(fn = obj.f, grad = gradient, hess = hessian, start = beta.in, pen =  lambda, tun =  alpha,
                 control = list(tol = tol, reltol = tol, gradtol = tol, iterlim = maxit))
    beta = mnr$estimate
    gr1 <- (0-inv.logit(x.p%*%beta))*(1-inv.logit(x.p%*%beta))^{1+alpha}+ (1-inv.logit(x.p%*%beta))*(inv.logit(x.p%*%beta))^{1+alpha}
    gr2 <-  (y-inv.logit(x.p%*%beta))*dbinom(y, size = 1 , prob = inv.logit( x.p%*%beta ))^{alpha}
    gr.t <- gr1 - gr2
    gr <- c((1+alpha)*gr.t)
    hess.f <- -mnr$hessian
    Q.m <- hess.f - 2*lambda*p.m
    hat.tr <- sum(diag(  solve(hess.f, Q.m, tol = 1e-25) ))/n
    return(list(beta = beta, hat.tr = hat.tr, gr = gr, hess.f = hess.f))
  }
  
  obj.f.np <- function(beta, alpha){
    t1 <- (1-inv.logit(x.p%*%beta))^{1+alpha} + inv.logit(x.p%*%beta)^{1+alpha}
    t2 <- (1+1/alpha)*dbinom(y, size = 1, prob = inv.logit(x.p%*%beta))^{alpha}
    contr <- t1 - t2
    obj <- as.numeric(mean(contr))
    return(obj)
  }
  
  Pen.cr <- function(lambda, beta.in, alpha){
    nr <- newraph(lambda, beta.in = beta.in, maxit = maxiter, alpha = alpha)
    beta <- nr$beta
    hat.tr <- nr$hat.tr
    IC <- 2*(obj.f.np(beta, alpha)) + 2*hat.tr
    return(IC)
  }
  
  lambda.cand <- c(1e-12, 6e-12, 1e-11, 6e-11, 1e-10, 6e-10, 1e-09, 5e-09, 9e-09, 1e-08, 5e-08, 9e-08, 1e-07, 5e-07, 9e-07, 1e-06, 5e-06, 9e-06, 1e-05, 5e-05, 9e-05, 1e-04, 5e-04, 9e-04,
                   1e-03, 5e-03, 9e-03, 1e-02, 5e-02, 9e-02, 1e-01, 5e-01, 9e-01, 5)
  lambda.e.in <- rep(0, length(lambda.cand))
  for(k in 1:length(lambda.e.in)){
    lambda.e.in[k] <-  try(Pen.cr(lambda.cand[k], beta.in = beta.in, alpha = 2),silent = TRUE)
  }
  lambda.opt.in <- lambda.cand[which.min(lambda.e.in)]
  beta.opt.in <- newraph(lambda.opt.in, beta.in = beta.in, maxit = maxiter, alpha = 2)$beta
  B.m <-  bsplinepen(b.sp, Lfdobj = 0)
  
  if(is.null(tuning)){
    comp.alpha <- function(alpha){
      lambda.e <- rep(0, length(lambda.cand))
      for(k in 1:length(lambda.e)){
        lambda.e[k] <-  try(Pen.cr(lambda.cand[k], beta.in = beta.in, alpha = alpha),silent = TRUE)
      }
      lambda.opt <- suppressWarnings(lambda.cand[which.min(lambda.e)])
      opt. <- newraph(lambda.opt, beta.in = beta.opt.in, maxit = maxiter, alpha = alpha)
      beta.opt  <- opt.$beta
      gr <- opt.$gr
      hessian.m <- -opt.$hess.f
      K.m <- scale(t(x.p), center = FALSE, scale = c(1/gr^2))%*%x.p/(n^2) 
      msq1 <- sum( diag( B.m%*%solve(hessian.m, K.m%*%solve(hessian.m, tol = 1e-25), tol = 1e-25)[2:(nbasis+1), 2:(nbasis+1)]   ) )
      return(list(msq1 = msq1,  beta.opt = beta.opt, lambda.opt = lambda.opt))
    }
    
    compts <-  vapply(alpha.cand, FUN = comp.alpha, FUN.VALUE = rep(list(1), 3) )
    ic <- 0
    istop <- 0
    alpha.in <- 2
    while(ic <= nsteps & istop == 0){
      ic = ic + 1 
      msqs <- rep(NA, length(alpha.cand))
      for(j in 1:length(alpha.cand)){
        msqs[j] <- compts[[1, j]] + t(compts[[2, j]][-1]-compts[[2, which(alpha.cand == alpha.in)]][-1])%*%B.m%*%(compts[[2, j]][-1]-compts[[2, which(alpha.cand == alpha.in)]][-1])
      }
      alpha.up <- alpha.cand[which.min(msqs)]
      check <- abs(alpha.up-alpha.in)
      if(check< 1e-06){istop = 1}
      alpha.in <- alpha.up
    }
    alpha.opt <- alpha.up
    beta.opt.f <- compts[[2, which(alpha.cand == alpha.opt)]]
    lambda.opt.f <- compts[[3, which(alpha.cand == alpha.up)]]
  } else {
    lambda.e <- rep(0, length(lambda.cand))
    for(k in 1:length(lambda.e)){
      lambda.e[k] <-  try(Pen.cr(lambda.cand[k], beta.in = beta.in, alpha = tuning),silent = TRUE)
    }
    lambda.opt.f <- suppressWarnings(lambda.cand[which.min(lambda.e)])
    opt. <- newraph(lambda.opt.f, beta.in = beta.opt.in, maxit = maxiter, alpha = tuning)
    beta.opt.f  <- opt.$beta
    gr <- opt.$gr
  }
  
  est <-  b.sp.e%*%beta.opt.f[-1]
  fitted <- inv.logit(beta.opt.f[1]+x%*%est/dim(x)[2])
  resids <- as.vector(y - inv.logit(fitted))
  p.resids <- as.vector(resids/sqrt(inv.logit(resids)*(1-inv.logit(resids))))
  
  ibeta <- function(x,a,b){ pbeta(x,a,b)*beta(a,b)}
  ansch.r <- function(y, mu){
    ansch.r <- (ibeta(y, 2/3, 2/3)-ibeta(mu, 2/3, 2/3))/(mu*(1-mu))^{1/6}
    return(ansch.r)
  }
  a.resids <- ansch.r(y, fitted)
  
  return(list(est = est, a.resids = a.resids, alpha = ifelse(is.null(tuning),alpha.opt, tuning), ic = ifelse(is.null(tuning),ic, 1), 
              lambda = lambda.opt.f))
}
###########################################################################################################################
###########################################################################################################################

# dpd.ffa <- function(x, y, m = 2, toler = 1e-08, maxiter = 1000, nbasis = NULL, norder = 4,  alpha = 1e-04){
#   
#   # Density power divergence with fixed tuning 
#   # alpha is the tuning parameter, by default a very small alpha corresponding to the penalized likelihood estimator
#   
#   x <- as.matrix(x)
#   n <- length(y)
#   y <- as.vector(y)
#   if(is.null(nbasis)){
#     nbasis = floor(min(n/4, 30))
#   } else{
#     nbasis = nbasis
#   }
#   
#   b.sp <- create.bspline.basis(c(0, 1), nbasis = nbasis, norder = norder)
#   b.sp.e <- eval.basis(seq(1/dim(x)[2], 1-1/dim(x)[2], len = dim(x)[2]), b.sp)
#   x.p.ni <- x%*%b.sp.e/dim(x)[2]
#   x.p <- cbind(rep(1, n), x.p.ni)
#   
#   p.m <- bsplinepen(b.sp, Lfdobj = m )
#   p.m <- rbind( rep(0, nbasis+1 ), cbind(rep(0, nbasis), p.m) )
#   
#   inv.logit <- function(x) 1/(1+exp(-x))
#   beta.il <- rep(0, dim(x.p)[2])
#   
#   obj.f <- function(beta, pen, tun){
#     lambda = pen
#     alpha = tun
#     t1 <- (1-inv.logit(x.p%*%beta))^{1+alpha} + inv.logit(x.p%*%beta)^{1+alpha}
#     t2 <- (1+1/alpha)*dbinom(y, size = 1, prob = inv.logit(x.p%*%beta))^{alpha}
#     contr <- t1 - t2
#     obj <- as.numeric(mean(contr) + lambda*t(beta)%*%(p.m%*%beta))
#     return(-obj)
#   }
#   
#   gradient <- function(beta, pen, tun){
#     lambda = pen
#     alpha = tun
#     gr1 <- (0-inv.logit(x.p%*%beta))*(1-inv.logit(x.p%*%beta))^{1+alpha}+ (1-inv.logit(x.p%*%beta))*(inv.logit(x.p%*%beta))^{1+alpha}
#     gr2 <-  (y-inv.logit(x.p%*%beta))*dbinom(y, size = 1 , prob = inv.logit( x.p%*%beta ))^{alpha}
#     gr.t <- gr1 - gr2
#     gr <- (1+alpha)*diag(c(gr.t))%*%x.p
#     gradient <- as.vector(colSums(gr)/n + 2*lambda*(p.m%*%beta))
#     return(-gradient)
#   }
#   
#   hessian <-  function(beta, pen, tun){
#     lambda = pen
#     alpha = tun
#     wt <- inv.logit(x.p%*%beta)^2*(1-inv.logit(x.p%*%beta))^{1+alpha} + (1-inv.logit(x.p%*%beta ) )^2*(inv.logit(x.p%*%beta))^{1+alpha}
#     wt <- (1+alpha)*wt
#     hs <- scale(t(x.p), center = FALSE, scale = c(1/wt))%*%x.p/n + 2*lambda*p.m
#     return(-hs)
#   }
#   
#   newraph <- function(lambda, beta.in, maxit = maxiter, tol = toler, alpha){
#     mnr <- maxNR(fn = obj.f, grad = gradient, hess = hessian, start = beta.in, pen =  lambda, tun =  alpha,
#                  control = list(tol = tol, reltol = tol, gradtol = tol, iterlim = maxit))
#     beta = mnr$estimate
#     gr1 <- (0-inv.logit(x.p%*%beta))*(1-inv.logit(x.p%*%beta))^{1+alpha}+ (1-inv.logit(x.p%*%beta))*(inv.logit(x.p%*%beta))^{1+alpha}
#     gr2 <-  (y-inv.logit(x.p%*%beta))*dbinom(y, size = 1 , prob = inv.logit( x.p%*%beta ))^{alpha}
#     gr.t <- gr1 - gr2
#     gr <- c((1+alpha)*gr.t)
#     hess.f <- -mnr$hessian
#     Q.m <- hess.f - 2*lambda*p.m
#     hat.tr <- sum(diag(  solve(hess.f, Q.m, tol = 1e-22) ))/n
#     return(list(beta = beta, hat.tr = hat.tr, gr = gr, hess.f = hess.f))
#   }
#   
#   obj.f.np <- function(beta, alpha){
#     t1 <- (1-inv.logit(x.p%*%beta))^{1+alpha} + inv.logit(x.p%*%beta)^{1+alpha}
#     t2 <- (1+1/alpha)*dbinom(y, size = 1, prob = inv.logit(x.p%*%beta))^{alpha}
#     contr <- t1 - t2
#     obj <- as.numeric(mean(contr))
#     return(obj)
#   }
#   
#   Pen.cr <- function(lambda, beta.in, alpha){
#     nr <- newraph(lambda, beta.in = beta.in, maxit = maxiter, alpha = alpha)
#     beta <- nr$beta
#     hat.tr <- nr$hat.tr
#     IC <- 2*(obj.f.np(beta, alpha)) + 2*hat.tr
#     return(IC)
#   }
#   
#   lambda.cand <- c(1e-12, 6e-12, 1e-11, 6e-11, 1e-10, 6e-10, 1e-09, 5e-09, 9e-09,  1e-08, 5e-08, 9e-08, 1e-07, 5e-07, 9e-07, 1e-06, 5e-06, 9e-06, 1e-05, 5e-05, 9e-05, 1e-04, 5e-04, 9e-04, 
#                    1e-03, 5e-03, 9e-03, 1e-02, 5e-02, 9e-02, 1e-01, 5e-01, 9e-01, 5)
#   lambda.e <-  vapply(lambda.cand, FUN = Pen.cr, beta.in = beta.il, alpha = alpha, FUN.VALUE = numeric(1))
#   lambda.opt.f <- lambda.cand[which.min(lambda.e)]
#   
#   beta.opt.f <- newraph(lambda.opt.f, beta.in = beta.il, maxit = maxiter, alpha = alpha)$beta
#   
#   est <-  b.sp.e%*%beta.opt.f[-1]
#   fitted <- inv.logit(beta.opt.f[1]+x%*%est/dim(x)[2])
#   resids <- as.vector(y - inv.logit(fitted))
#   p.resids <- as.vector(resids/sqrt(inv.logit(resids)*(1-inv.logit(resids))))
#   
#   ibeta <- function(x,a,b){ pbeta(x,a,b)*beta(a,b)}
#   ansch.r <- function(y, mu){
#     ansch.r <- (ibeta(y, 2/3, 2/3)-ibeta(mu, 2/3, 2/3))/(mu*(1-mu))^{1/6}
#     return(ansch.r)
#   }
#   a.resids <- ansch.r(y, fitted)
#   
#   return(list(est = est, a.resids = a.resids, lambda = lambda.opt.f))
# }