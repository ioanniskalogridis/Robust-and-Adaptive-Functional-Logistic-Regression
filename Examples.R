##################################################### Clean data #########################################################
##########################################################################################################################

nrep <- 1000
n <- 400
p <- 200
grid <- seq(1/p, 1-1/p, len = p)

m1 <- m2 <- m3 <- m4 <- matrix(NA, nrow = p, ncol = nrep)
msq1 = msq2 = msq03 = msq4 =  rep(NA, nrep)

for(k in 1:nrep){
  print(k)
  # f1 <- -sin(5*grid/1.2)/0.5-1
  # f1 <- 3*sin(3.4*grid^2)
  f1 <- 3*(grid-0.3)^2+1
  x <- matrix(0, n, p)
  for(i in 1:n){
    x[i, ] <- sqrt(2)*(1*pi-pi/2)^{-1}*rnorm(1, 0, 1)*sapply(grid, FUN= function(x) sin((1-1/2)*pi*x)  )
    for(j in 2:50){
      x[i, ] <- x[i, ] + (j*pi-pi/2)^{-1}*rnorm(1, 0, 1)*sqrt(2)*sapply(grid, FUN= function(x) sin((j-1/2)*pi*x)  )
    }
  }
  
  y0 = x%*%f1/p
  y <- rbinom(n, size = 1, prob = inv.logit(y0))
  
  fit1 <- dpd.f(x = x, y = y, norder = 4, m = 2)
  fit2 <- dpd.ffa(x, y)
  fit3 <- dpd.ffa(x = x, y = y, norder = 4, m = 2, alpha = 2.5)
  fit4 <- dpd.ffa(x = x, y = y, norder = 4, m = 2, alpha = 5)

  msq1[k] <- mean((f1-fit1$est)^2)
  msq2[k] <- mean((f1-fit2$est)^2)
  msq3[k] <- mean((f1-fit3$est)^2)
  msq4[k] <- mean((f1-fit4$est)^2)
  
  m1[, k] <- fit1$est
  m2[, k] <- fit2$est
  m3[, k] <- fit3$est
  m4[, k] <- fit4$est
  
}

matplot(grid, m0.o, type = "l", col = "gray", lwd = 3, cex = 2.5, cex.axis = 3, cex.lab = 2.5, xlab = "", ylab = ""); lines(grid, f1, lwd = 3, col = "black") ;grid()
matplot(grid, ma.o, type = "l", col = "gray", lwd = 3,  cex = 2.5, cex.axis = 3, cex.lab = 2.5, xlab = "", ylab = ""); lines(grid, f1, lwd = 3, col = "black") ; grid()

mean(msq, na.rm = TRUE)
median(msq, na.rm = TRUE) ; median(msq0, na.rm = TRUE) ; median(msq05, na.rm = TRUE) ; median(msq1, na.rm = TRUE)
median(alpha.sel, na.rm = TRUE)

boot.median <- function(x){
  b <- 10000
  med.v <- rep(0, length(b))
  for(i in 1:b){
    r.samp <- sample(x, size = 1000, replace = TRUE)
    med.v[i] <- median(r.samp)
  }
  sd.med <- sd(med.v)
  print(sd.med)
}

boot.median(msq)
boot.median(msq0)
boot.median(msq05)
boot.median(msq1)

plot(grid, f1, type = "l", lwd = 3)
lines(grid, fit$est, type = "l", lwd = 3, col = "blue")
lines(grid, fit1$est, lwd = 3)
lines(grid, fit0$est, lwd = 3, col = "red")

###################################################### Contamination #######################################################
############################################################################################################################

nrep <- 1000
n <- 400
p <- 200
eps <- 0.01
grid <- seq(1/p, 1-1/p, len = p)

msq = msq0 = msq05 = msq1 =  rep(NA, nrep)
alpha.sel <- rep(NA, nrep)
ic <- rep(NA, nrep)

ma <- m0 <- m05 <- m1 <- matrix(NA, nrow = p, ncol = nrep)

for(k in 1:nrep){
  print(k)
  # f1 <- 2*sin(40*grid*pi/3)
  # f1 <- -sin(5*grid/1.2)/0.5-1
  # f1 <- 3*exp(-grid^2)
  # f1 <- 1/(1+5*grid)
  # f1 <- 2*cos(4*grid*pi/3)
  f1 <- 2*cos(grid*pi/4)
  # f1 <- 3*sin(3.4*grid^2)
  # f1 <- 2*cos(4*pi*grid)
  # f1 <- 5*grid^2
  # x <- sqrt(2)*matrix(rnorm(n*p), nrow = n, ncol = p)
  # X0 <- sqrt(2)*matrix(rt(n*p, df = 3), nrow = n, ncol = p)
  x <- matrix(0, n, p)
  for(i in 1:n){
    x[i, ] <- sqrt(2)*(1*pi-pi/2)^{-1}*rnorm(1, 0, 1)*sapply(grid, FUN= function(x) sin((1-1/2)*pi*x)  )
    for(j in 2:50){
      x[i, ] <- x[i, ] + (j*pi-pi/2)^{-1}*rnorm(1, 0, 1)*sqrt(2)*sapply(grid, FUN= function(x) sin((j-1/2)*pi*x)  )
    }
  }
  samp <- sample(1:n, size = floor(n*eps))
  # for(i in samp){
  #   x[i, ] <- sqrt(2)*(1*pi-pi/2)^{-1}*rt(1, df = 1)*sapply(grid, FUN= function(x) sin((1-1/2)*pi*x)  )
  #   for(j in 2:10){
  #     x[i, ] <- x[i, ] + (j*pi-pi/2)^{-1}*rt(1, df = 1)*sqrt(2)*sapply(grid, FUN= function(x) sin((j-1/2)*pi*x)  ) 
  #   }
  # }
  x[samp, ] <- 5*x[samp, ]

  y0 = (x%*%f1)/p
  
  # sum(abs(y0)>3)
  # plot(inv.logit(y0))
  
  # plot(y0,inv.logit(y0))
  inv.logit <- function(x) 1/(1+exp(-x))
  y <- rbinom(n, size = 1, prob = inv.logit(y0))
  y[samp]<- 1- y[samp]
  
  fit <- dpd.f(x = x, y = y, norder = 4, m = 2)
  # fit0 <- dpd.ffa(x = x, y = y, norder = 4, m = 2, alpha = 1e-20)
  fit05 <- dpd.ffa(x = x, y = y, norder = 4, m = 2, alpha = 2.5)
  fit1 <- dpd.ffa(x = x, y = y, norder = 4, m = 2, alpha = 5)
  fit0 <- p.ml(x, y, norder = 4, nbasis = 30, m = 2)
  msq[k] <- mean((f1-fit$est)^2)
  alpha.sel[k] <- fit$alpha
  ic[k] <- fit$ic
  
  ma[, k] <- fit$est
  m0[, k] <- fit0$est
  m1[, k] <- fit1$est
  m05[, k] <- fit05$est
  
  msq0[k] <- mean((f1-fit0$est)^2)
  msq05[k] <- mean((f1-fit05$est)^2)
  msq1[k] <- mean((f1-fit1$est)^2)
}

matplot(grid, m0, type = "l", col = "gray", lwd = 3)
matplot(grid, ma, type = "l", col = "gray", lwd = 3); lines(grid, f1, lwd = 3, col = "black")

hist(alpha.sel)
# hist(ic)
mean(msq, na.rm = TRUE)
median(msq, na.rm = TRUE) ; median(msq0, na.rm = TRUE) ; median(msq05, na.rm = TRUE) ; median(msq1, na.rm = TRUE)
median(alpha.sel, na.rm = TRUE)

boot.median <- function(x){
  b <- 10000
  med.v <- rep(0, length(b))
  for(i in 1:b){
    r.samp <- sample(x, size = 1000, replace = TRUE)
    med.v[i] <- median(r.samp)
  }
  sd.med <- sd(med.v)
  print(sd.med)
}

boot.median(msq)
boot.median(msq0)
boot.median(msq05)
boot.median(msq1)



plot(grid, f1, type = "l", lwd = 3)
lines(grid, fit$est, type = "l", lwd = 3, col = "blue")
lines(grid, fit1$est, lwd = 3)
lines(grid, fit0$est, lwd = 3, col = "red")

#################################################################################################################################################################
#################################################################################################################################################################

nrep <- 30
n <- 300
p <- 250
grid <- seq(1/p, 1-1/p, len = p)

msq = msq0 = msq05 = msq1 =  rep(NA, nrep)
alpha.sel <- rep(NA, nrep)
ic <- rep(NA, nrep)

for(k in 1:nrep){
  print(k)
  # f1 <- -sin(5*grid/1.2)/0.8-1
  # f1 <- 1.8*sin(3.4*grid^2)
  # f1 <- cos(2*pi*grid)
  # f1 <- grid^2
  # x <- sqrt(2)*matrix(rnorm(n*p), nrow = n, ncol = p)
  # X0 <- sqrt(2)*matrix(rt(n*p, df = 3), nrow = n, ncol = p)
  f1 <- sin(10*grid*pi/3)
  f2 <- 2*sin(40*grid*grid/3)
  # x1 <- matrix(0, floor(0.9*n), p)
  pred.basis <- create.bspline.basis(c(0,1), nbasis = 13)
  pred.basis <- eval.basis(grid,pred.basis)
  C.m <- matrix(rnorm(n*13), nrow = n, ncol = 13)%*%matrix(runif(13*13), nrow = 13, ncol = 13)
  x <- t(pred.basis%*%t(C.m))  
  
  # for(i in 1:floor(0.9*n)){
  # x1[i, ] <- sqrt(2)*(1*pi-pi/2)^{-1}*rnorm(1, 0, 1)*sapply(grid, FUN= function(x) sin((1-1/2)*pi*x)  )
  # x[i, ] <- sqrt(2)*(1*pi-pi/2)^{-1}*rt(1, df = 3)*sapply(grid, FUN= function(x) sin((1-1/2)*pi*x)  )
  # for(j in 2:50){
  # x1[i, ] <- x1[i, ] + (j*pi-pi/2)^{-1}*rnorm(1, 0, 1)*sqrt(2)*sapply(grid, FUN= function(x) sin((j-1/2)*pi*x)  )
  # x[i, ] <- x[i, ] + (j*pi-pi/2)^{-1}*rt(1, df = 3)*sqrt(2)*sapply(grid, FUN= function(x) sin((j-1/2)*pi*x)  ) # For leverage contamination
  # }
  # }
  
  l1 = x1%*%f1/p
  
  x2 <- matrix(0, (n -floor(0.9*n)), p)
  for(i in 1:(n-floor(0.9*n))){
    # x2[i, ] <- sqrt(2)*(1*pi-pi/2)^{-1}*rnorm(1, 0, 1)*sapply(grid, FUN= function(x) sin((1-1/2)*pi*x)  )
    x2[i, ] <- sqrt(2)*(1*pi-pi/2)^{-1}*rt(1, df = 1)*sapply(grid, FUN= function(x) sin((1-1/2)*pi*x)  )
    for(j in 2:50){
      # x2[i, ] <- x2[i, ] +(j*pi-pi/2)^{-1}*rnorm(1, 0, 1)*sqrt(2)*sapply(grid, FUN= function(x) sin((j-1/2)*pi*x)  )
      # x2[i,] <- 2*x2[i, ]
      x2[i, ] <- x2[i, ] + (j*pi-pi/2)^{-1}*rt(1, df = 1)*sqrt(2)*sapply(grid, FUN= function(x) sin((j-1/2)*pi*x)  ) # For leverage contamination
    }
  }
  
  l2 <- x2%*%f2/p
  
  inv.logit <- function(x) 1/(1+exp(-x))
  y1 <- rbinom(floor(0.9*n), size = 1, prob = inv.logit(l1))
  y2 <- rbinom((n-floor(0.9*n)), size = 1, prob = inv.logit(l2))
  y2 <- 1-y2
  x <- rbind(x1, x2)
  y <- c(y1, y2)
  matplot(t(x), type = "l", lwd = 3, col = "gray")
  
  fit <- dpd.f(x = x, y = y, norder = 4, m = 2)
  fit0 <- dpd.ffa(x = x, y = y, norder = 4, m = 2, alpha = 1e-06)
  fit05 <- dpd.ffa(x = x, y = y, norder = 4, m = 2, alpha = 0.5)
  fit1 <- dpd.ffa(x = x, y = y, norder = 4, m = 2, alpha = 1)
  msq[k] <- mean((f1-fit$est)^2)
  alpha.sel[k] <- fit$alpha
  ic[k] <- fit$ic
  
  msq0[k] <- mean((f1-fit0$est)^2)
  msq05[k] <- mean((f1-fit05$est)^2)
  msq1[k] <- mean((f1-fit1$est)^2)
}

hist(alpha.sel)
# hist(ic)
mean(msq, na.rm = TRUE)
median(msq, na.rm = TRUE) ; median(msq0, na.rm = TRUE) ; median(msq05, na.rm = TRUE) ; median(msq1, na.rm = TRUE)
median(alpha.sel, na.rm = TRUE)

plot(grid, f1, type = "l", lwd = 3)
lines(grid, fit$est, type = "l", lwd = 3, col = "blue")
lines(grid, fit0$est)

