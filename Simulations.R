##################################################### Clean data #########################################################
##########################################################################################################################

nrep <- 1000
n <- 400
p <- 200
grid <- seq(1/p, 1-1/p, len = p)

m1 <- m2 <- m3 <- m4 <- matrix(NA, nrow = p, ncol = nrep)
msq1 = msq2 = msq3 = msq4 =  rep(NA, nrep)
inv.logit <- function(x) 1/(1+exp(-x))

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
  fit3 <- dpd.ffa(x = x, y = y, norder = 4, m = 2, alpha = 1)
  fit4 <- dpd.ffa(x = x, y = y, norder = 4, m = 2, alpha = 2)

  msq1[k] <- mean((f1-fit1$est)^2)
  msq2[k] <- mean((f1-fit2$est)^2)
  msq3[k] <- mean((f1-fit3$est)^2)
  msq4[k] <- mean((f1-fit4$est)^2)
  
  m1[, k] <- fit1$est
  m2[, k] <- fit2$est
  m3[, k] <- fit3$est
  m4[, k] <- fit4$est
  
}

matplot(grid, m1, type = "l", col = "gray", lwd = 3, cex = 2.5, cex.axis = 3, cex.lab = 2.5, xlab = "", ylab = ""); lines(grid, f1, lwd = 3, col = "black") ;grid()
matplot(grid, m2, type = "l", col = "gray", lwd = 3,  cex = 2.5, cex.axis = 3, cex.lab = 2.5, xlab = "", ylab = ""); lines(grid, f1, lwd = 3, col = "black") ; grid()

median(msq1) ; median(msq2) ; median(msq3) ; median(msq4)

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


###################################################### Contamination #######################################################
############################################################################################################################

nrep <- 1000
n <- 400
p <- 200
eps <- 0.01
grid <- seq(1/p, 1-1/p, len = p)

msq1 = msq2 = msq3 = msq4 =  rep(NA, nrep)
m1 <- m2 <- m3 <- m4 <- matrix(NA, nrow = p, ncol = nrep)
inv.logit <- function(x) 1/(1+exp(-x))

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
  samp <- sample(1:n, size = floor(n*eps))
  x[samp, ] <- 5*x[samp, ]

  y0 = (x%*%f1)/p
  y <- rbinom(n, size = 1, prob = inv.logit(y0))
  y[samp]<- 1- y[samp]
  
  fit1 <- dpd.f(x = x, y = y, norder = 4, m = 2)
  fit2 <- p.ml(x, y, norder = 4, nbasis = 30, m = 2)
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

matplot(grid, m1, type = "l", col = "gray", lwd = 3)
matplot(grid, m2, type = "l", col = "gray", lwd = 3); lines(grid, f1, lwd = 3, col = "black")

median(msq1) ; median(msq2) ; median(msq3) ; median(msq4)

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

boot.median(msq1)
boot.median(msq2)
boot.median(msq3)
boot.median(msq4)
