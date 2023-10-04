####################################################### Gait Analysis in Kalogridis (2023) #########################################
rm(list = ls())
# Import the two datasets
df1 <- read.table(file = file.choose(), sep = '\t', header = F)
df2 <- read.table(file = file.choose(), sep = '\t', header = F)

# Merge the training and testing datasets
df <- rbind(df1, df2)
y1 <- df[, 1]
x1 <- df[, -1]

par(mgp = c(3.5, 2.5, 0), mar = c(5.1, 5.1, 4.1, 2.1))
matplot(t(x1), lwd = 3, type = "l", col = "gray", lty = 1, xlab = "", ylab = "", cex.lab = 2.5, cex.axis = 3, ylim = c(-4, 5.5)) ; grid()

# Fit the penalized adaptive and ML estimators
fit.robust <- dpd.f(x = x1, y = y1)
fit.ml <- dpd.f(x1, y1, tuning = 1e-04)
fit.robust$alpha 
# Notice that a large tuning parameter is selected
# Plot the estimators
plot(fit.robust$est, type = "l", lwd = 3, col = "blue",  cex = 2.5, cex.axis = 3, cex.lab = 2.5, xlab = "", ylab = "") ; grid()
lines(fit.ml$est, type = "l", col = "red", lwd = 3, lty = 2)

u.tr <- function(x, alpha){
  n <- length(x)
  m <- round(alpha*n)
  x.s <- sort(x)
  u.tr <- mean(x.s[1:(n-m)])
  return(u.tr)
}

fitted.r <- predict.dpd(x1, fit.robust)
sum(fitted.r$fitted.prob)
u.tr(fitted.r$fitted.prob[y1==1&abs(fit.robust$a.resids) <= 2], alpha = 0.00001) # correctly classified as "normal", normal correctly classified as normal
u.tr(fitted.r$fitted.prob[y1==0&abs(fit.robust$a.resids) <= 2], alpha = 0.1) # abnormal incorrectly classified as "normal"


fitted.ml <- predict.dpd(x1, fit.ml)
sum(fitted.ml$fitted.prob)
u.tr(fitted.ml$fitted.prob[y1==1&abs(fit.robust$a.resids) <= 2], alpha = 0.1) # correctly classified as "normal", normal correctly classified as normal
u.tr(fitted.ml$fitted.prob[y1==0&abs(fit.robust$a.resids) <= 2], alpha = 0.1) # abnormal incorrectly classified as "normal"

which(y1==1)
mean(fitted.r$fitted.prob[which(y1==1)])


# Check Anscombe residuals
hist(fit.robust$a.resids, cex = 2.5, cex.axis = 3, cex.lab = 2.5, xlab = "", ylab = "", main = "") ; grid()
hist(fit.ml$a.resids,  cex = 2.5, cex.axis = 3, cex.lab = 2.5, xlab = "", ylab = "", main = "") ; grid()
# The robust estimator exhibits some very large residuals but the ML estimator only yields modest values of residuals
sum(abs(fit.robust$a.resids)>2) 
# 19 residuals larger than 2 according to the robust estimator
# Plot trajectories of observations with large residuals
matplot(t(x1[abs(fit.robust$a.resids)>2, ]), col = "gray", lwd = 3, type = "l", lty = 1,  cex = 2.5, cex.axis = 3, cex.lab = 2.5,
        xlab = "", ylab = "", ylim = c(-4, 5.5)) ; grid()

# Refit the estimators with outliers removed
y2 <- y1[abs(fit.robust$a.resids)<=2]
x2 <- x1[abs(fit.robust$a.resids)<=2,]
matplot(t(x2), lwd = 3, type = "l", col = "gray", lty = 1)
fit.robust.wo <- dpd.f(x = x2, y = y2 )
# fit.ml.wo <- dpd.ffa(x2, y2, alpha = 1e-04)
fit.ml.wo <- dpd.f(x2, y2, tuning = 1e-04)
fit.robust.wo$alpha
# Notice the low value of the tuning that is selected after the outliers have been removed

plot(fit.robust.wo$est, type = "l", lwd = 3, col = "blue",  cex = 2.5, cex.axis = 3, cex.lab = 2.5, xlab = "", ylab = "") ; grid()
lines(fit.ml.wo$est, type = "l", col = "red", lwd = 3, lty = 2)

# Check residuals again
hist(fit.robust.wo$a.resids, cex = 2.5, cex.axis = 3, cex.lab = 2.5, xlab = "", ylab = "")
hist(fit.ml.wo$a.resids, cex = 2.5, cex.axis = 3, cex.lab = 2.5, xlab = "", ylab = "")
# Identical results after the removal of the outliers

#################################5-fold CV ##################################################################

u.tr <- function(x, alpha){
  n <- length(x)
  m <- round(alpha*n)
  x.s <- sort(x)
  u.tr <- mean(x.s[1:(n-m)])
  return(u.tr)
}

cv.function <- function(x, y, nfolds, int){
  x <- as.matrix(x)
  remainder = dim(x)[1]%%nfolds
  splits <- split(1:dim(x)[1], sample( rep(1:nfolds, times = c(floor(dim(x)[1]/nfolds+remainder), rep( floor(dim(x)[1]/nfolds), nfolds-1 ))  ) ) )
  rmspe.r <- rep(NA, nfolds)
  rmpse.ls <- rep(NA, nfolds)
  
  for(j in 1:nfolds){
    x.test <- x[splits[[j]], ]
    x.train <- x[-splits[[j]], ]
    y.test <- y[splits[[j]]]
    y.train <- y[-splits[[j]]]
    fit.r <- dpd.f(x = x.train, y = y.train)
    fit.ml <- dpd.f(x = x.train, y = y.train, tuning = 1e-04)
    
    x.test <- as.matrix(x.test)
    pred.values.r <- predict.dpd( x.test, fit.r)
    pred.values.ml <- predict.dpd( x.test, fit.ml)
    rmspe.r[j] <- sqrt(u.tr( (y.test-pred.values.r$fitted.prob)^2, 0.1 ))
    rmpse.ls[j] <- sqrt(u.tr( (y.test-pred.values.ml$fitted.prob)^2,0.1 ))
  }
  rmspe.r <- mean(rmspe.r)
  rmpse.ls <- mean(rmpse.ls)
  return(list(rmspe.r, rmpse.ls))
}

