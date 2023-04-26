####################################################### Gait Analysis in Kalogridis (2023) #########################################

# Import the two datasets
#df1 <- read.table(file = file.choose(), sep = '\t', header = F)
#df2 <- read.table(file = file.choose(), sep = '\t', header = F)

# Merge the training and testing datasets
df <- rbind(df1, df2)
y1 <- df[, 1]
x1 <- df[, -1]

par(mgp = c(3.5, 2.5, 0), mar = c(5.1, 5.1, 4.1, 2.1))
matplot(t(x1), lwd = 3, type = "l", col = "gray", lty = 1, xlab = "", ylab = "")

# Fit the penalized adaptive and ML estimators
fit.in <- dpd.f(x = x1, y = y1)
fit.ml.in <- dpd.ffa(x1, y1)
fit.in$alpha # Notice that a high tuning is selected
# Plot the estimators
plot(fit.in$est, type = "l", lwd = 3, col = "blue",  cex = 2.5, cex.axis = 3, cex.lab = 2.5, xlab = "", ylab = "") ; grid()
lines(fit.ml.in$est, type = "l", col = "red", lwd = 3, lty = 2)

# Check residuals
hist(fit1$a.resids)
hist(fit1.ml$a.resids)
# The robust estimator exhibits some very large residuals but the ML estimator only yields modest values of residuals
sum(abs(fit1$a.resids)>2)
sum(abs(fit1.ml$a.resids)>2)


matplot(t(x1), lwd = 3, type = "l", col = "gray", lty = 1,  cex = 2.5, cex.axis = 3, cex.lab = 2.5, xlab = "", ylab = "") ; grid()
matplot(t(x1[abs(fit1$a.resids)>2, ]), col = "gray", lwd = 3, type = "l", lty = 1,  cex = 2.5, cex.axis = 3, cex.lab = 2.5, xlab = "", ylab = "") ; grid()
y1[abs(fit1$a.resids)>2]

y2 <- y1[abs(fit1$a.resids)<=2]
x2 <- x1[abs(fit1$a.resids)<=2,]
matplot(t(x2), lwd = 3, type = "l", col = "gray", lty = 1)
fit2 <- dpd.f(x = x2, y = y2, norder = 4, m = 2 )
fit.ml2 <- dpd.ffa(x2, y2, alpha = 1e-04)
fit2$alpha
fit2$lambda
plot(fit2$est, type = "l", lwd = 3, col = "blue",  cex = 2.5, cex.axis = 3, cex.lab = 2.5, xlab = "", ylab = "") ; grid()
lines(fit.ml2$est, type = "l", col = "red", lwd = 3, lty = 2)
# lines(fit.ml2$est, type = "l", col = "red", lwd = 3)
hist(fit2$a.resids)
hist(fit.ml2$a.resids)
sum(abs(fit2$a.resids)>2.6)