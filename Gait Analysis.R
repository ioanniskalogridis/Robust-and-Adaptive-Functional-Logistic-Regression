####################################################### Gait Analysis in Kalogridis (2023) #########################################

# Import the two datasets
df1 <- read.table(file = file.choose(), sep = '\t', header = F)
df2 <- read.table(file = file.choose(), sep = '\t', header = F)

# Merge the training and testing datasets
df <- rbind(df1, df2)
y1 <- df[, 1]
x1 <- df[, -1]

par(mgp = c(3.5, 2.5, 0), mar = c(5.1, 5.1, 4.1, 2.1))
matplot(t(x1), lwd = 3, type = "l", col = "gray", lty = 1, xlab = "", ylab = "")

# Fit the penalized adaptive and ML estimators
fit.robust <- dpd.f(x = x1, y = y1)
fit.ml <- dpd.ffa(x1, y1)
fit.robust$alpha # Notice that a high tuning is selected
# Plot the estimators
plot(fit.robust$est, type = "l", lwd = 3, col = "blue",  cex = 2.5, cex.axis = 3, cex.lab = 2.5, xlab = "", ylab = "", ylim = c(-100, 100)) ; grid()
lines(fit.ml$est, type = "l", col = "red", lwd = 3, lty = 2)

# Check Anscombe residuals
hist(fit.robust$a.resids, cex = 2.5, cex.axis = 3, cex.lab = 2.5, xlab = "", ylab = "")
hist(fit.ml$a.resids,  cex = 2.5, cex.axis = 3, cex.lab = 2.5, xlab = "", ylab = "")
# The robust estimator exhibits some very large residuals but the ML estimator only yields modest values of residuals
sum(abs(fit.robust$a.resids)>2) 
# 19 residuals larger than 2 according to the robust estimator
# Plot trajectories of observations with large residuals
matplot(t(x1[abs(fit.robust$a.resids)>2, ]), col = "gray", lwd = 3, type = "l", lty = 1,  cex = 2.5, cex.axis = 3, cex.lab = 2.5, xlab = "", ylab = "") ; grid()

# Refit the estimators with outliers removed
y2 <- y1[abs(fit.robust$a.resids)<=2]
x2 <- x1[abs(fit.robust$a.resids)<=2,]
matplot(t(x2), lwd = 3, type = "l", col = "gray", lty = 1)
fit.robust.wo <- dpd.f(x = x2, y = y2 )
fit.ml.wo <- dpd.ffa(x2, y2, alpha = 1e-04)
fit.robust.wo$alpha
# Notice the low value of the tuning that is selected after the outliers have been removed

plot(fit.robust.wo$est, type = "l", lwd = 3, col = "blue",  cex = 2.5, cex.axis = 3, cex.lab = 2.5, xlab = "", ylab = "",
     ylim = c(-100, 100)) ; grid()
lines(fit.ml.wo$est, type = "l", col = "red", lwd = 3, lty = 2)

# Check residuals again
hist(fit.robust.wo$a.resids, cex = 2.5, cex.axis = 3, cex.lab = 2.5, xlab = "", ylab = "")
hist(fit.ml.wo$a.resids, cex = 2.5, cex.axis = 3, cex.lab = 2.5, xlab = "", ylab = "")
# Identical results after the removal of the outliers