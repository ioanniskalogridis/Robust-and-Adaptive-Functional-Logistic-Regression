####################################################### Gait Analysis in Kalogridis (2023) #########################################

# Import the two datasets
#df1 <- read.table(file = file.choose(), sep = '\t', header = F)
#df2 <- read.table(file = file.choose(), sep = '\t', header = F)

# Merge the training and testing data
df <- rbind(df1, df2)
y1 <- df[, 1]
x1 <- df[, -1]
matplot(t(x1), lwd = 3, type = "l", col = "gray", lty = 1)
fit1 <- dpd.f(x = x1, y = y1, norder = 4, m = 2 )
fit1.ml <- dpd.ffa(x1, y1, alpha = 1e-04)
fit1$alpha
fit1$lambda
plot(fit1$est, type = "l", lwd = 3, col = "blue",  cex = 2.5, cex.axis = 3, cex.lab = 2.5, xlab = "", ylab = "") ; grid()
lines(fit1.ml$est, type = "l", col = "red", lwd = 3, lty = 2)
hist(fit1$a.resids, breaks = 20)
hist(fit1.ml$a.resids)
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