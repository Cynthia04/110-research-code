#-------------------------------------------------------------------------------
#plot 95% CI
#-------------------------------------------------------------------------------
TP <- as.vector(t(t3))
cpal <-function(x){
  CIs <- cbind(par[x,]-1.96*par[x+7,],par[x,]+1.96*par[x+7,])
  ID <- which(!(CIs[1:100, 1] <= TP[x] & TP[x] <= CIs[1:100, 2]))
  cp <- 1-length(ID)/100
  al <- mean(CIs[,2]-CIs[,1])
  return(c(cp,al))
} 
write.csv(apply(matrix(1:6,ncol = 1),1,cpal),"weibullcpal.csv")


par(mfrow=c(3,2))

###a0
CIs <- cbind(par[1,]-1.96*par[8,],par[1,]+1.96*par[8,])
ID <- which(!(CIs[1:100, 1] <= TP[1] & TP[1] <= CIs[1:100, 2]))
1-length(ID)/100
# initialize the plot
plot(0, 
     xlim = c(0.91,0.92), 
     ylim = c(1, 100), 
     ylab = "Sample", 
     xlab = expression(a[0]), 
     main = expression(paste("95% Confidence Intervals of ",a[0])))

# set up color vector
colors <- rep(gray(0.6), 100)
colors[ID] <- "red"

# draw reference line at mu=5
abline(v = TP[1], lty = 2)

# add horizontal bars representing the CIs
for(j in 1:100) {
  
  lines(c(CIs[j, 1], CIs[j, 2]), 
        c(j, j), 
        col = colors[j], 
        lwd = 2)
  
}
###   a1-------------------------
CIs <- cbind(par[2,]-1.96*par[9,],par[2,]+1.96*par[9,])
ID <- which(!(CIs[1:100, 1] <= TP[2] & TP[2] <= CIs[1:100, 2]))
# initialize the plot
plot(0, 
     xlim = c(0.005,0.015), 
     ylim = c(1, 100), 
     ylab = "Sample", 
     xlab = expression(a[1]), 
     main = expression(paste("95% Confidence Intervals of ",a[1])))

# set up color vector
colors <- rep(gray(0.6), 100)
colors[ID] <- "red"

# draw reference line at mu=5
abline(v = TP[2], lty = 2)

# add horizontal bars representing the CIs
for(j in 1:100) {
  
  lines(c(CIs[j, 1], CIs[j, 2]), 
        c(j, j), 
        col = colors[j], 
        lwd = 2)
  
}


### b0-------------------------------------------
CIs <- cbind(par[3,]-1.96*par[10,],par[3,]+1.96*par[10,])
ID <- which(!(CIs[1:100, 1] <= TP[3] & TP[3] <= CIs[1:100, 2]))
# initialize the plot
plot(0, 
     xlim = c(0.0205,0.0209), 
     ylim = c(1, 100), 
     ylab = "Sample", 
     xlab = expression(b[0]), 
     main = expression(paste("95% Confidence Intervals of ",b[0])))

# set up color vector
colors <- rep(gray(0.6), 100)
colors[ID] <- "red"

# draw reference line at mu=5
abline(v = TP[3], lty = 2)

# add horizontal bars representing the CIs
for(j in 1:100) {
  
  lines(c(CIs[j, 1], CIs[j, 2]), 
        c(j, j), 
        col = colors[j], 
        lwd = 2)
  
}


#b1
CIs <- cbind(par[4,]-1.96*par[11,],par[4,]+1.96*par[11,])
ID <- which(!(CIs[1:100, 1] <= TP[4] & TP[4] <= CIs[1:100, 2]))
# initialize the plot
plot(0, 
     xlim = c(0.0009,0.0014), 
     ylim = c(1, 100), 
     ylab = "Sample", 
     xlab = expression(b[1]), 
     main = expression(paste("95% Confidence Intervals of ",b[1])))

# set up color vector
colors <- rep(gray(0.6), 100)
colors[ID] <- "red"

# draw reference line at mu=5
abline(v = TP[4], lty = 2)

# add horizontal bars representing the CIs
for(j in 1:100) {
  
  lines(c(CIs[j, 1], CIs[j, 2]), 
        c(j, j), 
        col = colors[j], 
        lwd = 2)
  
}


#beta0
CIs <- cbind(par[5,]-1.96*par[12,],par[5,]+1.96*par[12,])
ID <- which(!(CIs[1:100, 1] <= TP[5] & TP[5] <= CIs[1:100, 2]))
min(CIs[,1]);max(CIs[,2])
# initialize the plot
plot(0, 
     xlim = c(min(CIs[,1]),max(CIs[,2])), 
     ylim = c(1, 100), 
     ylab = "Sample", 
     xlab = expression(beta[0]), 
     main = expression(paste("95% Confidence Intervals of ",beta[0])))

# set up color vector
colors <- rep(gray(0.6), 100)
colors[ID] <- "red"

# draw reference line at mu=5
abline(v = TP[5], lty = 2)

# add horizontal bars representing the CIs
for(j in 1:100) {
  
  lines(c(CIs[j, 1], CIs[j, 2]), 
        c(j, j), 
        col = colors[j], 
        lwd = 2)
  
}


#beta1
CIs <- cbind(par[6,]-1.96*par[13,],par[6,]+1.96*par[13,])
ID <- which(!(CIs[1:100, 1] <= TP[6] & TP[6] <= CIs[1:100, 2]))

# initialize the plot
plot(0, 
     xlim = c(min(CIs[,1]),max(CIs[,2])), 
     ylim = c(1, 100), 
     ylab = "Sample", 
     xlab = expression(beta[1]), 
     main = expression(paste("95% Confidence Intervals of ",beta[1])))

# set up color vector
colors <- rep(gray(0.6), 100)
colors[ID] <- "red"

# draw reference line at mu=5
abline(v = TP[6], lty = 2)

# add horizontal bars representing the CIs
for(j in 1:100) {
  
  lines(c(CIs[j, 1], CIs[j, 2]), 
        c(j, j), 
        col = colors[j], 
        lwd = 2)
  
}

