par1 <- read.csv("par1.csv")
par2 <- read.csv("par2.csv")
par3 <- read.csv("par3.csv")
par4 <- read.csv("par4.csv")
par5 <- read.csv("par5.csv")
par6 <- read.csv("par6.csv")

par <- cbind(par1[,-c(1,2)],par2[,-c(1,2)],par3[,-c(1,2)],par4[,-c(1,2)],par5[,-c(1,2)],par6[,-c(1,2,3)])

par(mfrow=c(3,2))
hist(as.numeric(par[1,]),main=expression(paste("Histogram of ",a[0],"(1000 simulation)")),
     xlab = expression(a[0]))
abline(v=3,col="red",lty=2,lwd=2)
labels <- "True=3, \n mean=3.00162, \n sd=0.06983, \n MSE=0.00487"
legend("topleft",y.intersp=0.01,cex=1.5,bty="n", 
       labels, col=colors)
hist(as.numeric(par[2,]),main=expression(paste("Histogram of ",a[1],"(1000 simulation)")),
     xlab = expression(a[1]))
abline(v=2,col="red",lty=2,lwd=2)
labels <- "True=2, \n mean=2.00355 , \n sd=0.14098, \n MSE=0.01987"
legend("topleft",y.intersp=0.01,cex=1.5,bty="n", 
       labels, col=colors)
hist(as.numeric(par[3,]),main=expression(paste("Histogram of ",b[0],"(1000 simulation)")),
     xlab = expression(b[0]))
abline(v=5,col="red",lty=2,lwd=2)
labels <- "True=5, \n mean=4.99938, \n sd=0.03337, \n MSE=0.00111 "
legend("topleft",y.intersp=0.01,cex=1.5,bty="n", 
       labels, col=colors)
hist(as.numeric(par[4,]),main=expression(paste("Histogram of ",b[1],"(1000 simulation)")),
     xlab = expression(b[1]))
abline(v=4,col="red",lty=2,lwd=2)
labels <- "True=4, \n mean=3.99790, \n sd=0.06987, \n MSE=0.00488 "
legend("topleft",y.intersp=0.01,cex=1.5,bty="n", 
       labels, col=colors)
hist(as.numeric(par[5,]),main=expression(paste("Histogram of ",beta[0],"(1000 simulation)")),
     xlab = expression(beta[0]))
abline(v=20,col="red",lty=2,lwd=2)
labels <- "True=20, \n mean=20.23702, \n sd=1.36473, \n MSE=1.91679 "
legend("topleft",x.intersp=0.00000001,y.intersp=0.00001,cex=1.5,bty="n", 
       labels, col=colors)
hist(as.numeric(par[6,]),main=expression(paste("Histogram of ",beta[1],"(1000 simulation)")),
     xlab = expression(beta[1]))
abline(v=-5,col="red",lty=2,lwd=2)
labels <- "True=-5, \n mean=-5.08808, \n sd=1.84041, \n MSE=3.39149"
legend("topleft",y.intersp=0.01,cex=1.5,bty="n", 
       labels, col=colors)
par(mfrow=c(1,1))

apply(par[1:6,],1,function(x){return(c(mean(x),sd(x),quantile(x,c(.025,.975))))})


################################################################################################


#-------------------------------------------------------------------------------
#plot 95% CI
#-------------------------------------------------------------------------------
par <- as.matrix(par)
t3 <- matrix(c(3,
               2,
               5,
               4,
               20,
               -5),ncol=2,byrow = T)
TP <- as.vector(t(t3))
cpal <-function(x){
  CIs <- cbind(par[x,]-1.96*par[x+7,],par[x,]+1.96*par[x+7,])
  ID <- which(!(CIs[1:1000, 1] <= TP[x] & TP[x] <= CIs[1:1000, 2]))
  cp <- 1-length(ID)/1000
  al <- mean(CIs[,2]-CIs[,1])
  return(c(cp,al))
} 

write.csv(apply(matrix(1:6,ncol = 1),1,cpal),"weibullcpal.csv")


par(mfrow=c(3,2))
###a0
CIs <- cbind(par[1,]-1.96*par[8,],par[1,]+1.96*par[8,])
ID <- which(!(CIs[1:1000, 1] <= TP[1] & TP[1] <= CIs[1:1000, 2]))
# initialize the plot
plot(0, 
     xlim = c(2.5,3.4), 
     ylim = c(1, 1000), 
     ylab = "Sample", 
     xlab = expression(a[0]), 
     main = expression(paste("95% Confidence Intervals of ",a[0])),
     cex.lab=2, cex.axis=1.5, cex.main=2, cex.sub=1.5)

# set up color vector
colors <- rep(gray(0.6), 1000)
colors[ID] <- "red"

for(j in 1:1000) {
  
  lines(c(CIs[j, 1], CIs[j, 2]), 
        c(j, j), 
        col = colors[j], 
        lwd = 2)
}
abline(v = TP[1], lty = 2,col="blue")

###   a1-------------------------
CIs <- cbind(par[2,]-1.96*par[9,],par[2,]+1.96*par[9,])
ID <- which(!(CIs[1:1000, 1] <= TP[2] & TP[2] <= CIs[1:1000, 2]))
# initialize the plot
plot(0, 
     xlim = c(1.2,2.8), 
     ylim = c(1, 1000), 
     ylab = "Sample", 
     xlab = expression(a[1]), 
     main = expression(paste("95% Confidence Intervals of ",a[1])),
     cex.lab=2, cex.axis=1.5, cex.main=2, cex.sub=1.5)

# set up color vector
colors <- rep(gray(0.6), 1000)
colors[ID] <- "red"

for(j in 1:1000) {
  
  lines(c(CIs[j, 1], CIs[j, 2]), 
        c(j, j), 
        col = colors[j], 
        lwd = 2)
}
abline(v = TP[2], lty = 2,col="blue")


### b0-------------------------------------------
CIs <- cbind(par[3,]-1.96*par[10,],par[3,]+1.96*par[10,])
ID <- which(!(CIs[1:1000, 1] <= TP[3] & TP[3] <= CIs[1:1000, 2]))
# initialize the plot
plot(0, 
     xlim = c(4.8,5.2), 
     ylim = c(1, 1000), 
     ylab = "Sample", 
     xlab = expression(b[0]), 
     main = expression(paste("95% Confidence Intervals of ",b[0])),
     cex.lab=2, cex.axis=1.5, cex.main=2, cex.sub=1.5)

# set up color vector
colors <- rep(gray(0.6), 1000)
colors[ID] <- "red"
for(j in 1:1000) {
  
  lines(c(CIs[j, 1], CIs[j, 2]), 
        c(j, j), 
        col = colors[j], 
        lwd = 2)
}
abline(v = TP[3], lty = 2,col="blue")



#b1
CIs <- cbind(par[4,]-1.96*par[11,],par[4,]+1.96*par[11,])
ID <- which(!(CIs[1:1000, 1] <= TP[4] & TP[4] <= CIs[1:1000, 2]))
# initialize the plot
plot(0, 
     xlim = c(3.6,4.4), 
     ylim = c(1, 1000), 
     ylab = "Sample", 
     xlab = expression(b[1]), 
     main = expression(paste("95% Confidence Intervals of ",b[1])),
     cex.lab=2, cex.axis=1.5, cex.main=2, cex.sub=1.5)

# set up color vector
colors <- rep(gray(0.6), 1000)
colors[ID] <- "red"
for(j in 1:1000) {
  
  lines(c(CIs[j, 1], CIs[j, 2]), 
        c(j, j), 
        col = colors[j], 
        lwd = 2)
}
abline(v = TP[4], lty = 2,col="blue")



#beta0
CIs <- cbind(par[5,]-1.96*par[12,],par[5,]+1.96*par[12,])
ID <- which(!(CIs[1:1000, 1] <= TP[5] & TP[5] <= CIs[1:1000, 2]))
min(CIs[,1]);max(CIs[,2])
# initialize the plot
plot(0, 
     xlim = c(14,28), 
     ylim = c(1, 1000), 
     ylab = "Sample", 
     xlab = expression(beta[0]), 
     main = expression(paste("95% Confidence Intervals of ",beta[0])),
     cex.lab=2, cex.axis=1.5, cex.main=2, cex.sub=1.5)

# set up color vector
colors <- rep(gray(0.6), 1000)
colors[ID] <- "red"

for(j in 1:1000) {
  
  lines(c(CIs[j, 1], CIs[j, 2]), 
        c(j, j), 
        col = colors[j], 
        lwd = 2)
}
abline(v = TP[5], lty = 2,col="blue")


#beta1
CIs <- cbind(par[6,]-1.96*par[13,],par[6,]+1.96*par[13,])
ID <- which(!(CIs[1:1000, 1] <= TP[6] & TP[6] <= CIs[1:1000, 2]))

# initialize the plot
plot(0, 
     xlim = c(-15,5), 
     ylim = c(1, 1000), 
     ylab = "Sample", 
     xlab = expression(beta[1]), 
     main = expression(paste("95% Confidence Intervals of ",beta[1])),
     cex.lab=2, cex.axis=1.5, cex.main=2, cex.sub=1.5)

# set up color vector
colors <- rep(gray(0.6), 1000)
colors[ID] <- "red"
for(j in 1:1000) {
  
  lines(c(CIs[j, 1], CIs[j, 2]), 
        c(j, j), 
        col = colors[j], 
        lwd = 2)
}
abline(v = TP[6], lty = 2,col="blue")



TP
apply(matrix(1:6,ncol=1),1,function(x){mean((par[x,]-TP[x])^2)})




