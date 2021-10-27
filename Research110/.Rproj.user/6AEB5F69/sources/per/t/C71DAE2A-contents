par <- read.csv("par0final.csv")
par<- par[,-1]
TP <- c(3,5,20)
par <- as.matrix(par)
table <- rbind(apply(par[1:3,],1,function(x){return(c(mean(x),sd(x),quantile(x,c(.025,.975))))}),
               apply(matrix(1:3,ncol=1),1,function(x){mean((par[x,]-TP[x])^2)}))
colnames(table) <- c("a0","b0","beta0")
rownames(table) <- c("mean","sd","2.5%","97.5%","MSE")
table <- round(t(table),5)
par(mfrow=c(3,1))
hist(as.numeric(par[1,]),main=expression(paste("Histogram of ",a[0],"(1000 simulation)")),
     xlab = expression(a[0]))
abline(v=3,col="red",lty=2,lwd=2)
print(paste("True=",TP[1],\n," mean=",round(table[1,1],5)))
print(1 )
labels <- "True=3, \n mean=3.00154, \n sd=0.06201, \n MSE=0.00381"
legend("topleft",y.intersp=0.01,cex=1.5,bty="n", 
       labels, col=colors)
hist(as.numeric(par[2,]),main=expression(paste("Histogram of ",b[0],"(1000 simulation)")),
     xlab = expression(b[0]))
abline(v=5,col="red",lty=2,lwd=2)
labels <- "True=5, \n mean=4.99933, \n sd=0.02931, \n MSE=0.00085 "
legend("topleft",y.intersp=0.01,cex=1.5,bty="n", 
       labels, col=colors)
hist(as.numeric(par[3,]),main=expression(paste("Histogram of ",beta[0],"(1000 simulation)")),
     xlab = expression(beta[0]))
abline(v=20,col="red",lty=2,lwd=2)
labels <- "True=20, \n mean=20.35515, \n sd=1.39275, \n MSE=2.04608 "
legend("topleft",x.intersp=0.00000001,y.intersp=0.00001,cex=1.5,bty="n", 
       labels, col=colors)
par(mfrow=c(1,1))



################################################################################################


#-------------------------------------------------------------------------------
#plot 95% CI
#-------------------------------------------------------------------------------
N=98
par <- as.matrix(par)
cpal <-function(x){
  CIs <- cbind(par[x,]-1.96*par[x+4,],par[x,]+1.96*par[x+4,])
  ID <- which(!(CIs[1:N, 1] <= TP[x] & TP[x] <= CIs[1:N, 2]))
  cp <- 1-length(ID)/N
  al <- mean(CIs[,2]-CIs[,1])
  return(c(cp,al))
} 

write.csv(apply(matrix(1:3,ncol = 1),1,cpal),"weibullcpal0.csv")

par(mfrow=c(3,1))
###a0
CIs <- cbind(par[1,]-1.96*par[5,],par[1,]+1.96*par[5,])
ID <- which(!(CIs[1:N, 1] <= TP[1] & TP[1] <= CIs[1:N, 2]))
# initialize the plot
plot(0, 
     xlim = c(2.7,3.4), 
     ylim = c(1, N), 
     ylab = "Sample", 
     xlab = expression(a[0]), 
     main = expression(paste("95% Confidence Intervals of ",a[0])),
     cex.lab=2, cex.axis=1.5, cex.main=2, cex.sub=1.5)

# set up color vector
colors <- rep(gray(0.6), N)
colors[ID] <- "red"


for(j in 1:N) {
  
  lines(c(CIs[j, 1], CIs[j, 2]), 
        c(j, j), 
        col = colors[j], 
        lwd = 2)
}
abline(v = TP[1], lty = 2,col="blue")

### b0-------------------------------------------
CIs <- cbind(par[2,]-1.96*par[6,],par[2,]+1.96*par[6,])
ID <- which(!(CIs[1:N, 1] <= TP[2] & TP[2] <= CIs[1:N, 2]))
# initialize the plot
plot(0, 
     xlim = c(4.8,5.2), 
     ylim = c(1, N), 
     ylab = "Sample", 
     xlab = expression(b[0]), 
     main = expression(paste("95% Confidence Intervals of ",b[0])),
     cex.lab=2, cex.axis=1.5, cex.main=2, cex.sub=1.5)

# set up color vector
colors <- rep(gray(0.6), N)
colors[ID] <- "red"
for(j in 1:N) {
  
  lines(c(CIs[j, 1], CIs[j, 2]), 
        c(j, j), 
        col = colors[j], 
        lwd = 2)
}
abline(v = TP[2], lty = 2,col="blue")



#beta0
CIs <- cbind(par[3,]-1.96*par[7,],par[3,]+1.96*par[7,])
ID <- which(!(CIs[1:N, 1] <= TP[3] & TP[3] <= CIs[1:N, 2]))
min(CIs[,1]);max(CIs[,2])
# initialize the plot
plot(0, 
     xlim = c(14,28), 
     ylim = c(1, N), 
     ylab = "Sample", 
     xlab = expression(beta[0]), 
     main = expression(paste("95% Confidence Intervals of ",beta[0])),
     cex.lab=2, cex.axis=1.5, cex.main=2, cex.sub=1.5)

# set up color vector
colors <- rep(gray(0.6), N)
colors[ID] <- "red"
# add horizontal bars representing the CIs
for(j in 1:N) {
  
  lines(c(CIs[j, 1], CIs[j, 2]), 
        c(j, j), 
        col = colors[j], 
        lwd = 2)
}
abline(v = TP[3], lty = 2,col="blue")


