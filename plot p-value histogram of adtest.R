#--------------------------- plot p-value of ADtest ----------------------------
#
result <- apply(rbind(ad.p1,ad.p2,ad.p3), 1, function(x) c(mean(x),sd(x),quantile(x, probs = c(0.025, 0.975)),
                                                           HPDinterval(as.mcmc(x),prob = 0.95)))
rownames(result) <- c("mean","sd","2.5%","97.5%",
                      "95% credible interval (lower bound)","95% credible interval (upper bound)")
colnames(result) <- c("1C","3C","5C")
write.csv(result,"adtest_result_0111_3chains_accept (第一種演算法).csv")
#-----------------------
# 各應力下的pvalue直方圖
#-----------------------
par(mfrow=c(3,1))
hist(ad.p1,prob = TRUE,main="p-value of AD test under 1C stress")
# lines(density(ad.p1, adjust =1.5),
#       lwd = 2, # thickness of line
#       col = "red")
abline(v=0.05, lwd = 2, lty=1, col = "blue")
mapply(function(x){abline(v=x,lwd=2,lty=2,col="green")},result[3:4,1])
mapply(function(x){abline(v=x,lwd=2,lty=3,col="red")},result[5:6,1])


hist(ad.p2,prob = TRUE,main="p-value of AD test under 3C stress")
# lines(density(ad.p2, adjust =1.5),
#       lwd = 2, # thickness of line
#       col = "red")
abline(v=0.05, lwd = 2, lty=1, col = "blue")
mapply(function(x){abline(v=x,lwd=2,lty=2,col="green")},result[3:4,2])
mapply(function(x){abline(v=x,lwd=2,lty=3,col="red")},result[5:6,2])


hist(ad.p3,prob = TRUE,main="p-value of AD test under 5C stress")
# lines(density(ad.p3, adjust =1.5),
#       lwd = 2, # thickness of line
#       col = "red")
abline(v=0.05, lwd = 2, lty=1, col = "blue")
mapply(function(x){abline(v=x,lwd=2,lty=2,col="green")},result[3:4,3])
mapply(function(x){abline(v=x,lwd=2,lty=3,col="red")},result[5:6,3])

legend("topright",c("p-value=0.05","95% quantile","95% credible interval"),
       col=c("blue","green","red"),lty=c(1,2,3),lwd=2,box.lty = 0)
par(mfrow=c(1,1))