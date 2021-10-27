reg <- read.csv("D:/1101/research/0824/table/regresult.csv")
reg <- reg[,-1]
dfreg <- data.frame(SE=matrix(apply(reg,1,function(x) sd(x)),ncol=1),
          ci=t(apply(reg,1,function(x) quantile(x,c(.025,.975)))))



t3 <- read.csv("D:/1101/research/0824/table/t3result.csv")
t3 <- t3[,-1]
dft3 <- data.frame(SE=matrix(apply(t3,1,function(x) sd(x)),ncol=1),
           ci=t(apply(t3,1,function(x) quantile(x,c(.025,.975)))))


result <- list(regression=dfreg,
               table3=dft3)
library(data.table)
outputfile <- "result.csv" #output file name
sep <- "," #define the separator (related to format of the output file)
for(i in names(result)){
  fwrite(list(i), file=outputfile, sep=sep, append=T) #write names of the list elements
  ele <- result[[i]]
  fwrite(ele, file=outputfile, sep=sep, append=T, col.names=T) 
  fwrite(list(NA), file=outputfile, append=T) #add an empty row to separate elements
}
