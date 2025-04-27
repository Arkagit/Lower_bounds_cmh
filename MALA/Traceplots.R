set.seed(1234)

source("CMH.R")

library(car)
library(ggplot2)
library(plotly)
#install.packages(c("plotly", "webshot", "htmlwidgets"))
#webshot::install_phantomjs()  # Run this once
library(htmlwidgets)
library(webshot)
library(calculus)



rep = 1e4 # number of repetition

################ Gibbs
X1 <- rep(2, length = rep)
X2 = rep(-2, length = rep) 

for(i in 2:rep){
  out1 = gibbs_output(X1[i-1], X2[i-1])
  X1[i] = out1[1]
  X2[i] = out1[2]
}
 
Out_gibbs = cbind(X1, X2)

pdf("TracePlot_X1G.pdf", height = 6, width = 8)
plot.ts(X1, type = "l", ylim = c(-2.5, 1.5))
dev.off()

pdf("TracePlot_X2G.pdf", height = 6, width = 8)
plot.ts(X2, type = "l")
dev.off()

################### CMH
h = c(0.001, 1, 100) # proposal standard deviation

X1 = matrix(2, nrow = rep, ncol = length(h))
X2 = matrix(-2, nrow = rep, ncol = length(h))
A = matrix(0, nrow = rep, ncol = length(h))

for(j in 1:length(h)){
  for(i in 2:rep){
    Out = CMH_output(X1[i-1, j], X2[i-1, j], h[j])
    X1[i, j] = Out[1]
    X2[i, j] = Out[2]
    A[i, j] = Out[3]
  }
}

data1 = data.frame(cbind(X1[,1], X2[,1]))
data2 = data.frame(cbind(X1[,2], X2[,2]))
data3 = data.frame(cbind(X1[,3], X2[,3]))

pdf("TracePlot_X1.pdf", height = 6, width = 8)
par(mfrow = c(1, 3))  ## set the layout to be 3 by 3
sapply(1:length(h), function(i) plot.ts(X1[,i], ylim = c(-3, 3), type = "l", ylab = "X1", main = paste("h = ", h[i])))
dev.off()

pdf("TracePlot_X2.pdf", height = 6, width = 8)
par(mfrow = c(1, 3))  ## set the layout to be 3 by 3
sapply(1:length(h), function(i) plot.ts(X2[,i], ylim = c(-3, 3), ylab = "X2", main = paste("h = ", h[i]), type = "l"))
dev.off()

############################################
############################################ CMH-2
h = c(0.001, 1, 100) # proposal standard deviation

X1 = matrix(2, nrow = rep, ncol = length(h))
X2 = matrix(-2, nrow = rep, ncol = length(h))
A = matrix(0, nrow = rep, ncol = length(h))

for(j in 1:length(h)){
  for(i in 2:rep){
    Out = CMH2_output(X1[i-1, j], X2[i-1, j], h[j])
    X1[i, j] = Out[1]
    X2[i, j] = Out[2]
    A[i, j] = Out[3]
  }
}

data1 = data.frame(cbind(X1[,1], X2[,1]))
data2 = data.frame(cbind(X1[,2], X2[,2]))
data3 = data.frame(cbind(X1[,3], X2[,3]))

pdf("TracePlot2_X1.pdf", height = 6, width = 8)
par(mfrow = c(1, 3))  ## set the layout to be 3 by 3
sapply(1:length(h), function(i) plot.ts(X1[,i], ylim = c(-3, 3), ylab = "X1", main = paste("h = ", h[i]), type = "l"))
dev.off()

pdf("TracePlot2_X2.pdf", height = 6, width = 8)
par(mfrow = c(1, 3))  ## set the layout to be 3 by 3
sapply(1:length(h), function(i) plot.ts(X2[,i], ylim = c(-3, 3), ylab = "X2", main = paste("h = ", h[i]), type = "l"))
dev.off()
