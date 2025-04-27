###################### 2-component CMH
######## (X|Y) ~ Normal(0, 1/(2(1+y2))) Gibbs step
######## (Y|X) ~ Normal(0, 1/(2(1+x2))) MH step by Normal(y, h)
set.seed(1234)

source("CMH_func.R")

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

fig = plot_ly(x = X1, y = X2, type = "histogram2dcontour") %>%
  layout(
    title = list(text = "Gibbs", x = 0.5, font = list(size = 18, color = "black")),
    xaxis = list(title = "X1"),
    yaxis = list(title = "X2")
  )
saveWidget(fig, paste("temp_plotG.html"), selfcontained = FALSE)
webshot(paste("temp_plotG.html"), paste("contour_gibbs.pdf"))

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


##################################################
################################################## CMH -updated
sequence = 1e2
h = c(0.001, 1, 100) # proposal standard deviation

l2 = seq(-5, 5, length.out = sequence)# Grid

lb_list = list()

for(i in 1:length(h)){
  lb_list[[i]] = rep(0, length = sequence)
}

for(i2 in 1:length(l2)){
  x_2 = l2[i2]
  for(k in 1:length(h)){
    sd = h[k]
    lb_list[[k]][i2] = f(x_2, sd, sequence) 
  }
}

pdf(paste("Plot_lb_2.pdf"), height = 6, width = 8)
par(mfrow = c(3, 2))
for(i in 1:length(h)){
  plot(l2, lb_list[[i]], type = "l", ylab = "A*", xlab = "X2", main = paste("h = ",h[i]))
  }
dev.off()


###### Histogram
pdf(paste("Contour_lb_2.pdf"), height = 6, width = 8)
par(mfrow = c(3, 2))
for(i in 1:length(h)){
  hist(lb_list[[i]], bins = 100,  xlab = "X2", ylab = "density", main = paste("h = ",h[i]))
}
dev.off()

for(i in 1:length(h)){
  df = data.frame(lb_list[[i]])
  ggplot(df, aes(x = l1, y = l2)) +
    geom_density_2d_filled() +
    xlab("X1") +
    ylab("X2") +
    ggtitle(paste("h = ", h[i]))
  ggsave(paste("Contour_lb_",i,".pdf"))
}

k= 4
lowb_lb(l1[99], l2[99], h[k])



##################################################
################################################## CMH -updated (A1)
sequence = 1e2
h = c(0.001, 0.01, 1, 10, 50, 100) # proposal standard deviation

l1 = seq(-5, 5, length.out = sequence)
l2 = seq(-5, 5, length.out = sequence)
lb_list = list()

for(i in 1:length(h)){
  lb_list[[i]] = matrix(0, nrow = length(l1), length(l2))
}


for(i1 in 1:length(l1)){
for(i2 in 1:length(l2)){
  for(k in 1:length(h)){
    x_2 = l2[i2]
    sd = h[k]
    lb_list[[k]][i1, i2] = integral(f, bounds = list(x = c(-5, 5), y = c(-5, 5)))$value
  }
}
}












