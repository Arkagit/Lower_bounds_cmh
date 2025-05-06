set.seed(1234)

library(ggplot2)
library(dplyr)
library(cowplot)
library(foreach)
library(doParallel)
library(latex2exp)
library(patchwork)
################################################## Function
f1 = function(x_1, x_2, sd, r1){
  #foo1 = rep(0, length(r1))
  foo1 = pmin(exp(dnorm(r1, 0, sqrt(1/(2*(1 + x_2^2))), log = TRUE) - dnorm(x_1, 0, sqrt(1/(2*(1 + x_2^2))), log = TRUE) +
       dnorm(x_1, r1 - r1*sd^2 - r1*x_2^2*sd^2, sd, log = TRUE) - dnorm(r1, x_1 - x_1*sd^2 - x_1*x_2^2*sd^2, sd, log = TRUE)), 1)
  dum1 = mean(foo1)
	return(dum1)
}

f2 = function(x_1, x_2, sd, r2){
  foo2 = pmin(exp(dnorm(r2, 0, sqrt(1/(2*(1 + x_1^2))), log = TRUE) - dnorm(x_2, 0, sqrt(1/(2*(1 + x_1^2))), log = TRUE) +
       dnorm(x_2, r2 - r2*sd^2 - r2*x_1^2*sd^2, sd, log = TRUE) - dnorm(r2, x_2 - x_2*sd^2 - x_2*x_1^2*sd^2, sd, log = TRUE)), 1)
 	dum2 = mean(foo2)
	return(dum2)
}

##################################################
################################################## CMH -updated (A1)
sequence = 1e2
mc.size = 1e4
h = c(0.001, 1, 5) # proposal standard deviation

l1 = seq(-5, 5, length.out = sequence)
l2 = seq(-5, 5, length.out = sequence)

lb_list1 = list()
lb_list2 = list()

for(i in 1:length(h)){
  lb_list1[[i]] = matrix(0, nrow = length(l1), ncol = length(l2))
  lb_list2[[i]] = matrix(0, nrow = length(l1), ncol = length(l2))
}

parallel::detectCores()
n.cores <- 4
doParallel::registerDoParallel(cores = n.cores)


lb_list1 = foreach(k = 1:length(h))%dopar%{
  print(k)
  sd = h[k]
  dumat = matrix(0, nrow = length(l1), ncol = length(l2))
	for(i2 in 1:length(l2)){
	  for(i1 in 1:length(l1)){
	  	x_1 = l1[i1]
	    x_2 = l2[i2]
	    r1 = rnorm(mc.size, x_1 - x_1*sd^2 - x_1*x_2^2*sd^2, sd)
	    dumat[i1, i2] = f1(x_1, x_2, sd, r1)
	  }
	}
  dumat
}

lb_list2 = foreach(k = 1:length(h))%dopar%{
  print(k)
  sd = h[k]
  dumat = matrix(0, nrow = length(l1), ncol = length(l2))
  for(i2 in 1:length(l2)){
    for(i1 in 1:length(l1)){
      x_1 = l1[i1]
      x_2 = l2[i2]
      r2 = rnorm(mc.size, x_2 - x_2*sd^2 - x_2*x_1^2*sd^2, sd)
      dumat[i1, i2] = f2(x_1, x_2, sd, r2)
    }
  }
  dumat
}


for(i in 1:length(h)){
pdf(paste("A1_3d_",i,".pdf"))
  persp(l1, l2, lb_list1[[i]], xlab = "X1", ylab = "X2", zlab = "A1", theta = 45,
  main = paste("h = ",h[i]), phi = 35, shade = 0.4)
  dev.off()
}


for(i in 1:length(h)){
pdf(paste("A2_3d_",i,".pdf"))
  persp(l1, l2, lb_list2[[i]], xlab = "X1", ylab = "X2", zlab = "A2", theta = 45,
  main = paste("h = ",h[i]), phi = 35, shade = 0.4)
  dev.off()
}


####### Plots for A1
levels = seq(-100, 0, by = 0.01)

df = as.data.frame(matrix(0, nrow = sequence^2, ncol = 3))
df[,c(1,2)] = expand.grid(x = l1, y = l2)
df[,3] = log((1 - c(lb_list1[[1]]))*(1 - c(lb_list2[[1]])))
colnames(df) = c("X1", "X2", "A1")
plot1 = ggplot(df, aes(x = X1, y = X2, z = A1)) +
  geom_contour(aes(color = ..level..), breaks = levels) +
  scale_color_viridis_c(
    limits = range(levels),   
    #breaks = levels,          
    guide = guide_colorbar(title = "Contour levels")
  ) +
  labs(
    title = paste("h = ", h[1]),
    x = TeX(r'($X_1$)'),
    y = TeX(r'($X_2$)'),
    color = "Elevation"
  ) +
  theme_minimal()

df = as.data.frame(matrix(0, nrow = sequence^2, ncol = 3))
df[,c(1,2)] = expand.grid(x = l1, y = l2)
df[,3] = log((1 - c(lb_list1[[2]]))*(1 - c(lb_list2[[2]])))
colnames(df) = c("X1", "X2", "A1")
  plot2 = ggplot(df, aes(x = X1, y = X2, z = A1)) +
  geom_contour(aes(color = ..level..), breaks = levels) +
  scale_color_viridis_c(
    limits = range(levels),   
    #breaks = levels,          
    guide = guide_colorbar(title = "Contour levels")
  ) +
  labs(
    title = paste("h = ", h[2]),
    x = TeX(r'($X_1$)'),
    y = TeX(r'($X_2$)'),
    color = "Elevation"
  ) +
  theme_minimal()

df = as.data.frame(matrix(0, nrow = sequence^2, ncol = 3))
df[,c(1,2)] = expand.grid(x = l1, y = l2)
df[,3] = log((1 - c(lb_list1[[3]]))*(1 - c(lb_list2[[3]])))
colnames(df) = c("X1", "X2", "A1")
  plot3 = ggplot(df, aes(x = X1, y = X2, z = A1)) +
  geom_contour(aes(color = ..level..), breaks = levels) +
  scale_color_viridis_c(
    limits = range(levels),   
    #breaks = levels,          
    guide = guide_colorbar(title = "Contour levels")
  ) +
  labs(
    title = paste("h = ", h[3]),
    x = TeX(r'($X_1$)'),
    y = TeX(r'($X_2$)'),
    color = "Elevation"
  ) +
  theme_minimal()

combined_plot = (plot1/plot2/plot3)
ggsave(paste("A12_MALA.pdf"), combined_plot, height = 18, width = 8, units = "in")




####### Plots for A1
df = cbind(expand.grid(x = l1, y = l2), c(lb_list1[[1]]))
colnames(df) = c("X1", "X2", "A1")
plot1 = ggplot(df, aes(x = X1, y = X2, z = A1)) +
geom_contour_filled(bins = 4) +
xlab(TeX(r'($X_1$)')) +
ylab(TeX(r'($X_2$)'))+
ggtitle(paste("h = ", h[1]))

df <- cbind(expand.grid(x = l1, y = l2), c(lb_list1[[2]]))
colnames(df) = c("X1", "X2", "A1")
  plot2 = ggplot(df, aes(x = X1, y = X2, z = A1)) +
    geom_contour_filled(bins = 4) +
    xlab(TeX(r'($X_1$)')) +
    ylab(TeX(r'($X_2$)'))+
    ggtitle(paste("h = ", h[2]))

df <- cbind(expand.grid(x = l1, y = l2), c(lb_list1[[3]]))
colnames(df) = c("X1", "X2", "A1")
  plot3 = ggplot(df, aes(x = X1, y = X2, z = A1)) +
    xlab(TeX(r'($X_1$)')) +
    ylab(TeX(r'($X_2$)'))+
    geom_contour_filled(bins = 4) +
    ggtitle(paste("h = ", h[3]))

combined_plot <- plot_grid(plot1, plot2, plot3, ncol = 1, nrow = 3)
combined_plot
ggsave(paste("A1_MALA.pdf"), width = 8, height = 18, units = "in")

####### Plots for A2
df = cbind(expand.grid(x = l1, y = l2), c(lb_list2[[1]]))
colnames(df) = c("X1", "X2", "A2")
plot1 = ggplot(df, aes(x = X1, y = X2, z = A2)) +
geom_contour_filled(bins = 4) +
xlab(TeX(r'($X_1$)')) +
ylab(TeX(r'($X_2$)'))+
ggtitle(paste("h = ", h[1]))

df <- cbind(expand.grid(x = l1, y = l2), c(lb_list2[[2]]))
colnames(df) = c("X1", "X2", "A2")
  plot2 = ggplot(df, aes(x = X1, y = X2, z = A2)) +
    geom_contour_filled(bins = 4) +
    xlab(TeX(r'($X_1$)')) +
    ylab(TeX(r'($X_2$)'))+
    ggtitle(paste("h = ", h[2]))

df <- cbind(expand.grid(x = l1, y = l2), c(lb_list2[[3]]))
colnames(df) = c("X1", "X2", "A2")
  plot3 = ggplot(df, aes(x = X1, y = X2, z = A2)) +
    geom_contour_filled(bins = 4) +
    xlab(TeX(r'($X_1$)')) +
    ylab(TeX(r'($X_2$)'))+
    ggtitle(paste("h = ", h[3]))

combined_plot <- plot_grid(plot1, plot2, plot3, ncol = 1, nrow = 3)
combined_plot
ggsave(paste("A2_MALA.pdf"), width = 8, height = 18, units = "in")



