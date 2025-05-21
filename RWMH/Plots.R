set.seed(1234)

library(ggplot2)
library(dplyr)
library(cowplot)
library(foreach)
library(doParallel)
library(latex2exp)
library(patchwork)

load("lower_bound.Rdata")

################################################## Plots - CMH1
pdf("As_RWMH.pdf", height = 8, width = 25)

par(mfrow = c(1, 3), mar = c(7, 7, 6, 4)) 

for (i in 1:length(h)) {
  plot(l2, 1 - lb_list[[i]],
       type = "l",
       ylim = c(0, 1),
       ylab = "(1 - A*)",
       xlab = TeX(r'($X_2$)'),
       main = paste("h = ", h[i]),
       cex.lab = 5,     
       cex.main = 5.0    
  )
}

dev.off()


pdf(paste("hist2.pdf"), height = 6, width = 8)
par(mfrow = c(1,3))
hist(lb_list[[1]], xlim = c(0,1))
hist(lb_list[[2]], xlim = c(0,1))
hist(lb_list[[3]], xlim = c(0,1))
dev.off()

##################################################
################################################## Plots - CMH2
for(i in 1:length(h)){
pdf(paste("A1_3d_",i,".pdf"))
  pmat = persp(l1, l2, lb_list1[[i]], xlab = "X1", ylab = "X2", zlab = "A1", theta = 45,
  main = paste("h = ",h[i]), phi = 35, shade = 0.4, cex.main = 2.2, cex.lab = 2.2)

  z.ticks <- pretty(range(lb_list1[[i]]), n = 5)

  x0 <- max(l1)
  y0 <- max(l2)

  for (z in z.ticks) {
    coords <- trans3d(x = x0, y = y0, z = z, pmat = pmat)
    points(coords, pch = 20, col = "black")
    text(coords$x, coords$y, labels = round(z, 2), pos = 4, cex = 1.1)
  }

  z.line <- trans3d(x = rep(x0, 2), y = rep(y0, 2), z = range(z.ticks), pmat = pmat)
  lines(z.line, lwd = 1.5)
  dev.off()
}


for(i in 1:length(h)){
pdf(paste("A2_3d_",i,".pdf"))
  pmat = persp(l1, l2, lb_list2[[i]], xlab = "X1", ylab = "X2", zlab = "A2", theta = 45,
  main = paste("h = ",h[i]), phi = 35, shade = 0.4, cex.main = 2.2, cex.lab = 2.2)

  z.ticks <- pretty(range(lb_list2[[i]]), n = 5)

  x0 <- max(l1)
  y0 <- max(l2)

  for (z in z.ticks) {
    coords <- trans3d(x = x0, y = y0, z = z, pmat = pmat)
    points(coords, pch = 20, col = "black")
    text(coords$x, coords$y, labels = round(z, 2), pos = 4, cex = 1.1)
  }

  z.line <- trans3d(x = rep(x0, 2), y = rep(y0, 2), z = range(z.ticks), pmat = pmat)
  lines(z.line, lwd = 1.5)
  dev.off()
}

for(i in 1:length(h)){
pdf(paste("A1A2_3d_",i,".pdf"))
  pmat = persp(l1, l2, (1-lb_list1[[i]])*(1 - lb_list2[[i]]), xlab = "X1", ylab = "X2", zlab = "(1 -A1)(1 - A2)", theta = 45,
  main = paste("h = ",h[i]), phi = 35, shade = 0.4, cex.main = 2.2, cex.lab = 2.2)

  z.ticks <- pretty(range((1-lb_list1[[i]])*(1 - lb_list2[[i]])), n = 5)

  x0 <- max(l1)
  y0 <- max(l2)

  for (z in z.ticks) {
    coords <- trans3d(x = x0, y = y0, z = z, pmat = pmat)
    points(coords, pch = 20, col = "black")
    text(coords$x, coords$y, labels = round(z, 2), pos = 4, cex = 1.1)
  }

  z.line <- trans3d(x = rep(x0, 2), y = rep(y0, 2), z = range(z.ticks), pmat = pmat)
  lines(z.line, lwd = 1.5)
  dev.off()
}

####### Plots for A1-A2
levels = seq(-50, 0, by = 0.01)

df = as.data.frame(matrix(0, nrow = sequence^2, ncol = 3))
df[,c(1,2)] = expand.grid(x = l1, y = l2)
df[,3] = log((1 - c(lb_list1[[1]]))*(1 - c(lb_list2[[1]])))
colnames(df) = c("X1", "X2", "A1")
plot1 = ggplot(df, aes(x = X1, y = X2, z = A1)) +
  geom_contour(aes(color = after_stat(level)), breaks = levels) +
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
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 40),  # X-axis label font size
    axis.title.y = element_text(size = 40),   # Y-axis label font size
    plot.title = element_text(size = 40, face = "bold", hjust = 0.5)
  )

df = as.data.frame(matrix(0, nrow = sequence^2, ncol = 3))
df[,c(1,2)] = expand.grid(x = l1, y = l2)
df[,3] = log((1 - c(lb_list1[[2]]))*(1 - c(lb_list2[[2]])))
colnames(df) = c("X1", "X2", "A1")
  plot2 = ggplot(df, aes(x = X1, y = X2, z = A1)) +
  geom_contour(aes(color = after_stat(level)), breaks = levels) +
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
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 40),  # X-axis label font size
    axis.title.y = element_text(size = 40),   # Y-axis label font size
    plot.title = element_text(size = 40, face = "bold", hjust = 0.5)
  )

df = as.data.frame(matrix(0, nrow = sequence^2, ncol = 3))
df[,c(1,2)] = expand.grid(x = l1, y = l2)
df[,3] = log((1 - c(lb_list1[[3]]))*(1 - c(lb_list2[[3]])))
colnames(df) = c("X1", "X2", "A1")
  plot3 = ggplot(df, aes(x = X1, y = X2, z = A1)) +
  geom_contour(aes(color = after_stat(level)), breaks = levels) +
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
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 40),  # X-axis label font size
    axis.title.y = element_text(size = 40),   # Y-axis label font size
    plot.title = element_text(size = 40, face = "bold", hjust = 0.5)
  )

combined_plot = (plot1/plot2/plot3)
ggsave(paste("A12_RWMH.pdf"), combined_plot, height = 20, width = 8, units = "in")

######################