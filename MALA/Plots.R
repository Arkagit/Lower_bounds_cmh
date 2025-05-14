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
pdf(paste("As_MALA.pdf"), height = 6, width = 18)
par(mfrow = c(1,3))
for(i in 1:length(h)){
  plot(l2, 1 - lb_list[[i]], type = "l", ylim = c(0, 1),
       ylab = "(1 - A*)", xlab = TeX(r'($X_2$)'), main = paste("h = ",h[i]))
}
dev.off()

pdf(paste("Hist_al2.pdf"), height = 6, width = 8)
par(mfrow = c(1,3))
for(i in 1:length(h)){
  hist(lb_list[[i]], xlim = c(0, 1), main = paste("h = ",h[i]), xlab = "A*")
}
dev.off()


################################################# Plots - CMH2
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
levels = seq(0, 100, by = 0.01)

df = as.data.frame(matrix(0, nrow = sequence^2, ncol = 3))
df[,c(1,2)] = expand.grid(x = l1, y = l2)
df[,3] = -log((1 - c(lb_list1[[1]]))*(1 - c(lb_list2[[1]])))
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
df[,3] = -log((1 - c(lb_list1[[2]]))*(1 - c(lb_list2[[2]])))
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
df[,3] = -log((1 - c(lb_list1[[3]]))*(1 - c(lb_list2[[3]])))
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


