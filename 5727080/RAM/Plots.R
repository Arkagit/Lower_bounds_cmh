############################################################

load("RAM_data.Rdata")

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

RAM = matrix(0, nrow = nsim, ncol = 1)
Metro = matrix(0, nrow = nsim, ncol = 1)

for(i in 1:nsim){
  RAM[i, ] = prod(1 - colMeans(alpha_list[[i]][[1]]))
  Metro[i, ] = prod(1 - colMeans(alpha_list[[i]][[2]]))
}

pdf(paste("lb_hist.pdf"), height = 6, width = 16)
par(mfrow = c(1, 2))
hist(RAM, xlab = paste("Lower Bound"), main  = "RAM", breaks = 20, ylim = c(0, 20))
hist(Metro, xlab = paste("Lower Bound"), main  = "MH", breaks = 20, ylim = c(0, 20))
dev.off()
