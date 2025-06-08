############################################################

load("RAM_data.Rdata")

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

RAM_d = list(matrix(0, nrow = nsim, ncol = repet),matrix(0, nrow = nsim, ncol = repet),
			matrix(0, nrow = nsim, ncol = repet),matrix(0, nrow = nsim, ncol = repet))
RAM_test = list(matrix(0, nrow = nsim, ncol = repet),matrix(0, nrow = nsim, ncol = repet),
			matrix(0, nrow = nsim, ncol = repet),matrix(0, nrow = nsim, ncol = repet))
Metro_d = list(matrix(0, nrow = nsim, ncol = repet),matrix(0, nrow = nsim, ncol = repet),
			matrix(0, nrow = nsim, ncol = repet),matrix(0, nrow = nsim, ncol = repet))
RAM_J_Barker_d = list(matrix(0, nrow = nsim, ncol = repet),matrix(0, nrow = nsim, ncol = repet),
			matrix(0, nrow = nsim, ncol = repet),matrix(0, nrow = nsim, ncol = repet))


for(k in 1:4){
	for(i in 1:nsim){
		for(j in 1:repet){
			RAM_d[[k]][i, j] = alpha_list[[i]][[2]][[j]][,k]
            Metro_d[[k]][i, j] = alpha_list[[i]][[4]][[j]][,k]
            RAM_J_Barker_d[[k]][i, j] = alpha_list[[i]][[1]][[j]][,k]
            RAM_test[[k]][i, j] = alpha_list[[i]][[3]][[j]][,k]
	    }
	}
}


RAM = (1 - rowMeans(RAM_d[[1]]))*(1 - rowMeans(RAM_d[[2]]))*
				(1 - rowMeans(RAM_d[[3]]))*(1 - rowMeans(RAM_d[[4]]))

RAM_test = (1 - rowMeans(RAM_test[[1]]))*(1 - rowMeans(RAM_test[[2]]))*
				(1 - rowMeans(RAM_test[[3]]))*(1 - rowMeans(RAM_test[[4]]))
	
Metro = (1 - rowMeans(Metro_d[[1]]))*(1 - rowMeans(Metro_d[[2]]))*
				(1 - rowMeans(Metro_d[[3]]))*(1 - rowMeans(Metro_d[[4]]))

RAM_J_Barker = (1 - rowMeans(RAM_J_Barker_d[[1]]))*(1 - rowMeans(RAM_J_Barker_d[[2]]))*
				(1 - rowMeans(RAM_J_Barker_d[[3]]))*(1 - rowMeans(RAM_J_Barker_d[[4]]))


pdf(paste("lb_RAM_Metro.pdf"), height = 8, width = 8)
#par(mfrow = c(1, 2))
plot(rowMeans(jumping.scale), RAM, ylab = "Lower Bounds", xlab = "Jump Scale", type = "l", col = "blue", lwd = "3")
lines(rowMeans(jumping.scale), Metro)
lines(rowMeans(jumping.scale), RAM_J_Barker, col = "green")
lines(rowMeans(jumping.scale), RAM_test, col = "brown")
legend("bottomright", legend = c("RAM", "RAM (Barker with same chain)", "MH", "RAM (Joint Barker rate)"),
 col = c("blue", "brown", "black", "green"), lty = 1, cex = 1.5)
dev.off()


