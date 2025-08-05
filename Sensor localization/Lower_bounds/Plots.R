############################################################

load("RAM_data.Rdata")

RAM_d = list(matrix(0, nrow = nsim, ncol = repet),matrix(0, nrow = nsim, ncol = repet),
			matrix(0, nrow = nsim, ncol = repet),matrix(0, nrow = nsim, ncol = repet))
RAM_Barker_d = list(matrix(0, nrow = nsim, ncol = repet),matrix(0, nrow = nsim, ncol = repet),
			matrix(0, nrow = nsim, ncol = repet),matrix(0, nrow = nsim, ncol = repet))
Metro_d = list(matrix(0, nrow = nsim, ncol = repet),matrix(0, nrow = nsim, ncol = repet),
			matrix(0, nrow = nsim, ncol = repet),matrix(0, nrow = nsim, ncol = repet))
RAM_2coin_d = list(matrix(0, nrow = nsim, ncol = repet),matrix(0, nrow = nsim, ncol = repet),
			matrix(0, nrow = nsim, ncol = repet),matrix(0, nrow = nsim, ncol = repet))


for(k in 1:4){
	for(i in 1:nsim){
		for(j in 1:repet){
			RAM_d[[k]][i, j] = alpha_list[[i]][[1]][[j]][,k]
            RAM_Barker_d[[k]][i, j] = alpha_list[[i]][[2]][[j]][,k]
            Metro_d[[k]][i, j] = alpha_list[[i]][[3]][[j]][,k]
            RAM_2coin_d[[k]][i, j] = alpha_list[[i]][[4]][[j]][,k]
	    }
	}
}


RAM = (1 - rowMeans(RAM_d[[1]]))*(1 - rowMeans(RAM_d[[2]]))*
				(1 - rowMeans(RAM_d[[3]]))*(1 - rowMeans(RAM_d[[4]]))

RAM_Barker = (1 - rowMeans(RAM_Barker_d[[1]]))*(1 - rowMeans(RAM_Barker_d[[2]]))*
				(1 - rowMeans(RAM_Barker_d[[3]]))*(1 - rowMeans(RAM_Barker_d[[4]]))
	
Metro = (1 - rowMeans(Metro_d[[1]]))*(1 - rowMeans(Metro_d[[2]]))*
				(1 - rowMeans(Metro_d[[3]]))*(1 - rowMeans(Metro_d[[4]]))

RAM_2coin = (1 - rowMeans(RAM_2coin_d[[1]]))*(1 - rowMeans(RAM_2coin_d[[2]]))*
				(1 - rowMeans(RAM_2coin_d[[3]]))*(1 - rowMeans(RAM_2coin_d[[4]]))


pdf(paste("lb_RAM_Metro.pdf"), height = 10, width = 10)
plot(rowMeans(jumping.scale), RAM, ylab = "Lower Bounds", xlab = "h", type = "l",
 col = "blue", lwd = 2, cex.lab = 1.3)
#lines(rowMeans(jumping.scale), RAM_Barker, col = "brown")
lines(rowMeans(jumping.scale), Metro, col = "red", lwd = 2)
lines(rowMeans(jumping.scale), RAM_2coin, col = "black", lwd = 2)
legend("bottomright", legend = c("CRAM (aux)", "CMH", "CRAB (2 coin Barker)"),
 col = c("blue", "red", "black"), lty = 1, cex = 1.5, lwd = 2)
dev.off()

log_jump = log(rowMeans(jumping.scale))
pdf(paste("Log_lb_RAM_Metro.pdf"), height = 8, width = 8)
plot(log_jump, RAM, ylab = "Lower Bounds", xlab = "log(h)", type = "l", col = "blue")
lines(log_jump, RAM_Barker, col = "brown")
lines(log_jump, Metro)
lines(log_jump, RAM_2coin, col = "purple")
legend("topleft", legend = c("RAM", "RAM (Barker)", "RWM", "RAM (2 coin Barker)"),
 col = c("blue", "brown", "black", "purple"), lty = 1, cex = 1.5)
dev.off()


