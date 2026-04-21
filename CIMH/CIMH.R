set.seed(1000)

library(latex2exp)
library(dplyr)

x_1 <- seq(0, 3, length.out = 200)
x_2 <- seq(0, 6, length.out = 200)


bv_pdf = function(x_1, x_2){
  if (!is.numeric(x_1) || !is.numeric(x_2)) {
    stop("Inputs x_1 and x_2 must be numeric.")
  }

  out <- numeric(length(x_1))

  condition <- (x_1 > 0) & (x_2 > x_1)

  out[condition] <- exp(- x_2[condition] + x_1[condition])
  
  return(out)
}

grid_points <- expand.grid(x = x_1, y = x_2)

z_values = bv_pdf(grid_points[,1], grid_points[,2])

z_matrix <- matrix(z_values, nrow = length(x_1), ncol = length(x_2), byrow = FALSE)

pdf("CIMH_contour_plot.pdf", width = 4, height = 4)

contour(x = x_1,
        y = x_2,
        z = z_matrix,
        col = "black",
        lwd = 1,
        nlevels = 10,
        xlab = TeX(r'($X$)'),
        ylab = TeX(r'($Y$)'),
        cex = 2) 

dev.off()


###########################################################
x_2 = seq(1, 6, length.out = 200)

lower_bound = function(x_2){
	out = 1 - exp(-x_2^2 + x_2)
	return(out)
}

pdf("CIMH_lower_bound.pdf", width = 4, height = 4)

plot(x_2, lower_bound(x_2), xlab = TeX(r'($Y$)'), 
  ylab = "Lower Bound", type = "l", cex = 2)

dev.off()





