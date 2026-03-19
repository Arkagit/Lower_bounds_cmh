set.seed(1234)

lb_TV = function(x, y){
  b = 1 - (3 - 4*exp(-1))/(1260* y^4 * (1-y)^4 * ((1 - exp(-1)) + y*(1 - 2*exp(-1))))
  if(b < 0){
    b = 0
  }
  return(b)
}

ub_TV = function(x, y){
  b1 = 1 - 2*(3 - 4*exp(-1))/(5*(6 - 9*exp(-1)))
  return(b1)
}

y_grid = seq(1, 1e4, 1)/1e4
x_grid = seq(1, 1e4, 1)/1e4

lb_grid = rep(0, 1e4)
ub_grid = rep(0, 1e4)

for(i in 1: 1e4){
  lb_grid[i] = lb_TV(x_grid[i], y_grid[i])
  ub_grid[i] = ub_TV(x_grid[i], y_grid[i])
}
sum(lb_grid > ub_grid)

pdf("Lower_Upper_bound_2.pdf", width = 5, height = 5)

plot(y_grid, lb_grid, xlab = "Y", ylab = "Total Variation Distance",
 ylim = c(0, 1), lwd = 1, type = "l", col = "blue")
lines(y_grid, ub_grid, col = "brown", lwd = 1)

legend("topright", bty = "n", legend = c("Upper Bound", "Lower Bound"),
 col = c("brown", "blue"), lty = 1, cex = 1)

dev.off()