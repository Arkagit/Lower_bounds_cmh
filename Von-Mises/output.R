set.seed(1234)

source("functions.R")

c0 = 0.68
a = 0.1
b = 1
h_choices = 4

dd_psi(100)

n = c(5e2, 7e2, 1e3, 3e3, 5e3, 1e4, 5e4, 1e5)
h = matrix(0, nrow = h_choices, ncol = length(n))
lb_out = matrix(0, nrow = h_choices, ncol = length(n))

h[1, ] = 5
h[2, ] = 5/2
h[3, ] = 5/log(n)
h[4, ] = 5/n

for(i in 1:4){
    lb_out[i,] = vm_lb(n, h[i,], a, b)
}

lb_min = min(lb_out)
lb_max = max(lb_out)

pdf("VM_lower_bound.pdf", height = 8, width = 9)

# Increase bottom margin for x-axis labels + legend
par(mar = c(8, 6, 5, 4))

x <- seq_along(n)

plot(x, lb_out[1, ],
     type = "b",
     xaxt = "n",
     pch = 16,
     col = "steelblue",
     lwd = 2.5,
     cex = 1.1,
     ylim = c(lb_min, lb_max),
     xlab = "Data Size",
     ylab = "Lower Bounds",
     cex.lab = 1.5,
     cex.axis = 1.3)

axis(1,
     at = x,
     labels = c("500", "700", "1K", "3K",
                "5K", "10K", "50K", "100K"),
     cex.axis = 1.3)

lines(x, lb_out[2, ],
      type = "b",
      pch = 16,
      col = "darkolivegreen",
      lwd = 2.5,
      cex = 1.1)

lines(x, lb_out[3, ],
      type = "b",
      pch = 16,
      col = "purple",
      lwd = 2.5,
      cex = 1.1)

lines(x, lb_out[4, ],
      type = "b",
      pch = 16,
      col = "brown",
      lwd = 2.5,
      cex = 1.1)

#grid(lty = 2, col = "gray")

legend("right",
       legend = c(
         expression(h == 5),
         expression(h == 5/2),
         expression(h == 5/log(n)),
         expression(h == 5/n)
       ),
       col = c("steelblue", "darkolivegreen", "purple", "brown"),
       pch = 16,
       lty = 1,
       lwd = 2.5,
       bty = "n",
       cex = 1.1)

dev.off()
