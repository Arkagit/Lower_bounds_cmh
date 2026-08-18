set.seed(1234)

source("functions.R")

# Example
n <- c(1e3, 3e3, 5e3, 7e3, 1e4, 3e4, 5e4)
d <- c(2, 4, 6, 8, 10, 12, 14)
zeta_star = 1
b = 0.001

num_h = 4




lb1 = matrix(0, ncol = num_h, nrow = length(d))
lb2 = matrix(0, ncol = num_h, nrow = length(d))


for(j in 1:length(d)){
	xi_star = rnorm(d[j], mean = 0, sd = 1) 
	h = matrix(0, nrow = num_h, ncol = length(d[j]))
	h[1,] = 2
	h[2,] = 2/d[j]
	h[3,] = 2/(d[j]* n[j])
	h[4,] = 2/(d[j]* log(n[j]))

	W1 = design_mat_poor_cond(n[j], d[j])$W
	W2 = design_mat_good_cond(n[j], d[j])$W
	

	for(i in 1:num_h){
		lb1[j, i] = VM_lower_bound(n[j], d[j], h[i], b, zeta_star, xi_star, W1)
		lb2[j, i] = VM_lower_bound(n[j], d[j], h[i], b, zeta_star, xi_star, W2)
	}
}

lb1
lb2


num_pairs <- length(d)  # Number of (d,n) pairs

# Create labels for each (d, n) pair
pair_labels <- paste0("(", d, ",", format(n, scientific = FALSE), ")")


# Create PDF with wider dimensions for better label spacing
pdf("BLR_lower_bound_logn.pdf", width = 14, height = 7)

# Increase bottom margin significantly to fit horizontal labels
par(mar = c(8, 7, 5, 4), mgp = c(4, 1, 0))

# X-axis positions (equally spaced)
x_pos <- 1:num_pairs

# Get y-axis limits
y_min <- min(lb1, na.rm = TRUE)
y_max <- max(lb1, na.rm = TRUE)

# Create base plot with first h value (blue line)
plot(x_pos, lb1[, 1],
     type = "b",
     pch = 16,
     col = "blue",
     lwd = 2.5,
     cex = 1.2,
     ylim = c(y_min, y_max),
     xaxt = "n",
     xlab = "",
     ylab = "Lower Bounds",
     cex.lab = 2,
     cex.axis = 1.2,
     cex.main = 1.3)

# Add the remaining 3 lines with different colors
lines(x_pos, lb1[, 2],
      type = "b",
      pch = 17,
      col = "darkolivegreen",
      lwd = 2.5,
      cex = 1.2)

lines(x_pos, lb1[, 3],
      type = "b",
      pch = 15,
      col = "brown",
      lwd = 2.5,
      cex = 1.2)

lines(x_pos, lb1[, 4],
      type = "b",
      pch = 18,
      col = "purple",
      lwd = 2.5,
      cex = 1.2)

# Add X-axis with (d, n) pair labels - HORIZONTAL
axis(1,
     at = x_pos,
     labels = pair_labels,
     cex.axis = 1.5,
     las = 1)  # Keep labels horizontal (las = 1)

# Add X-axis title below the labels
mtext("(Dimension, Data Size)", side = 1, line = 4.5, cex = 2)


legend("bottomright", bty = "n",
       legend = c(
         expression(h == 1),
         expression(h == 1/d),
         expression(h == 1/d*n),
         expression(h == 1/d*log(n))
       ),
       col = c("blue", "darkolivegreen", "brown", "purple"),
       pch = c(16, 17, 15, 18),
       lty = 1,
       lwd = 2.5,
       cex = 2,
       bg = "white")

dev.off()




##############################################################

pdf("BLR_lower_bound_n.pdf", width = 14, height = 7)

# Increase bottom margin significantly to fit horizontal labels
par(mar = c(8, 7, 5, 4), mgp = c(4, 1, 0))

# X-axis positions (equally spaced)
x_pos <- 1:num_pairs

# Get y-axis limits
y_min <- min(lb2, na.rm = TRUE)
y_max <- max(lb2, na.rm = TRUE)

# Create base plot with first h value (blue line)
plot(x_pos, lb2[, 1],
     type = "b",
     pch = 16,
     col = "blue",
     lwd = 2.5,
     cex = 1.2,
     ylim = c(y_min, y_max),
     xaxt = "n",
     xlab = "",
     ylab = "Lower Bounds",
     cex.lab = 2,
     cex.axis = 1.5,
     cex.main = 1.3)

# Add the remaining 3 lines with different colors
lines(x_pos, lb2[, 2],
      type = "b",
      pch = 17,
      col = "darkolivegreen",
      lwd = 2.5,
      cex = 1.2)

lines(x_pos, lb2[, 3],
      type = "b",
      pch = 15,
      col = "brown",
      lwd = 2.5,
      cex = 1.2)

lines(x_pos, lb2[, 4],
      type = "b",
      pch = 18,
      col = "purple",
      lwd = 2.5,
      cex = 1.2)


axis(1,
     at = x_pos,
     labels = pair_labels,
     cex.axis = 1.5,
     las = 1)  

mtext("(Dimension, Data Size)", side = 1, line = 4.5, cex = 2)


legend("bottomright", bty = "n",
       legend = c(
         expression(h == 1),
         expression(h == 1/d),
         expression(h == 1/d*n),
         expression(h == 1/d*log(n))
       ),
       col = c("blue", "darkolivegreen", "brown", "purple"),
       pch = c(16, 17, 15, 18),
       lty = 1,
       lwd = 2.5,
       cex = 2,
       bg = "white")

dev.off()