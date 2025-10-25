set.seed(1234)

source("BVN_functions.R")

library(MASS)
library(latex2exp)
# library(foreach)
# library(doParallel)
# library(matrixcalc)

# rho values
rho = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 0.999)

## Number of iteration of A's
n = 1e4

## Auxiliary variances
lambda = 2 #c(0.1, 1, 10)

## Target mean of X
tar_mean = c(0, 0)

## Initial points
x_init = c(0, 0)
z_init = c(0, 0)

## Initializations for acceptance probs and avg acceptance probs
A_block = matrix(0, nrow = n, ncol = length(rho))
A_star_cmh = matrix(0, nrow = n, ncol = length(rho))
A_star_proximal = matrix(0, nrow = n, ncol = length(rho))

## Proposal covariance matrix
h = c(0.1, 1, 10)

#prop_cov = matrix(c(2, 0, 0, 2), nrow = 2)

output = list()

out_rho = list()

for(i1 in 1:length(h)){

	prop_cov = matrix(c(h[i1], 0, 0, h[i1]), nrow = 2)

	for(i in 1:length(rho)){

		tar_cov = matrix(c(1, rho[i], rho[i], 1), nrow = 2)

		for(j in 1:n){

			A_block[j, i] = alpha_block(x_init, tar_mean, tar_cov, prop_cov)

			x1 = rnorm(1, x_init[1], prop_cov[1, 1])

			A_star_cmh[j, i] = alpha_cmh(x_init, x1, tar_mean, tar_cov) #mean(A_mh)

			z = rnorm(2, x_init, sqrt(lambda))

			A_star_proximal[j, i] = alpha_proximal(x_init, z, tar_mean,
			                                tar_cov, prop_cov, lambda)

		}

	}

	output[[i1]] = list(A_block, A_star_cmh, A_star_proximal)
}


lb_block = matrix(0, nrow = length(rho), ncol = length(h))
lb_cmh = matrix(0, nrow = length(rho), ncol = length(h))
lb_proximal = matrix(0, nrow = length(rho), ncol = length(h))


for(i in 1:length(h))
{
	lb_block[, i] = 1 - colMeans(output[[i]][[1]])
	lb_cmh[, i] = 1 - colMeans(output[[i]][[2]])
	lb_proximal[, i] = 1 - colMeans(output[[i]][[3]])
}


pdf("Gaussian_Lower_bound.pdf", height = 10, width = 20)

par(mfrow = c(1, 3))

for(i in 1:length(h)){

	par(mar = c(5.1, 4.8, 4.1, 2.1))
	plot(rho, lb_block[, i], type = 'b', col = "black", ylim = c(0,1),
	      ylab = "Estimated Lower Bounds", xlab = expression(rho), main = paste("h = ", h[i]),
	      cex.main = 2, cex.lab = 2)
	lines(rho, lb_proximal[, i], type = 'b', col = "blue")
	lines(rho, lb_cmh[, i], type = 'b', col = "brown")
	legend("bottomright", bty = "n", legend = c("Complete Block MH", "CMH1",
	 expression("Augmented CMH ("* lambda* "= 2)")), col = c("black", "brown", "blue"), lty = 1, cex = 1.6)

}

dev.off()


