set.seed(1234)

source("BVN_functions.R")

library(MASS)
library(latex2exp)
library(foreach)
library(doParallel)
library(matrixcalc)

# rho values
rho = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 0.999)

## Number of iteration of A's
n = 1e3

## Number of iteration of \alpha's
m = 1e2

## Initialization
A_mh = rep(0, m)
A_p = rep(0, m)

## Auxiliary variances
lambda = c(1, 10)

## Target mean of X
tar_mean = c(0, 0)

## Initial points
x_init = c(0, 0)
z_init = c(0, 0)

## Initializations for acceptance probs and avg acceptance probs
A_B = matrix(0, nrow = n, ncol = length(rho))
A_s_mh = matrix(0, nrow = n, ncol = length(rho))
A_s_p = list(matrix(0, nrow = n, ncol = length(rho)), matrix(0, nrow = n, ncol = length(rho)))

## Proposal covariance matrix
prop_cov = matrix(c(2, 0, 0, 2), nrow = 2)

output = list()


for(i in 1:length(rho)){

	tar_cov = matrix(c(1, rho[i], rho[i], 1), nrow = 2)

	for(j in 1:n){

		A_B[j, i] = alpha_B(x_init, tar_mean, tar_cov, prop_cov)

		#A_mh[j, i] = alpha_mh(x_init, tar_mean, tar_cov, prop_cov)

		for(k in 1:m){
			x1 = rnorm(1, x_init[1] + tar_cov[2,1]*x_init[2]/tar_cov[2, 2], sqrt(tar_cov[1,1] - (tar_cov[2, 1]^(2))/tar_cov[2,2]))
			A_mh[k] = alpha_mh(x_init, x1, tar_mean, tar_cov, prop_cov)
		}

		A_s_mh[j, i] = mean(A_mh)

		for(k in 1:length(lambda)){

			for(l in 1:m){
				z = rnorm(2, x_init, sqrt(lambda[k]))
				A_p[k] = alpha_p(x_init, z, tar_cov, prop_cov, lambda[k])
			}

			A_s_p[[k]][j, i] = mean(A_p)
		}
		print(i)
		print(j)
	}
	list(A_B, A_s_mh, A_s_p)
}


pdf("Gaussian_Lower_bound.pdf", height = 6, width = 8)
par(mar = c(5.1, 4.8, 4.1, 2.1))
plot(rho, 1 - colMeans(A_B), type = 'b', col = "black", ylim = c(0,1), ylab = "Estimated Lower Bounds", xlab = expression(rho))
lines(rho, 1 - colMeans(A_s_p[[1]]), type = 'b', col = "blue")
lines(rho, 1 - colMeans(A_s_p[[2]]), type = 'b', col = "brown")
lines(rho, 1 - colMeans(A_s_mh), type = 'b', col = "red")
legend("bottomright", bty = "n", legend = c("Complete Block MH", "CMH1", expression("Augmented CMH ("* lambda* "= 1)"), expression("Augmented CMH ("* lambda* "= 10)")), 
	col = c("black", "red", "blue", "brown"), lty = 1, cex=0.65)
dev.off()

