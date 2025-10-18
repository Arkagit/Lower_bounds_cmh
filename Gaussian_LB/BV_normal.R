set.seed(1234)

source("BVN_functions.R")

library(MASS)
library(latex2exp)

rho = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 0.999)
n = 1e4
m = 1e2
lambda = c(1, 10)
tar_mean = c(0, 0)
x_init = c(0, 0)
z_init = c(0, 0)
A_B = matrix(0, nrow = n, ncol = length(rho))
A_mh = matrix(0, nrow = n, ncol = length(rho))
A_p = list(matrix(0, nrow = n, ncol = length(rho)), matrix(0, nrow = n, ncol = length(rho)))
prop_cov = matrix(c(2, 0, 0, 2), nrow = 2)

output = list()


for(i in 1:length(rho)){

	tar_cov = matrix(c(1, rho[i], rho[i], 1), nrow = 2)

	for(j in 1:n){

		A_B[j, i] = alpha_B(x_init, tar_mean, tar_cov, prop_cov)

		A_mh[j, i] = alpha_mh(x_init, tar_mean, tar_cov, prop_cov)

		for(k in 1:length(lambda)){

			A_p[[k]][j, i] = alpha_p(x_init, tar_mean, tar_cov, prop_cov, lambda[k])
		}
	}
}


pdf("Gaussian_Lower_bound.pdf", height = 6, width = 8)
par(mar = c(5.1, 4.8, 4.1, 2.1))
plot(rho, 1 - colMeans(A_B), type = 'b', col = "black", ylim = c(0,1), ylab = "Estimated Lower Bounds", xlab = expression(rho))
lines(rho, 1 - colMeans(A_p[[1]]), type = 'b', col = "blue")
lines(rho, 1 - colMeans(A_p[[2]]), type = 'b', col = "brown")
lines(rho, 1 - colMeans(A_mh), type = 'b', col = "red")
legend("bottomright", bty = "n", legend = c("Complete Block MH", "CMH1", expression("Augmented CMH ("* lambda* "= 1)"), expression("Augmented CMH ("* lambda* "= 10)")), 
	col = c("black", "red", "blue", "brown"), lty = 1, cex=0.65)
dev.off()


# #### Block updates and CMH1 updates

# for(i in 1:length(rho)){

# 	tar_cov = matrix(c(1, rho[i], rho[i], 1), nrow = 2)
# 	prop_cov = matrix(c(2, 0, 0, 2), nrow = 2)

# 	# proposal pdf for block MH
# 	q = function(y1, x1){
# 		out1 = exp(-log(2* pi) - log(det(prop_cov))/2 - t(y1 - x1)%*%solve(prop_cov)%*%(y1 - x1)/2)
# 		return(out1)
# 	}

# 	# target pdf for block MH
# 	f = function(x2){
# 		out2 = exp(-log(2* pi) - log(det(tar_cov))/2 - t(x2)%*%solve(tar_cov)%*%(x2)/2)
# 		return(out2)
# 	}


# 	# A_B simulation
# 	x = c(0, 0)
# 	A_B = 0

# 	for(j in 1:n){

# 		y = c(rnorm(1, x[1], sqrt(prop_cov[1,1])), rnorm(1, x[2], sqrt(prop_cov[2,2])))

# 		r = exp(log(f(y)) + log(q(x, y)) - log(f(x)) - log(q(y, x)))

# 		if(r > 1){
# 			alpha = 1
# 		}else{  
# 			alpha = r
# 		}

# 		A_B = A_B + alpha/n
# 	}
# 	print(1)

# 	#A_B simulation
# 	A_2 = 0

# 	for(j in 1:n){
# 	y = rnorm(1, x[2], sqrt(prop_cov[2,2]))

# 	r = exp(log(dnorm(y, x[2] + rho[i]*x[1]/tar_cov[1, 1], sqrt(tar_cov[2,2] - (rho[i]^(2))/tar_cov[1,1]))) + 
# 	    log(dnorm(x[2], y, sqrt(prop_cov[2,2]))) -
# 	 log(dnorm(x[2], x[2] + rho[i]*x[1]/tar_cov[1,1], sqrt(tar_cov[2,2] - (rho[i]^(2))/tar_cov[1,1]))) - 
# 	 log(dnorm(y, x[2], sqrt(prop_cov[2,2]))))

# 	if(r > 1){
# 		alpha = 1
# 	}else{
# 		alpha = r
# 	}

# 	A_2 = A_2 + alpha/n
# 	}

# 	print(2)

# 	####################### Proximal Aug CARB1
# 	x_init = c(0, 0)

# 	z_init = c(0, 0)

# 	A_aug = rep(0, length(lambda))

# 	for(k in 1:length(lambda)){

# 		Omega = solve(solve(tar_cov) + diag(c(1, 1))/lambda[k])

# 		z = matrix(0, nrow = n, ncol = 2)

# 		U = matrix(0, nrow = n, ncol = 2)

# 		r_mh = rep(0, n)

# 		for(j in 1:n){

# 			z[j, ] = mvrnorm(1, x_init, lambda[k]*diag(c(1, 1)))

# 			pop_mean = Omega%*%z[j, ]/lambda[k]

# 			dummy = mvrnorm(1, x_init, 2*diag(c(1, 1)))

# 			r_mh[j] = exp(- 0.5*t(dummy - pop_mean)%*%solve(Omega)%*%(dummy - pop_mean) +
# 				0.5*t(x_init - x_init)%*%solve(Omega)%*%(x_init - x_init) - 
# 				0.5*t(x_init - pop_mean)%*%solve(2*diag(c(1, 1)))%*%(x_init - pop_mean) +
# 				0.5*t(dummy - x_init)%*%solve(2*diag(c(1, 1)))%*%(dummy - x_init))
# 		}

# 		A_aug[k] = mean(r_mh)
# }

# 	output[[i]] = c(A_B, A_2, A_aug)

# }

# A_B = rep(0, length(rho))
# A_2 = rep(0, length(rho))
# A_aug = matrix(0, nrow = length(lambda), ncol = length(rho))

# for(i in 1:length(rho)){
# 	A_B[i] = output[[i]][1]
# 	A_2[i] = output[[i]][2]
# 	for(j in 1: length(lambda)){
# 		A_aug[j, i] = output[[i]][j+2]
# 	}
# }


# ####################### Proximal Aug CARB1

# # h= c(1, 10)
# # x_init = c(0, 0)

# # z_init = c(0, 0)

# # A_aug = matrix(0, nrow = length(rho), ncol = length(h))

# # for(k in 1:length(h)){

# # 	for(j in 1:length(rho)){

# # 		tar_cov = matrix(c(1, rho[j], rho[j], 1), nrow = 2)

# # 		z = matrix(0, nrow = n, ncol = 2)

# # 		U = matrix(0, nrow = n, ncol = 2)

# # 		r_mh = rep(0, n)

# # 		for(i in 1:n){
# # 			z[i,] = mvrnorm(1, x_init, h[k]*diag(c(1, 1)))

# # 			Omega = solve(solve(tar_cov) + diag(c(1, 1))/h[k])

# # 			pop_mean = Omega%*%z[i,]/h[k]

# # 			dummy = mvrnorm(1, pop_mean, 2*diag(c(1, 1)))

# # 			r_mh[i] = exp(- 0.5*t(dummy - pop_mean)%*%solve(Omega)%*%(dummy - pop_mean) +
# # 				0.5*t(c(0, 0) - x_init)%*%solve(Omega)%*%(c(0, 0) - x_init) - 
# # 				0.5*t(x_init - pop_mean)%*%solve(2*diag(c(1, 1)))%*%(x_init - pop_mean) +
# # 				0.5*t(dummy - x_init)%*%solve(2*diag(c(1, 1)))%*%(dummy - x_init))
# # 		}

# # 		A_aug[j, k] = mean(r_mh)

# # 		print(j)
# # 	}
# # }



# pdf("Gaussian_Lower_bound.pdf", height = 6, width = 8)
# par(mar = c(5.1, 4.8, 4.1, 2.1))
# plot(rho, 1 - A_B, type = 'b', col = "black", ylim = c(0,1), ylab = "Estimated Lower Bounds")
# lines(rho, 1 - A_aug[1, ], type = 'b', col = "blue")
# lines(rho, 1 - A_aug[2, ], type = 'b', col = "brown")
# lines(rho, 1 - A_2, type = 'b', col = "red")
# legend("bottomright", bty = "n", legend = c("Complete Block", "CMH1", "CMH (Augmented) (h = 1)", "CMH (Augmented) (h = 10)"), 
# 	col = c("black", "red", "blue", "brown"), lty = 1, cex=0.65)
# dev.off()
