set.seed(1234)

rho = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 0.999)
n = 1e5
output = list()

for(i in 1:length(rho)){

	tar_cov = matrix(c(1, rho[i], rho[i], 1), nrow = 2)
	prop_cov = matrix(c(2, 0, 0, 2), nrow = 2)

	q = function(y1, x1){
		out1 = exp(-log(2* pi) - log(det(prop_cov))/2 - t(y1 - x1)%*%solve(prop_cov)%*%(y1 - x1)/2)
		return(out1)
	}

	f = function(x2){
		out2 = exp(-log(2* pi) - log(det(tar_cov))/2 - t(x2)%*%solve(tar_cov)%*%(x2)/2)
		return(out2)
	}

	x = c(0, 0)
	A_B = 0
	A_1 = 0
	A_2 = 0



	### acceptance for block updates
	for(j in 1:n){
	y = c(rnorm(1, x[1], sqrt(prop_cov[2,2])), rnorm(1, x[2], sqrt(prop_cov[2,2])))

	r = exp(log(f(y)) + log(q(x, y)) - log(f(x)) - log(q(y, x)))

	if(r > 1){
		alpha = 1
	}else{  
		alpha = r
	}

	A_B = A_B + alpha/n
	}
	print(1)



	### acceptance for 1st compnent
	for(j in 1:n){
	y1 = rnorm(1, x[1], sqrt(prop_cov[1,1]))

	r1 = exp(log(dnorm(y1, x[1] + rho[i]*x[2]/tar_cov[1, 1], sqrt(tar_cov[1,1] - (rho[i]^(2))/tar_cov[2,2]))) + log(dnorm(x[1], y1, sqrt(prop_cov[1,1]))) -
	 log(dnorm(x[1], x[1] + rho[i]*x[2]/tar_cov[1,1], sqrt(tar_cov[1,1] - (rho[i]^(2))/tar_cov[2,2]))) - log(dnorm(y1, x[1], sqrt(prop_cov[1,1]))))

	if(r1 > 1){
		alpha1 = 1
	}else{
		alpha1 = r1
	}

	A_1 = A_1 + alpha1/n
	}
	print(2)

	### acceptance for 2nd compnent
	for(j in 1:n){
	y2 = rnorm(1, x[2], sqrt(prop_cov[2,2]))

	r2 = exp(log(dnorm(y2, x[2] + rho[i]*x[1]/tar_cov[1, 1], sqrt(tar_cov[2,2] - (rho[i]^(2))/tar_cov[1,1]))) + log(dnorm(x[2], y2, sqrt(prop_cov[2,2]))) -
	 log(dnorm(x[2], x[2] + rho[i]*x[1]/tar_cov[1,1], sqrt(tar_cov[2,2] - (rho[i]^(2))/tar_cov[1,1]))) - log(dnorm(y2, x[2], sqrt(prop_cov[2,2]))))

	if(r2 > 1){
		alpha2 = 1
	}else{
		alpha2 = r2
	}

	A_2 = A_2 + alpha2/n
	}
	print(3)

	output[[i]] = c(A_B, A_1, A_2)
}

A_B = c(); A_1 = c(); A_2 = c()

for(i in 1:length(rho)){
	A_B = c(A_B, output[[i]][1])
	A_1 = c(A_1, output[[i]][2])
	A_2 = c(A_2, output[[i]][3])
}



pdf("Gaussian_Lower_bound_CMH2.pdf", height = 6, width = 8)
par(mar = c(5.1, 4.8, 4.1, 2.1))
plot(rho, 1 - A_B, type = 'b', col = "black", ylim = c(0,1), ylab = "Estimated Lower Bounds")
lines(rho, (1 - A_1)*(1- A_2), type = 'b', col = "red")
legend("bottomright", bty = "n", legend = c("(1 - A_B)", "(1 - A_1)(1 - A_2)"), 
	col = c("black", "red"), lty = 1, cex=0.65)
dev.off()

