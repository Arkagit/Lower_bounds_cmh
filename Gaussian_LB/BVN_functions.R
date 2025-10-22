# log proposal pdf for block MH
log.q = function(y1, x1, prop_cov){
	out1 = - t(y1 - x1)%*%solve(prop_cov)%*%(y1 - x1)/2
	return(out1)
}

# target pdf for block MH
log.f = function(x2, tar_mean, tar_cov){
	out2 = - t(x2 - tar_mean)%*%solve(tar_cov)%*%(x2 - tar_mean)/2
	return(out2)
}

# Block \alpha_block
alpha_block = function(x_init, tar_mean, tar_cov, prop_cov){
		
	y = rnorm(2, x_init, sqrt(diag(prop_cov)))

	ratio = exp(log.f(y, tar_mean, tar_cov) - log.f(x_init, tar_mean, tar_cov))

	if(ratio < 1){
		alpha = ratio
	}else{  
		alpha = 1
	}

	return(alpha)
}


# CMH1 \alpha_cmh
alpha_cmh = function(x_init, x1, tar_mean, tar_cov){
		
	curr = c(x1, x_init[2])

	prop = c(x1, rnorm(1, x_init[2], sqrt(prop_cov[2,2])))

	ratio = exp(log.f(prop, tar_mean, tar_cov) - log.f(curr, tar_mean, tar_cov))

	if(ratio < 1){
		alpha = ratio
	}else{  
		alpha = 1
	}

	return(alpha)
}


# Proximal Augmented \alpha_p

log.f2 = function(x_init, z, tar_mean, tar_cov, lambda){

	out2 = log.f(x_init, tar_mean, tar_cov) - (1/(2 * lambda))*(sum((z - x_init)^2))
	
	return(out2)
}

alpha_proximal = function(x_init, z, tar_mean, tar_cov, prop_cov, lambda){

	y = rnorm(2, x_init, sqrt(diag(prop_cov)))

	ratio = exp(log.f2(y, z, tar_mean, tar_cov, lambda) - log.f2(x_init, z, tar_mean, tar_cov, lambda))

	if(ratio < 1){
		alpha = ratio
	}else{  
		alpha = 1
	}

	return(alpha)
}






