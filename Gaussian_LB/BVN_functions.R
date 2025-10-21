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

# Block \alpha_B
alpha_B = function(x_init, tar_mean, tar_cov, prop_cov){
		
	#y = rnorm(1, x_init[1], sqrt(prop_cov[1,1])), rnorm(1, x_init[2], sqrt(prop_cov[2,2])))
	y = rnorm(2, 0, sqrt(2))

	r = exp(log.f(y, tar_mean, tar_cov) + log.q(x_init, x_init, prop_cov) - log.f(x_init, tar_mean, tar_cov) - log.q(y, x_init, prop_cov))

	if(r < 1){
		alpha = r
	}else{  
		alpha = 1
	}

	return(r)
}


# CMH1 \alpha_mh
alpha_mh = function(x_init, x1, tar_mean, tar_cov, prop_cov){
		
	y = rnorm(1, x_init[2], sqrt(prop_cov[2,2]))

	r = exp(dnorm(y, x_init[2] + tar_cov[1,2]*x_init[1]/tar_cov[1, 1], sqrt(tar_cov[2,2] - (tar_cov[1, 2]^(2))/tar_cov[1,1]), log = TRUE) +
		     dnorm(x_init[2], y, sqrt(prop_cov[2,2]), log = TRUE) -
			 dnorm(x_init[2], x_init[2] + tar_cov[1, 2]*x_init[1]/tar_cov[1,1], sqrt(tar_cov[2,2] - (tar_cov[1, 2]^(2))/tar_cov[1,1]), log = TRUE) - 
			 dnorm(y, x_init[2], sqrt(prop_cov[2,2]), log = TRUE))

	if(r < 1){
		alpha = r
	}else{  
		alpha = 1
	}

	return(r)
}


# Proximal Augmented \alpha_p

log.f2 = function(x_init, z, inv_tar_cov, prop_cov, lambda){
	out2 = - t(x_init)%*%(inv_tar_cov + diag(c(1, 1))/ lambda)%*%x_init/2 + t(z)%*%x_init/lambda
	return(out2)
}

alpha_p = function(x_init, z, tar_cov, prop_cov, lambda){

	inv_tar_cov = solve(tar_cov)

	#dummy = c(rnorm(1, x_init[1], sqrt(prop_cov[1,1])), rnorm(1, x_init[2], sqrt(prop_cov[2,2])))
	y = rnorm(2, 0, sqrt(2))

	r = exp(log.f2(y, z, inv_tar_cov, prop_cov, lambda) + log.q(x_init, x_init, prop_cov) - 
		log.f2(x_init, z, inv_tar_cov, prop_cov, lambda) - log.q(y, x_init, prop_cov))

	# r = exp(- 0.5*t(y - pop_mean)%*%Omega_inv%*%(y - pop_mean) +
	# 			0.5*t(x_init - x_init)%*%solve(prop_cov)%*%(x_init - x_init) - 
	# 			0.5*t(x_init - pop_mean)%*%solve(prop_cov)%*%(x_init - pop_mean) +
	# 			0.5*t(y - x_init)%*%Omega_inv%*%(y - x_init))

	if(r < 1){
		alpha = r
	}else{  
		alpha = 1
	}

	return(r)
}






