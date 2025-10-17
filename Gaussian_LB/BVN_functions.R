# proposal pdf for block MH
q = function(y1, x1, prop_cov){
	out1 = exp(-log(2* pi) - log(det(prop_cov))/2 - t(y1 - x1)%*%solve(prop_cov)%*%(y1 - x1)/2)
	return(out1)
}

# target pdf for block MH
f = function(x2, tar_mean, tar_cov){
	out2 = exp(-log(2* pi) - log(det(tar_cov))/2 - t(x2 - tar_mean)%*%solve(tar_cov)%*%(x2 - tar_mean)/2)
	return(out2)
}

# Block \alpha_B
alpha_B = function(x_init, tar_mean, tar_cov, prop_cov){
		
	y = c(rnorm(1, x_init[1], sqrt(prop_cov[1,1])), rnorm(1, x_init[2], sqrt(prop_cov[2,2])))

	r = exp(log(f(y, tar_mean, tar_cov)) + log(q(x_init, x_init, prop_cov)) - log(f(x_init, tar_mean, tar_cov)) - log(q(y, x_init, prop_cov)))

	if(r < 1){
		alpha = r
	}else{  
		alpha = 1
	}

	return(alpha)
}


# CMH1 \alpha_mh
alpha_mh = function(x_init, tar_mean, tar_cov, prop_cov){
		
	y = rnorm(1, x_init[2], sqrt(prop_cov[2,2]))

	r = exp(log(dnorm(y, x_init[2] + tar_cov[1,2]*x_init[1]/tar_cov[1, 1], sqrt(tar_cov[2,2] - (tar_cov[1, 2]^(2))/tar_cov[1,1]))) +
	    log(dnorm(x_init[2], y, sqrt(prop_cov[2,2]))) -
		 log(dnorm(x_init[2], x_init[2] + tar_cov[1, 2]*x_init[1]/tar_cov[1,1], sqrt(tar_cov[2,2] - (tar_cov[1, 2]^(2))/tar_cov[1,1]))) - 
		 log(dnorm(y, x_init[2], sqrt(prop_cov[2,2]))))

	if(r < 1){
		alpha = r
	}else{  
		alpha = 1
	}

	return(alpha)
}


# Proximal Augmented \alpha_p

alpha_p = function(x_init, tar_mean, tar_cov, prop_cov, lambda){

	Omega = solve(solve(tar_cov) + diag(c(1, 1))/lambda)

	z = c(rnorm(1, x_init[1], sqrt(lambda)), rnorm(1, x_init[2], sqrt(lambda)))

	pop_mean = Omega%*%z/lambda

	dummy = c(rnorm(1, x_init[1], sqrt(prop_cov[1,1])), rnorm(1, x_init[2], sqrt(prop_cov[2,2])))

	r = exp(- 0.5*t(dummy - pop_mean)%*%solve(Omega)%*%(dummy - pop_mean) +
				0.5*t(x_init - x_init)%*%solve(Omega)%*%(x_init - x_init) - 
				0.5*t(x_init - pop_mean)%*%solve(prop_cov)%*%(x_init - pop_mean) +
				0.5*t(dummy - x_init)%*%solve(prop_cov)%*%(dummy - x_init))

	if(r < 1){
		alpha = r
	}else{  
		alpha = 1
	}

	return(alpha)
}

