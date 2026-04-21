# log proposal pdf for block MH
log.q = function(eval, curr, prop_cov){
	out1 = - t(eval - curr)%*%solve(prop_cov)%*%(eval - curr)/2
	return(out1)
}

# target pdf for block MH
log.f = function(w_eval, tar_mean, tar_cov)
{
	out2 = - t(w_eval - tar_mean)%*%solve(tar_cov)%*%(w_eval - tar_mean)/2
	return(out2)
}

# Block \alpha_block
alpha_block = function(w_init, tar_mean, tar_cov, prop_cov)
{		
	prop = rnorm(2, w_init, sqrt(diag(prop_cov)))

	ratio = exp(log.f(prop, tar_mean, tar_cov) - log.f(w_init, tar_mean, tar_cov))
	alpha <- min(1, ratio)
	return(alpha)
}


# CMH1 \alpha_cmh
alpha_cmh = function(w_init, x1, tar_mean, tar_cov)
{	
	x_init <- w_init[1]
	y_init <- w_init[2]
	curr = c(x1, y_init)

	prop = c(x1, rnorm(1, y_init, sqrt(prop_cov[2,2])) )

	ratio = exp(log.f(prop, tar_mean, tar_cov) - log.f(curr, tar_mean, tar_cov))
	alpha <- min(1, ratio)
	return(alpha)
}


# Proximal Augmented \alpha_p
log.f2 = function(w_init, z, tar_mean, tar_cov, lambda)
{
	out2 = log.f(w_init, tar_mean, tar_cov) - (1/(2 * lambda))*(sum((z - w_init)^2))
	return(out2)
}

alpha_proximal = function(w_init, z, tar_mean, tar_cov, prop_cov, lambda)
{
	prop = rnorm(2, w_init, sqrt(diag(prop_cov)))

	ratio = exp(log.f2(prop, z, tar_mean, tar_cov, lambda) - log.f2(w_init, z, tar_mean, tar_cov, lambda))
	alpha <- min(1, ratio)
	return(alpha)
}






