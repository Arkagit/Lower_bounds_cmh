########### All the Bivariate Normal related CMH functions
######### Bivariate Gibbs
gibbs_output = function(x1, x2){
	new_x1 = rnorm(1, 0, sqrt(1/(2*(1 + x2^(2)))))
	new_x2 = rnorm(1, 0, sqrt(1/(2*(1 + new_x1^(2)))))
	return(c(new_x1, new_x2))
}



###################### 2-component CMH
######## (X1|X2) ~ Normal(0, 1/(2(1+x2^2))) Gibbs step
######## (X2|X1) ~ Normal(0, 1/(2(1+x1^2))) MH step by Normal(x1, h^2)
CMH_output = function(x1, x2, sd){
  new_x1 = rnorm(1, 0, sqrt(1/(2*(1 + x2^(2)))))
  z = rnorm(1, x2, sd)
  log_alpha = dnorm(z, 0, sqrt(1/(2*(1 + new_x1^(2)))), log = TRUE) + dnorm(x2, x2, sd, log = TRUE) -
              dnorm(x2, 0, sqrt(1/(2*(1 + new_x1^(2)))), log = TRUE) - dnorm(z, x2, sd, log = TRUE)
  u = runif(1)
  new_x2 = x2
  if(u <= exp(log_alpha)){
    new_x2 = z
  }
  alpha = min(1, exp(log_alpha))
  return(c(new_x1, new_x2, alpha))
}

###################### 2-component CMH
######## (X1|X2) ~ Normal(0, 1/(2(1+x2^2))) MH step by Normal(x2, h^2)
######## (X2|X1) ~ Normal(0, 1/(2(1+x1^2))) MH step by Normal(x1, h^2)
CMH2_output = function(x1, x2, sd){
	#
	  z1 = rnorm(1, x1, sd)
	  log_alpha1 = dnorm(z1, 0, sqrt(1/(2*(1 + x2^(2)))), log = TRUE) + dnorm(x1, x1, sd, log = TRUE) -
	              dnorm(x1, 0, sqrt(1/(2*(1 + x2^(2)))), log = TRUE) - dnorm(z1, x1, sd, log = TRUE)
	  u1 = runif(1)
	  new_x1 = x1
	  if(u1 <= exp(log_alpha1)){
	    new_x1 = z1
	  }
	  alpha1 = min(1, exp(log_alpha1))
	#
	  z2 = rnorm(1, x2, sd)
	  log_alpha2 = dnorm(z2, 0, sqrt(1/(2*(1 + new_x1^(2)))), log = TRUE) + dnorm(x2, x2, sd, log = TRUE) -
	              dnorm(x2, 0, sqrt(1/(2*(1 + new_x1^(2)))), log = TRUE) - dnorm(z2, x2, sd, log = TRUE)
	  u2 = runif(1)
	  new_x2 = x2
	  if(u2 <= exp(log_alpha2)){
	    new_x2 = z2
	  }
	  alpha2 = min(1, exp(log_alpha2))
	  return(c(new_x1, new_x2, alpha1, alpha2))
}

##########################################
################################# lower bound of lower bound function





