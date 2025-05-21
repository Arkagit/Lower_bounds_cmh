set.seed(1234)

q = function(y1, x1){
	out1 = exp(-log(2* pi) - log(det(prop_cov))/2 - t(y1 - x1)%*%solve(prop_cov)%*%(y1 - x1)/2)
	return(out1)
}

f = function(x2){
	out2 = exp(-log(2* pi) - log(det(tar_cov))/2 - t(x2)%*%solve(tar_cov)%*%(x2)/2)
	return(out2)
}

rho = c(-1.99, -1.75, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 1.75, 1,99)
n = 1e5
output = list()
init_point = c(-2, 2)
prop_cov = matrix(c(3, 0, 0, 3), nrow = 2, byrow = TRUE)

for(i in 1:length(rho)){
	cov_mat = matrix(c(1, rho[i], rho[i], 4), nrow = 3, byrow = TRUE)
