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


# A_c simulationa
x = c(0, 0)
A_B = 0
A_2 = 0

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

#A_B simulation
for(j in 1:n){
y = rnorm(1, x[2], sqrt(prop_cov[2,2]))

r = exp(log(dnorm(y, x[2] + rho[i]*x[1]/tar_cov[1, 1], sqrt(tar_cov[2,2] - (rho[i]^(2))/tar_cov[1,1]))) + 
    log(dnorm(x[2], y, sqrt(prop_cov[2,2]))) -
 log(dnorm(x[2], x[2] + rho[i]*x[1]/tar_cov[1,1], sqrt(tar_cov[2,2] - (rho[i]^(2))/tar_cov[1,1]))) - 
 log(dnorm(y, x[2], sqrt(prop_cov[2,2]))))

if(r > 1){
	alpha = 1
}else{
	alpha = r
}

A_2 = A_2 + alpha/n
}
print(2)

#Upper bound of A_B
ub_B = (det(tar_cov)/det(prop_cov))^(1/2) *
		exp(t(x)%*%solve(tar_cov)%*%x / 2)


#Upper bound of A_2
ub_2 = sqrt((tar_cov[2,2] - tar_cov[1,2]*tar_cov[2,1]/tar_cov[1,1])/(prop_cov[2,2])) *
		exp(((x[2] - x[1]*tar_cov[2,1]/tar_cov[1,1])^(2))/(2*prop_cov[2,2]))

#Upper bound of A_c/A_B
ub_B2 = (max(eigen(tar_cov[1,1])$values)/min(eigen(prop_cov[1,1])$values))^(1/2) *
		exp(t(x)%*%solve(tar_cov)%*%x / 2)

output[[i]] = c(A_B, A_2, ub_B, ub_2, ub_B2)

}

A_B = c();A_2 = c(); ub_B = c(); ub_2 = c(); ub_B2 = c()

for(i in 1:length(rho)){
	A_B = c(A_B, output[[i]][1])
	A_2 = c(A_2, output[[i]][2])
	ub_B = c(ub_B, output[[i]][3])
	ub_2 = c(ub_2, output[[i]][4])
	ub_B2 = c(ub_B2, output[[i]][5])
}

pdf("Gaussian_Lower_bound.pdf", height = 6, width = 8)
par(mar = c(5.1, 4.8, 4.1, 2.1))
plot(rho, 1 - A_B, type = 'b', col = "black", ylim = c(0,1), ylab = "Estimated Lower Bounds")
lines(rho, 1 - A_2, type = 'b', col = "red")
legend("bottomright", bty = "n", legend = c("Complete Block", "CMH1"), 
	col = c("black", "red"), lty = 1, cex=0.65)
dev.off()
