set.seed(1234)

rho = c(-1.98, -1.75, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 1.75, 1.98)
n = 1e5
output = list()

for(i in 1:length(rho)){
tar_cov = matrix(c(1, rho[i], rho[i], 4), nrow = 2)
prop_cov = matrix(c(5, 0, 0, 5), nrow = 2)

q = function(y1, x1){
	out1 = exp(-log(2* pi) - log(det(prop_cov))/2 - t(y1 - x1)%*%solve(prop_cov)%*%(y1 - x1)/2)
	return(out1)
}

f = function(x2){
	out2 = exp(-log(2* pi) - log(det(tar_cov))/2 - t(x2)%*%solve(tar_cov)%*%(x2)/2)
	return(out2)
}


# A_c simulation
x = c(0, 0)
A_B = 0
A_2 = 0

for(j in 1:n){
y = c(rnorm(1, x[1], sqrt(5)), rnorm(1, x[2], sqrt(5)))

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
y = rnorm(1, x[2], sqrt(5))

r = exp(log(dnorm(y, x[2] - rho[i]*x[1]/tar_cov[1, 1], sqrt(4 - (rho[i]^(2))))) + log(dnorm(x[2], y, sqrt(5))) -
 log(dnorm(x[2], x[2] - rho[i]*x[1]/tar_cov[1,1], sqrt(4 - (rho[i]^(2))))) - log(dnorm(y, x[2], sqrt(5))))

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
#lines(rho, 1 - ub_B,type = 'b', col = "purple")
#lines(rho, 1 - ub_2,type = 'b', col = "skyblue")

legend("bottomright", bty = "n", legend = c("(1 - A_B)", "(1 - A_2)"), 
	col = c("black", "red"), lty = 1, cex=0.65)
dev.off()
#########################################

set.seed(1234)

rho = 0.5
sigma = c(0.1, 1, 5, 10, 50, 100)
n = 1e5
output = list()

for(i in 1:length(sigma)){
tar_cov = matrix(c(1, rho, rho, 4), nrow = 2)
prop_cov = matrix(c(sigma[i], 0, 0, sigma[i]), nrow = 2)

q = function(y1, x1){
    out1 = exp(-log(2* pi) - log(det(prop_cov))/2 - t(y1 - x1)%*%solve(prop_cov)%*%(y1 - x1)/2)
    return(out1)
}

f = function(x2){
    out2 = exp(-log(2* pi) - log(det(tar_cov))/2 - t(x2)%*%solve(tar_cov)%*%(x2)/2)
    return(out2)
}


# A_c simulation
x = c(0, 0)
A_B = 0
A_2 = 0

for(j in 1:n){
y = c(rnorm(1, x[1], sqrt(sigma[i])), rnorm(1, x[2], sqrt(sigma[i])))

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
y = rnorm(1, x[2], sqrt(5))

r1 = exp(log(dnorm(y, x[2] - rho*x[1]/tar_cov[1, 1], sqrt(4 - (rho^(2))))) + log(dnorm(x[2], y, sqrt(sigma[i]))) -
 log(dnorm(x[2], x[2] - rho*x[1]/tar_cov[1,1], sqrt(4 - (rho^(2))))) - log(dnorm(y, x[2], sqrt(sigma[i]))))

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

for(i in 1:length(sigma)){
    A_B = c(A_B, 1 - output[[i]][1])
    A_2 = c(A_2, 1 - output[[i]][2])
    ub_B = c(ub_B, 1 - output[[i]][3])
    ub_2 = c(ub_2, 1 - output[[i]][4])
    ub_B2 = c(ub_B2, 1 - output[[i]][5])
}

pdf("Acceptance_Probability_1.pdf", height = 6, width = 8)
par(mar = c(5.1, 4.8, 4.1, 2.1))
plot(sigma, A_2, type = 'b', col = "black", ylim = range(A_2), ylab = "Lower Bound")


#legend("topright", bty = "n", legend = c("A_B", "A_2", "Upper Bound (Block)", "Upper Bound A_2"), 
#    col = c("black", "red", "purple", "skyblue"), lty = 1, cex=0.65)
dev.off()




#lines(rho, A_B/A_2,type = 'b', col = "orange")
#lines(rho, ub_B2,type = 'b', col = "brown")
###########################


load("localanc_French_2way_1-3_18-22_148_60_0.99_100.Rdata")
str(localanc)
localanc[[5]][1,1:6,60]
Average = rep(0.2, 5)
str(c(localanc[[5]][,3,7]))


load("MOSAIC_RESULTS_TRAIN_noEM/localanc_UKB_500_Train__5way_1-500_11-11_5574_60_0.99_100.RData")
 
HAP_train=localanc
DIP_train=dip(localanc)

#DIP_train
str(DIP_train)
#lapply(DIP_train, mean,2)



Average=list()
for (i in c(1:42))
     {

        Average[[i]]=apply(HAP_train[[1]][,,i],1,mean)
     }
str(Average)
s

#Average_matrix=matrix(unlist(Average), ncol = length(Average), byrow = TRUE)
#Average_matrix

rm(localanc)
load("MOSAIC_RESULTS_TEST_noEM/localanc_UKB_500Test__5way_1-500_11-11_5574_60_0.99_100.RData")
HAP_test=localanc
DIP_test=dip(localanc)
#DIP_test[[1]][,1,1]
HAP_test[[1]][,,20]
str(hap_iid)

library(philentropy)
JSD=matrix(0,nrow=42, ncol = 6)
EU_dist=matrix(NA,nrow=42,ncol=6)
for (j in c(1:6))
{
        for (i in c(1:42))

        {
                avg=Average[[i]]
                m = matrix(NA, nrow = 2, ncol = length(avg))
                hap_iid=c(HAP_test[[1]][,j,i])
                m[1,]=avg
                m[2,]=hap_iid
#               print(m)
                m = as.matrix(m)
        #       m1=apply(m,2,c)
        #       print(str(m1))
                print(1)
                JSD[i,j]=JSD(m, unit = "log2")
                print(2)
                EU_dist[i,j]=sqrt(sum((m[1,]-m[2,])**2))
        }

}
"get_lc.R" 55L, 1078C 


JSD


#################################################
set.seed(1234)

AR = function(nsamp, rho){
	x = numeric(nsamp)
	x[1] = 0
	eps = rnorm(nsamp)
	for (i in 2:nsamp) {
		x[i] = rho*x[i-1] + eps[i-1]
	}
	return(x)
}

n = c(1e2, 1e3, 1e4, 1e5, 1e6)

rho = 0.5

data = AR(max(n), rho)

for(j in 1:length(n)){
	non_fft_acf[j] = acf(data, lag.max = n[j]-1)
	fft_acf[j] = 
}

non_fft = acf(data, lag.max = )



