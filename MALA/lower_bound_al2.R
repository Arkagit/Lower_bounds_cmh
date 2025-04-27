set.seed(1234)

library(latex2exp)
library(doParallel)
library(latex2exp)

############## Function
f = function(x_2, sd, mc.size, r1){
  dum = 0
  foo = rep(0, mc.size)
  for(t1 in 1:mc.size)
  {
    x_1 = r1[t1] 
    z <- rnorm(mc.size, x_2 - x_2*sd^2 - x_2*x_1^2*sd^2, sd)
    foo = exp(dnorm(z, 0, sqrt(1/(2*(1 + x_1^2))), log = TRUE) - dnorm(x_2, 0, sqrt(1/(2*(1 + x_1^2))), log = TRUE) +
      dnorm(z, x_2 - x_2*sd^2 - x_2*x_1^2*sd^2, sd, log = TRUE) - dnorm(x_2, z - z*sd^2 - z*x_1^2*sd^2, sd, log = TRUE))
    dum <- dum + mean(pmin(foo, 1))
  }
  return(dum/(mc.size))
}


############ Output and Plot
set.seed(1234)
sequence = 1e2
h = c(0.001, 1, 100) # proposal standard deviation

l2 = seq(-5, 5, length.out = sequence)# Grid

lb_list = list()

for(i in 1:length(h)){
  lb_list[[i]] = rep(0, length = sequence)
}

mc.size <- 1e4

parallel::detectCores()
n.cores <- 4
doParallel::registerDoParallel(cores = n.cores)

lb_list = foreach(k = 1:length(h))%dopar%{
  sd = h[k]
  dumat = rep(0, length(l2))
  for(i2 in 1:length(l2)){
      print(c(k, i2))
      x_2 = l2[i2]
      r1 = rnorm(mc.size, 0, sqrt(1/(2*(1+x_2^2))))
      dumat[i2] = f(x_2, sd, mc.size, r1)
  }
  dumat
}


pdf(paste("As_MALA.pdf"), height = 6, width = 18)
par(mfrow = c(1,3))
for(i in 1:length(h)){
  plot(l2, lb_list[[i]], type = "l", ylim = c(0, 1),
       ylab = "A*", xlab = TeX(r'($X_2$)'), main = paste("h = ",h[i]))
}
dev.off()

pdf(paste("Hist_al2.pdf"), height = 6, width = 8)
par(mfrow = c(1,3))
for(i in 1:length(h)){
  hist(lb_list[[i]], xlim = c(0, 1), main = paste("h = ",h[i]), xlab = "A*")
}
dev.off()