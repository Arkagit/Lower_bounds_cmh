set.seed(1234)

library(latex2exp)
library(doParallel)
library(latex2exp)

############## Function
f = function(x_2, sd, sequence, r1, gauss){
  dum = 0
  mc.size <- length(gauss)
  for(t1 in 1:mc.size)
  {
    x_1 = r1[t1] 
    z <- gauss
    foo = rbind(rep(1, mc.size), exp(-(1 + x_1^2)*(z^2 - x_2^2) -
            (z - x_2 + x_2*sd^2 + x_1*x_2^2*sd^2)^2/(2*sd^2) +
            (x_2 - z + z*sd^2 + z*x_1^2*sd^2)^2/(2*sd^2)) )
    dum <- dum + mean(apply(foo, 2, min))
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

mc.size <- 1e3

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
      gauss <- rnorm(mc.size, x_2, sd)
      dumat[i2] = f(x_2, sd, sequence, r1, gauss)
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