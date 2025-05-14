set.seed(1234)

library(ggplot2)
library(dplyr)
library(cowplot)
library(foreach)
library(doParallel)
library(latex2exp)
library(patchwork)
################################################## Function
f = function(x_2, sd, mc.size, r1){
  dum = 0
  foo = rep(0, mc.size)
  for(t1 in 1:mc.size)
  {
    x_1 = r1[t1] 
    z <- rnorm(mc.size, x_2 - x_2*sd^2 - x_2*x_1^2*sd^2, sd)
    foo = exp(dnorm(z, 0, sqrt(1/(2*(1 + x_1^2))), log = TRUE) - dnorm(x_2, 0, sqrt(1/(2*(1 + x_1^2))), log = TRUE) -
      dnorm(z, x_2 - x_2*sd^2 - x_2*x_1^2*sd^2, sd, log = TRUE) + dnorm(x_2, z - z*sd^2 - z*x_1^2*sd^2, sd, log = TRUE))
    dum <- dum + mean(pmin(foo, 1))
  }
  return(dum/(mc.size))
}

f1 = function(x_1, x_2, sd, r1){
  #foo1 = rep(0, length(r1))
  foo1 = pmin(exp(dnorm(r1, 0, sqrt(1/(2*(1 + x_2^2))), log = TRUE) - dnorm(x_1, 0, sqrt(1/(2*(1 + x_2^2))), log = TRUE) +
       dnorm(x_1, r1 - r1*sd^2 - r1*x_2^2*sd^2, sd, log = TRUE) - dnorm(r1, x_1 - x_1*sd^2 - x_1*x_2^2*sd^2, sd, log = TRUE)), 1)
  dum1 = mean(foo1)
  return(dum1)
}

f2 = function(x_1, x_2, sd, r2){
  foo2 = pmin(exp(dnorm(r2, 0, sqrt(1/(2*(1 + x_1^2))), log = TRUE) - dnorm(x_2, 0, sqrt(1/(2*(1 + x_1^2))), log = TRUE) +
       dnorm(x_2, r2 - r2*sd^2 - r2*x_1^2*sd^2, sd, log = TRUE) - dnorm(r2, x_2 - x_2*sd^2 - x_2*x_1^2*sd^2, sd, log = TRUE)), 1)
  dum2 = mean(foo2)
  return(dum2)
}



#################################################### CMH1 
sequence = 1e2
h = c(0.001, 1, 5) # proposal standard deviation

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


################################################## CMH2
sequence = 1e2
mc.size = 1e4
h = c(0.001, 1, 5) # proposal standard deviation

l1 = seq(-5, 5, length.out = sequence)
l2 = seq(-5, 5, length.out = sequence)

lb_list1 = list()
lb_list2 = list()

for(i in 1:length(h)){
  lb_list1[[i]] = matrix(0, nrow = length(l1), ncol = length(l2))
  lb_list2[[i]] = matrix(0, nrow = length(l1), ncol = length(l2))
}

parallel::detectCores()
n.cores <- 4
doParallel::registerDoParallel(cores = n.cores)


lb_list1 = foreach(k = 1:length(h))%dopar%{
  print(k)
  sd = h[k]
  dumat = matrix(0, nrow = length(l1), ncol = length(l2))
  for(i2 in 1:length(l2)){
    for(i1 in 1:length(l1)){
      x_1 = l1[i1]
      x_2 = l2[i2]
      r1 = rnorm(mc.size, x_1 - x_1*sd^2 - x_1*x_2^2*sd^2, sd)
      dumat[i1, i2] = f1(x_1, x_2, sd, r1)
    }
  }
  dumat
}

lb_list2 = foreach(k = 1:length(h))%dopar%{
  print(k)
  sd = h[k]
  dumat = matrix(0, nrow = length(l1), ncol = length(l2))
  for(i2 in 1:length(l2)){
    for(i1 in 1:length(l1)){
      x_1 = l1[i1]
      x_2 = l2[i2]
      r2 = rnorm(mc.size, x_2 - x_2*sd^2 - x_2*x_1^2*sd^2, sd)
      dumat[i1, i2] = f2(x_1, x_2, sd, r2)
    }
  }
  dumat
}



save(lb_list, lb_list1, lb_list2, h, l1, l2, sequence, mc.size, file = "lower_bound.Rdata")

#################################################################


