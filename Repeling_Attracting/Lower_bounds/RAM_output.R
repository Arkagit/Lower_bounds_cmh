############### RAM transition used at each iteration
############### target function is defined as a function "target"
############### See line 156 for example.
set.seed(11)

library(foreach)
library(doParallel)

source("RAM_functions.R")
############################ Example 3: Sensor location problem

######## Data

# Observation indicators from the fifth sensor (1st column) to the first four sensors
# and those from the sixth sensor (2nd column) to the first four sensors.
Ob <- matrix(c(1, 0, 1, 0, 1, 0, 1, 0), ncol = 2)

# Observation indicators among the first four sensors. 
Os <- matrix(c(0, 0, 0, 1,
               0, 0, 1, 1,
               0, 1, 0, 0,
               1, 1, 0, 0), ncol = 4)

# Each row indicates the location of the known sensors (5th and 6th).
Xb <- matrix(c(0.5, 0.3, 0.3, 0.7), ncol = 2)

# Each row indicates the location of the unknown sensors (1st, 2nd, 3rd, and 4th).
Xs <- matrix(c(0.5748, 0.0991, 0.2578, 0.8546, 
               0.9069, 0.3651, 0.1350, 0.0392), ncol = 2)

# The observed distances from the fifth sensor (1st column) to the first four sensors
# and those from the sixth sensor (2nd column) to the first four sensors.
Yb <- matrix(c(0.6103, 0, 0.2995, 0, 
               0.3631, 0, 0.5656, 0), ncol = 2)

# Observed distances among the first four sensors.
Ys <- matrix(c(0, 0, 0, 0.9266,
               0, 0, 0.2970, 0.8524,
               0, 0.2970, 0, 0,
               0.9266, 0.8524, 0, 0), ncol = 4)

eps = 10^(-308)

###################### Experimental Output

repet = 1e4
nsim = 2e1
simnum = 1e2

low_lim = 0.0001
up_lim = 2
jumping.scale = matrix(c(seq(low_lim, up_lim, length.out = nsim), seq(low_lim, up_lim, length.out = nsim),
                  seq(low_lim, up_lim, length.out = nsim), seq(low_lim, up_lim, length.out = nsim)), 
                    nrow = nsim, byrow = FALSE)

parallel::detectCores()
n.cores <- 4
doParallel::registerDoParallel(cores = n.cores)

initial.loc = c(0.5748, 0.0991, 0.2578, 0.8546, 
               0.9069, 0.3651, 0.1350, 0.0392)

alpha_list = list(0)

alpha_list = foreach(k = 1:nsim)%dopar%{
  print(k)
  j.scale = jumping.scale[k,]

  mat_RAM = list()
  mat_RAM_bark = list()
  mat_MH = list()
  mat_RAM_2coin = list()

  for(i in 1:repet){

    test = MHwG.RAM(initial.loc, runif(8), jump.scale = j.scale, 
                                  Ob, Os, Xb, Xs, Yb, Ys, 
                                  n.sample = 1, n.burn = 0)
    mat_RAM[[i]] = test$alpha_RAM

    mat_RAM_bark[[i]] = test$alpha_bark

    mat_MH[[i]] = MHwG.Metro(initial.loc, jump.scale = j.scale, 
                                   Ob, Os, Xb, Xs, Yb, Ys, n.sample = 1, n.burn = 0)$alpha_Metro
    
    mat_RAM_2coin[[i]] = MHwG.RAM.2coin.Barker(initial.loc, jump.scale = j.scale, 
                                  Ob, Os, Xb, Xs, Yb, Ys, 
                                  n.sample = 1, n.burn = 0)$alpha_RAM_2coin

    print(i)
  }

  out = list(mat_RAM, mat_RAM_bark, mat_MH, mat_RAM_2coin) 
  out
}





save(repet, nsim, simnum, jumping.scale, alpha_list, file = "RAM_data.Rdata")


