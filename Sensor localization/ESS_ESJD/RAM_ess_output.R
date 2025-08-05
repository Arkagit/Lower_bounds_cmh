############### RAM transition used at each iteration
############### target function is defined as a function "target"
############### See line 156 for example.
set.seed(1234)  

library(foreach)
library(doParallel)
library(matrixcalc)


source("RAM_ess_func_rcpp.R")
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

eps = 10^(-100)

###################### Experimental Output
samp_size = c(1e4, 5e4)
repet = 1e2
nsim = 2e1

parallel::detectCores()
n.cores <- 8
doParallel::registerDoParallel(cores = n.cores)

initial.loc = c(0.5748, 0.0991, 0.2578, 0.8546, 
               0.9069, 0.3651, 0.1350, 0.0392)

j.scale = rep(1.08, 4)

  mat_RAM = list()
  mat_MH = list()
  mat_RAB_2coin = list()

  ess_RAM = list()
  ess_MH = list()
  ess_RAB_2coin = list()

  comp_ess_RAM = list()
  comp_ess_MH = list()
  comp_ess_RAB_2coin = list()

  esjd_RAM = list()
  esjd_MH = list()
  esjd_RAB_2coin = list()

alpha_list = list()

alpha_list = foreach(k = 1:repet)%dopar%{
  
    print("a")
    test_RAM = MHwG.RAM(initial.loc, runif(8), jump.scale = j.scale, 
                                Ob, Os, Xb, Xs, Yb, Ys, 
                                n.sample = max(samp_size), n.burn = 0)
    print("b")

    test_MH = MHwG.Metro(initial.loc = initial.loc, jump.scale = j.scale, 
                                 Ob, Os, Xb, Xs, Yb, Ys, n.sample = max(samp_size), n.burn = 0)

    print("c")

    test_RAB_2coin = MHwG.RAM.2coin.Barker(initial.loc, j.scale, 
             Ob, Os, Xb, Xs, Yb, Ys, n.sample = max(samp_size), n.burn = 0)

    for(k1 in 1:length(samp_size)){
      
      mat_RAM[[k1]] = mcse.multi(test_RAM$x[1:samp_size[k1], ], method = "bm", r = 1)$cov

      if (is.positive.definite(mat_RAM[[k1]]+ t(mat_RAM[[k1]]))) {

        ess_RAM[[k1]] = multiESS(test_RAM$x[1:samp_size[k1], ], mat_RAM[[k1]])

        comp_ess_RAM[[k1]] = component_ess(test_RAM$x[1:samp_size[k1], ], mat_RAM[[k1]])

      } else {
        ess_RAM[[k1]] = 0
      }

      esjd_RAM[[k1]] = esjd(test_RAM$x[1:samp_size[k1], ])

      print("d")

      mat_MH[[k1]] = mcse.multi(test_MH$x[1:samp_size[k1], ], method = "bm", r = 1)$cov

      if (is.positive.definite(mat_MH[[k1]]+t(mat_MH[[k1]]))) {

        ess_MH[[k1]] = multiESS(test_MH$x[1:samp_size[k1], ], mat_MH[[k1]])

        comp_ess_MH[[k1]] = component_ess(test_MH$x[1:samp_size[k1], ], mat_MH[[k1]])

      } else {

        ess_MH[[k1]] = 0

      }

      esjd_MH[[k1]] = esjd(test_MH$x[1:samp_size[k1], ])

      print("e")

      mat_RAB_2coin[[k1]] = mcse.multi(test_RAB_2coin$x[1:samp_size[k1], ], method = "bm", r = 1)$cov

      if (is.positive.definite(mat_RAB_2coin[[k1]]+ t(mat_RAB_2coin[[k1]]))) {

        ess_RAB_2coin[[k1]] = multiESS(test_RAB_2coin$x[1:samp_size[k1], ], mat_RAB_2coin[[k1]])

        comp_ess_RAB_2coin[[k1]] = component_ess(test_RAB_2coin$x[1:samp_size[k1], ], mat_RAB_2coin[[k1]])

      } else {
        ess_RAB_2coin[[k1]] = 0
      }

      esjd_RAB_2coin[[k1]] = esjd(test_RAB_2coin$x[1:samp_size[k1], ])

      print("f")
    }

    print(k)

    out = list(mat_RAM, ess_RAM, comp_ess_RAM, esjd_RAM, mat_MH, ess_MH, comp_ess_MH, esjd_MH, mat_RAB_2coin, ess_RAB_2coin, comp_ess_RAB_2coin, esjd_RAB_2coin) 
    out

  }






save(repet, j.scale, samp_size, nsim, simnum, jumping.scale, alpha_list, file = "RAB_ess_data.Rdata")


