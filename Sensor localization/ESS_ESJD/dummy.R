set.seed(11)  

library(foreach)
library(doParallel)
library(matrixcalc)


source("RAM_ess_func_rcpp.R")

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


samp_size = 1e6
initial.loc = c(0.5748, 0.0991, 0.2578, 0.8546, 
               0.9069, 0.3651, 0.1350, 0.0392)

j.scale = rep(1.08, 4)

parallel::detectCores()
n.cores <- 1
doParallel::registerDoParallel(cores = n.cores)

data_list = list()

data_list = foreach(k = 1:3)%dopar%{

	if(k == 1){

		mcdata = MHwG.RAM(initial.loc, runif(8), jump.scale = j.scale, 
                                Ob, Os, Xb, Xs, Yb, Ys, 
                                n.sample = samp_size, n.burn = 0)

	}else if(k == 2){

		mcdata = MHwG.Metro(initial.loc = initial.loc, jump.scale = j.scale, 
                                 Ob, Os, Xb, Xs, Yb, Ys, n.sample = samp_size, n.burn = 0)

	}else{

		mcdata = MHwG.RAM.2coin.Barker(initial.loc, j.scale, 
             Ob, Os, Xb, Xs, Yb, Ys, n.sample = samp_size, n.burn = 0)
	}
mcdata

}

save(data_list, samp_size, j.scale, file = "location_data.Rdata")



CRAM_data = data_list[[1]][[1]]

CMH_data = data_list[[2]][[1]]

CRAB_data = data_list[[3]][[1]]

pdf("Location_CRAM.pdf", width = 20, height = 20)
par(mfrow = c(2,2), mar = c(6, 8, 6, 4))  

plot(CRAM_data[,1], CRAM_data[,2], type = "p", xlab = "x11", ylab = "x12"
  , main = "Location of x1", xlim = range(CRAM_data[,1]), ylim = range(CRAM_data[,2]),
  cex.lab = 3,cex.main = 3,cex.axis = 2, cex = 2)

plot(CRAM_data[,3], CRAM_data[,4], type = "p", xlab = "x21", ylab = "x22"
  , main = "Location of x2", xlim = range(CRAM_data[,3]), ylim = range(CRAM_data[,4]),
  cex.lab = 3,cex.main = 3,cex.axis = 2, cex = 2)

plot(CRAM_data[,5], CRAM_data[,6], type = "p", xlab = "x31", ylab = "x32"
  , main = "Location of x3", xlim = range(CRAM_data[,5]), ylim = range(CRAM_data[,6]),
  cex.lab = 3,cex.main = 3,cex.axis = 2, cex = 2)

plot(CRAM_data[,7], CRAM_data[,8], type = "p", xlab = "x41", ylab = "x42"
  , main = "Location of x4", xlim = range(CRAM_data[,7]), ylim = range(CRAM_data[,8]),
  cex.lab = 3,cex.main = 3,cex.axis = 2, cex = 2)

dev.off()


pdf("Location_CMH.pdf", width = 20, height = 20)
par(mfrow = c(2,2), mar = c(6, 8, 6, 4))  

plot(CMH_data[,1], CMH_data[,2], type = "p", xlab = "x11", ylab = "x12"
  , main = "Location of x1", xlim = range(CMH_data[,1]), ylim = range(CMH_data[,2]),
  cex.lab = 3,cex.main = 3,cex.axis = 2, cex = 2)

plot(CMH_data[,3], CMH_data[,4], type = "p", xlab = "x21", ylab = "x22"
  , main = "Location of x2", xlim = range(CMH_data[,3]), ylim = range(CMH_data[,4]),
  cex.lab = 3,cex.main = 3,cex.axis = 2, cex = 2)

plot(CMH_data[,5], CMH_data[,6], type = "p", xlab = "x31", ylab = "x32"
  , main = "Location of x3", xlim = range(CMH_data[,5]), ylim = range(CMH_data[,6]),
  cex.lab = 3,cex.main = 3,cex.axis = 2, cex = 2)

plot(CMH_data[,7], CMH_data[,8], type = "p", xlab = "x41", ylab = "x42"
  , main = "Location of x4", xlim = range(CMH_data[,7]), ylim = range(CMH_data[,8]),
  cex.lab = 3,cex.main = 3,cex.axis = 2, cex = 2)

dev.off()


pdf("Location_CRAB.pdf", width = 20, height = 20)
par(mfrow = c(2,2), mar = c(6, 8, 6, 4))  

plot(CRAB_data[,1], CRAB_data[,2], type = "p", xlab = "x11", ylab = "x12"
  , main = "Location of x1", xlim = range(CRAB_data[,1]), ylim = range(CRAB_data[,2]),
  cex.lab = 3,cex.main = 3,cex.axis = 2, cex = 2)

plot(CRAB_data[,3], CRAB_data[,4], type = "p", xlab = "x21", ylab = "x22"
  , main = "Location of x2", xlim = range(CRAB_data[,3]), ylim = range(CRAB_data[,4]),
  cex.lab = 3,cex.main = 3,cex.axis = 2, cex = 2)

plot(CRAB_data[,5], CRAB_data[,6], type = "p", xlab = "x31", ylab = "x32"
  , main = "Location of x3", xlim = range(CRAB_data[,5]), ylim = range(CRAB_data[,6]),
  cex.lab = 3,cex.main = 3,cex.axis = 2, cex = 2)

plot(CRAB_data[,7], CRAB_data[,8], type = "p", xlab = "x41", ylab = "x42"
  , main = "Location of x4", xlim = range(CRAB_data[,7]), ylim = range(CRAB_data[,8]),
  cex.lab = 3,cex.main = 3,cex.axis = 2, cex = 2)

dev.off()



