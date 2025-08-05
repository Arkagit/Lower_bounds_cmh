######## Target joint posterior density

norm2 <- function(loca, locb) {
  sqrt(sum((loca - locb)^2))
}

l.target <- function(loc, R = 0.3, sigma = 0.02, Ob, Os, Xb, Xs, Yb, Ys) {

  First.term <- NULL
  for (i in 1 : 2) {
    TEMP <- sapply(1 : 4, function(j) {
      exp(-norm2(Xb[i, ], loc[(2 * j -1) : (2 * j)])^2 / 2 / R^2 * Ob[j, i]) *
      (1 - exp(-norm2(Xb[i, ], loc[(2 * j -1) : (2 * j)])^2 / 2 / R^2))^(1 - Ob[j, i]) 
    })
    First.term <- c(First.term, TEMP)
  }

  Second.term <- NULL
  for (i in 1 : 3) {
    TEMP <- sapply((i + 1) : 4, function(j) {
      exp(-norm2(loc[(2 * i -1) : (2 * i)], 
                 loc[(2 * j -1) : (2 * j)])^2 / 2 / R^2 * Os[i, j]) *
      (1 - exp(-norm2(loc[(2 * i -1) : (2 * i)], 
                      loc[(2 * j -1) : (2 * j)])^2 / 2 / R^2))^(1 - Os[i, j]) 
    })
    Second.term <- c(Second.term, TEMP)
  }

  First.obs.term <- NULL
  for (i in 1 : 2) {
    TEMP <- sapply(1 : 4, function(j) {
      dnorm(Yb[j, i], mean = norm2(Xb[i, ], loc[(2 * j -1) : (2 * j)]), 
                      sd = sigma)^Ob[j, i]
    })
    First.obs.term <- c(First.obs.term, TEMP)
  }

  Second.obs.term <- NULL
  for (i in 1 : 3) {
    TEMP <- sapply((i + 1) : 4, function(j) {
      dnorm(Ys[i, j], mean = norm2(loc[(2 * i -1) : (2 * i)], 
                                   loc[(2 * j -1) : (2 * j)]), 
                      sd = sigma)^Os[i, j]
    })
    Second.obs.term <- c(Second.obs.term, TEMP)
  }

  log.lik <- sum(log(c(First.term, Second.term, First.obs.term, Second.obs.term)))
  post <- log.lik + sum(dnorm(loc, mean = rep(0, 8), sd = rep(10, 8), log = TRUE))
  post

}


######## RAM

ram.kernel <- function(current.location, current.aux, loc.number, scale) {

  eps <- 10^(-100)
  accept <- 0
  x.c <- current.location 
  log.x.c.den <- l.target(x.c, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  x.c.den <- exp(log.x.c.den)
  z.c <- current.aux
  log.z.c.den <- l.target(z.c, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  z.c.den <- exp(log.z.c.den)


  # downhill
  x.p1 <- x.c
  x.p1[(2 * loc.number - 1) : (2 * loc.number)] <- x.p1[(2 * loc.number - 1) : (2 * loc.number)] + 
                                                   rnorm(2, 0, scale)
  log.x.p1.den <- l.target(x.p1, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  x.p1.den <- exp(log.x.p1.den)
  N.d <- 1
  while (-rexp(1) > log(x.c.den + eps) - log(x.p1.den + eps)) {
    x.p1 <- x.c
    x.p1[(2 * loc.number - 1) : (2 * loc.number)] <- x.p1[(2 * loc.number - 1) : (2 * loc.number)] + 
                                                     rnorm(2, 0, scale)
    log.x.p1.den <- l.target(x.p1, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
    x.p1.den <- exp(log.x.p1.den)
    N.d <- N.d + 1
  }

  # uphill
  x.p2 <- x.p1
  x.p2[(2 * loc.number - 1) : (2 * loc.number)] <- x.p2[(2 * loc.number - 1) : (2 * loc.number)] + 
                                                   rnorm(2, 0, scale)
  log.x.p2.den <- l.target(x.p2, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  x.p2.den <- exp(log.x.p2.den)
  
  N.u <- 1
  while (-rexp(1) > log(x.p2.den + eps) - log(x.p1.den + eps)) {
    x.p2 <- x.p1
    x.p2[(2 * loc.number - 1) : (2 * loc.number)] <- x.p2[(2 * loc.number - 1) : (2 * loc.number)] + 
                                                     rnorm(2, 0, scale)
    log.x.p2.den <- l.target(x.p2, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
    x.p2.den <- exp(log.x.p2.den)
    N.u <- N.u + 1
  }
  
  # downhill for N.d
  N.dz <- 1     # number of total downhill trials for estimate
  z <- x.p2
  z[(2 * loc.number - 1) : (2 * loc.number)] <- z[(2 * loc.number - 1) : (2 * loc.number)] + 
                                                rnorm(2, 0, scale)
  log.z.den <- l.target(z, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  z.den <- exp(log.z.den)
  while (-rexp(1) > log(x.p2.den + eps) - log(z.den + eps)) {
    z <- x.p2
    z[(2 * loc.number - 1) : (2 * loc.number)] <- z[(2 * loc.number - 1) : (2 * loc.number)] + 
                                                  rnorm(2, 0, scale)
    log.z.den <- l.target(z, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
    z.den <- exp(log.z.den)
    N.dz <- N.dz + 1
  }

  # accept or reject the proposal
  min.nu <- min(1, (x.c.den + eps) / (z.c.den + eps))
  min.de <- min(1, (x.p2.den + eps) / (z.den + eps))

  l.mh <- log.x.p2.den + log(min.nu) - log.x.c.den - log(min.de)

  #if (is.infinite(log.x.p2.den) && is.infinite(log.x.c.den) && log.x.p2.den == log.x.c.den) {
   # l.mh <- log(min.nu) - log(min.de)  # or 0, depending on your context
   # print(1)
   # } else {
   # l.mh <- log.x.p2.den + log(min.nu) - log.x.c.den - log(min.de)
  #}

  if (l.mh > - rexp(1)) {
    x.c <- x.p2
    z.c <- z
    accept <- 1
  }

  acceptance_rate = min(1, exp(l.mh))
  accp_barker = exp(l.mh)/(1 + exp(l.mh))

  c(x.c, z.c, N.d, N.u, N.dz, accept, acceptance_rate, accp_barker)
}

MHwG.RAM <- function(initial.loc, initial.aux, jump.scale, 
             Ob, Os, Xb, Xs, Yb, Ys, n.sample = 10, n.burn = 10) {

  
  n.total <- n.sample + n.burn
  accept <- matrix(0, nrow = n.total, ncol = 4)
  alpha = matrix(0, nrow = n.total, ncol = 4)
  alpha_bark = matrix(0, nrow = n.total, ncol = 4)
  out <- matrix(NA, nrow = n.total, ncol = 8)
  loc.t <- initial.loc
  aux.t <- initial.aux
  Nd <- matrix(NA, nrow = n.total, ncol = 4)
  Nu <- matrix(NA, nrow = n.total, ncol = 4)
  Nz <- matrix(NA, nrow = n.total, ncol = 4)
  
  for (i in 1 : n.total) {
    for (j in 1 : 4) {
      TEMP <- ram.kernel(loc.t, aux.t, j, jump.scale[j])
      loc.t <- TEMP[1 : 8]
      aux.t <- TEMP[9 : 16]
      Nd[i, j] <- TEMP[17]
      Nu[i, j] <- TEMP[18]
      Nz[i, j] <- TEMP[19]
      accept[i, j] <- TEMP[20]
      alpha[i, j] = TEMP[21]
      alpha_bark[i, j] = TEMP[22]
    }
    out[i, ] <- loc.t
  }
  
  list(x = out[-c(1 : n.burn), ], 
       accept = accept[-c(1 : n.burn), ],
       N.d = Nd[-c(1 : n.burn), ],
       N.u = Nu[-c(1 : n.burn), ],
       N.z = Nz[-c(1 : n.burn), ],
       alpha_RAM = alpha,
       alpha_bark = alpha_bark
       )

}

######## Metropolis
Metro.kernel <- function(current.location, loc.number, jump.scale) {

  accept <- 0
  x.p <- current.location
  x.p[(2 * loc.number - 1) : (2 * loc.number)] <- x.p[(2 * loc.number - 1) : (2 * loc.number)] + 
                                                  rnorm(2, 0, jump.scale)
  

  l.new = l.target(x.p, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  l.current = l.target(current.location, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)

  #if (is.infinite(l.new) && is.infinite(l.current) && l.new == l.current) {
  #  l.metro <- 0  # or 0, depending on your context
  #  } else {
  #  l.metro <- l.new - l.current
  #}

  l.metro <- l.new - l.current

  if (l.metro >  -rexp(1)) {
    current.location <- x.p
    accept <- 1
  }
  acceptance_rate = min(1, exp(l.metro))
  c(current.location, accept, acceptance_rate)

}


MHwG.Metro <- function(initial.loc, jump.scale, Ob, Os, Xb, Xs, Yb, Ys, n.sample = 10, n.burn = 10) {

  n.total <- n.sample + n.burn
  accept <- matrix(0, nrow = n.total, ncol = 4)
  out <- matrix(NA, nrow = n.total, ncol = 8)
  loc.t <- initial.loc
  
  for (i in 1 : n.total) {

    for (j in 1 : 4) {
      TEMP <- Metro.kernel(loc.t, j, jump.scale[j])
      loc.t <- TEMP[1 : 8]
      accept[i, j] <- TEMP[9]
    }

    out[i, ] <- loc.t

  }

  list(x = out[-c(1 : n.burn), ], accept = accept[-c(1 : n.burn), ], alpha_Metro = accept)

}

############################# RAM BARKER
two_coin_algorithm <- function(c.xy, mu, mu.p, loc.number, scale) {
  repeat {
    # Step 1: Draw C1 ~ Bern(cy / (cx + cy))
    C1 <- rbinom(1, 1, 1/(1 + c.xy))
    
    if (C1 == 1) {
      # Step 3: Draw C2 ~ Bern(py)
      M <- mu
      M[(2 * loc.number - 1) : (2 * loc.number)] <- M[(2 * loc.number - 1) : (2 * loc.number)] + 
                                                   rnorm(2, 0, scale)
      py = min(1, (exp(l.target(mu, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)) + eps) /
             (exp(l.target(M, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)) + eps))

      C2 <- rbinom(1, 1, py)
      if (C2 == 1) {
        return(1)  # Output 1
      } 
      # else repeat from beginning
    } else {
      M <- mu.p
      M[(2 * loc.number - 1) : (2 * loc.number)] <- M[(2 * loc.number - 1) : (2 * loc.number)] + 
                                                   rnorm(2, 0, scale)
      px = min(1, (exp(l.target(mu.p, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)) + eps) /
             (exp(l.target(M, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)) + eps))
      # Step 9: Draw C2 ~ Bern(px)
      C2 <- rbinom(1, 1, px)
      if (C2 == 1) {
        return(0)  # Output 0
      }
      # else repeat from beginning
    }
  }
}


ram.2coin.barker.kernel = function(current.location, loc.number, scale, simnum){
  eps <- 10^(-100)
  check <- rep(0, simnum)
  x.c <- current.location 
  log.x.c.den <- l.target(x.c, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  x.c.den <- exp(log.x.c.den)
  


  # downhill
  x.p <- x.c
  x.p[(2 * loc.number - 1) : (2 * loc.number)] <- x.p[(2 * loc.number - 1) : (2 * loc.number)] + 
                                                   rnorm(2, 0, scale)
  log.x.p.den <- l.target(x.p, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  x.p.den <- exp(log.x.p.den)

  c.xy <- exp(log.x.c.den - log.x.p.den)
  
  for(s in 1:simnum){
    check[s] = two_coin_algorithm(c.xy, x.c, x.p, loc.number, scale)
  }
  
  check_prob = mean(check); check_prob

  if(rbinom(1, 1, check_prob) == 1){
    x.c = x.p
  }

  return(c(x.c, check_prob))
}



MHwG.RAM.2coin.Barker <- function(initial.loc, jump.scale, 
             Ob, Os, Xb, Xs, Yb, Ys, n.sample = 10, n.burn = 10){

  n.total <- n.sample + n.burn
  alpha = matrix(0, nrow = n.total, ncol = 4)
  out <- matrix(NA, nrow = n.total, ncol = 8)
  loc.t <- initial.loc
  
  for (i in 1 : n.total) {
    for (j in 1 : 4) {
      TEMP <- ram.2coin.barker.kernel(loc.t, j, jump.scale[j], simnum)
      loc.t <- TEMP[1 : 8]
      alpha[i, j] = TEMP[9]
    }
    out[i, ] <- loc.t
  }
  
  list(x = out[-c(1 : n.burn), ], 
       alpha_RAM_2coin = alpha
       )
}