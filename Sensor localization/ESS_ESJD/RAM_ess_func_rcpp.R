############################ Example 3: Sensor location problem
library(mcmcse)
library(mcmc)
library(Rcpp)
library(phonTools)
library(RcppArmadillo)
library(RcppEigen)

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



######## Target joint posterior density





cppFunction('
  double norm2(NumericVector loca, NumericVector locb) {
    int n = loca.size();
    double sum_sq = 0.0;

    for (int i = 0; i < n; ++i) {
      double diff = loca[i] - locb[i];
      sum_sq += diff * diff;
    }

    return std::sqrt(sum_sq);
  }
')


cppFunction('
#include <Rcpp.h>
using namespace Rcpp;

// Helper: Euclidean norm
double norm2(NumericVector a, NumericVector b) {
  double sum = 0.0;
  for (int i = 0; i < a.size(); ++i) {
    double d = a[i] - b[i];
    sum += d * d;
  }
  return std::sqrt(sum);
}

// [[Rcpp::export]]
double l_target(NumericVector loc, double R, double sigma, 
                NumericMatrix Ob, NumericMatrix Os, 
                NumericMatrix Xb, NumericMatrix Xs, 
                NumericMatrix Yb, NumericMatrix Ys) {

  std::vector<double> First_term;

  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 4; ++j) {
      NumericVector subloc = loc[Range(2 * j, 2 * j + 1)];
      double d = norm2(Xb(i, _), subloc);
      double term = std::exp(-d * d / (2 * R * R)) * Ob(j, i) +
                    std::pow(1 - std::exp(-d * d / (2 * R * R)), 1 - Ob(j, i));
      First_term.push_back(term);
    }
  }

  std::vector<double> Second_term;
  for (int i = 0; i < 3; ++i) {
    for (int j = i + 1; j < 4; ++j) {
      NumericVector loc_i = loc[Range(2 * i, 2 * i + 1)];
      NumericVector loc_j = loc[Range(2 * j, 2 * j + 1)];
      double d = norm2(loc_i, loc_j);
      double term = std::exp(-d * d / (2 * R * R)) * Os(i, j) +
                    std::pow(1 - std::exp(-d * d / (2 * R * R)), 1 - Os(i, j));
      Second_term.push_back(term);
    }
  }

  std::vector<double> First_obs_term;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 4; ++j) {
      NumericVector subloc = loc[Range(2 * j, 2 * j + 1)];
      double d = norm2(Xb(i, _), subloc);
      double val = Yb(j, i);
      double pdf = (1.0 / (sigma * std::sqrt(2 * M_PI))) * 
                   std::exp(-0.5 * std::pow((val - d) / sigma, 2));
      First_obs_term.push_back(std::pow(pdf, Ob(j, i)));
    }
  }

  std::vector<double> Second_obs_term;
  for (int i = 0; i < 3; ++i) {
    for (int j = i + 1; j < 4; ++j) {
      NumericVector loc_i = loc[Range(2 * i, 2 * i + 1)];
      NumericVector loc_j = loc[Range(2 * j, 2 * j + 1)];
      double d = norm2(loc_i, loc_j);
      double val = Ys(i, j);
      double pdf = (1.0 / (sigma * std::sqrt(2 * M_PI))) * 
                   std::exp(-0.5 * std::pow((val - d) / sigma, 2));
      Second_obs_term.push_back(std::pow(pdf, Os(i, j)));
    }
  }

  double log_lik = 0.0;
  for (double v : First_term) log_lik += std::log(v);
  for (double v : Second_term) log_lik += std::log(v);
  for (double v : First_obs_term) log_lik += std::log(v);
  for (double v : Second_obs_term) log_lik += std::log(v);

  // Prior: standard normal N(0, 10^2) for each loc element
  double prior = 0.0;
  for (int i = 0; i < loc.size(); ++i) {
    prior += R::dnorm(loc[i], 0.0, 10.0, true);  // log = true
  }

  return log_lik + prior;
}
')

########## RAM

ram_kernel <- function(current.location, current.aux, loc.number, scale) {

  eps <- 10^(-308)
  accept <- 0
  x.c <- current.location 
  log.x.c.den <- l_target(x.c, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  x.c.den <- exp(log.x.c.den)
  z.c <- current.aux
  log.z.c.den <- l_target(z.c, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  z.c.den <- exp(log.z.c.den)

  # downhill
  x.p1 <- x.c
  x.p1[(2 * loc.number - 1) : (2 * loc.number)] <- x.p1[(2 * loc.number - 1) : (2 * loc.number)] + 
                                                   rnorm(2, 0, scale)
  log.x.p1.den <- l_target(x.p1, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  x.p1.den <- exp(log.x.p1.den)
  N.d <- 1
  while (-rexp(1) > log(x.c.den + eps) - log(x.p1.den + eps)) {
    x.p1 <- x.c
    x.p1[(2 * loc.number - 1) : (2 * loc.number)] <- x.p1[(2 * loc.number - 1) : (2 * loc.number)] + 
                                                     rnorm(2, 0, scale)
    log.x.p1.den <- l_target(x.p1, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
    x.p1.den <- exp(log.x.p1.den)
    N.d <- N.d + 1
  }

  # uphill
  x.p2 <- x.p1
  x.p2[(2 * loc.number - 1) : (2 * loc.number)] <- x.p2[(2 * loc.number - 1) : (2 * loc.number)] + 
                                                   rnorm(2, 0, scale)
  log.x.p2.den <- l_target(x.p2, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  x.p2.den <- exp(log.x.p2.den)
  N.u <- 1
  while (-rexp(1) > log(x.p2.den + eps) - log(x.p1.den + eps)) {
    x.p2 <- x.p1
    x.p2[(2 * loc.number - 1) : (2 * loc.number)] <- x.p2[(2 * loc.number - 1) : (2 * loc.number)] + 
                                                     rnorm(2, 0, scale)
    log.x.p2.den <- l_target(x.p2, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
    x.p2.den <- exp(log.x.p2.den)
    N.u <- N.u + 1
  }

  # downhill for N.d
  N.dz <- 1     # number of total downhill trials for estimate
  z <- x.p2
  z[(2 * loc.number - 1) : (2 * loc.number)] <- z[(2 * loc.number - 1) : (2 * loc.number)] + 
                                                rnorm(2, 0, scale)
  log.z.den <- l_target(z, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  z.den <- exp(log.z.den)
  while (-rexp(1) > log(x.p2.den + eps) - log(z.den + eps)) {
    z <- x.p2
    z[(2 * loc.number - 1) : (2 * loc.number)] <- z[(2 * loc.number - 1) : (2 * loc.number)] + 
                                                  rnorm(2, 0, scale)
    log.z.den <- l_target(z, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
    z.den <- exp(log.z.den)
    N.dz <- N.dz + 1
  }

  # accept or reject the proposal
  min.nu <- min(1, (x.c.den + eps) / (z.c.den + eps))
  min.de <- min(1, (x.p2.den + eps) / (z.den + eps))
  l.mh <- log.x.p2.den - log.x.c.den + log(min.nu) - log(min.de)

  if (l.mh > -rexp(1)) {
    x.c <- x.p2
    z.c <- z
    accept <- 1
  }

  c(x.c, z.c, N.d, N.u, N.dz, accept)
}




MHwG.RAM <- function(initial.loc, initial.aux, jump.scale, 
             Ob, Os, Xb, Xs, Yb, Ys, n.sample = 10, n.burn = 10) {

  n.total <- n.sample + n.burn
  accept <- matrix(0, nrow = n.total, ncol = 4)
  out <- matrix(NA, nrow = n.total, ncol = 8)
  loc.t <- initial.loc
  aux.t <- initial.aux
  Nd <- matrix(NA, nrow = n.total, ncol = 4)
  Nu <- matrix(NA, nrow = n.total, ncol = 4)
  Nz <- matrix(NA, nrow = n.total, ncol = 4)
  
  for (i in 1 : n.total) {
    for (j in 1 : 4) {
      TEMP <- ram_kernel(loc.t, aux.t, j, jump.scale[j])
      loc.t <- TEMP[1 : 8]
      aux.t <- TEMP[9 : 16]
      Nd[i, j] <- TEMP[17]
      Nu[i, j] <- TEMP[18]
      Nz[i, j] <- TEMP[19]
      accept[i, j] <- TEMP[20]
    }
    out[i, ] <- loc.t
  }
  
  list(x = out, 
       accept = accept[-c(1 : n.burn), ])

}

################## MH

Metro_kernel <- function(current.location, loc.number, jump.scale) {

  accept <- 0
  x.p <- current.location
  x.p[(2 * loc.number - 1) : (2 * loc.number)] <- x.p[(2 * loc.number - 1) : (2 * loc.number)] + 
                                                  rnorm(2, 0, jump.scale)
  l.metro <- l_target(x.p, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys) - 
             l_target(current.location, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)

  if (l.metro >  -rexp(1)) {
    current.location <- x.p
    accept <- 1
  }
  c(current.location, accept)

}


MHwG.Metro <- function(initial.loc, jump.scale, Ob, Os, Xb, Xs, Yb, Ys, n.sample = 10, n.burn = 0) {

  #print(Sys.time())
  n.total <- n.sample + n.burn
  accept <- matrix(0, nrow = n.total, ncol = 4)
  out <- matrix(NA, nrow = n.total, ncol = 8)
  loc.t <- initial.loc
  
  for (i in 1 : n.total) {

    for (j in 1 : 4) {
      TEMP <- Metro_kernel(loc.t, j, jump.scale[j])
      loc.t <- TEMP[1 : 8]
      accept[i, j] <- TEMP[9]
    }

    out[i, ] <- loc.t

  }

  #print(Sys.time())
  list(x = out, accept = accept[-c(1 : n.burn), ])

}


############### 2 coin Barker

library(Rcpp)

cppFunction('
#include <Rcpp.h>
using namespace Rcpp;

// Compute Euclidean norm squared
double norm2(NumericVector a, NumericVector b) {
  double sum = 0.0;
  for (int i = 0; i < a.size(); ++i) {
    double diff = a[i] - b[i];
    sum += diff * diff;
  }
  return std::sqrt(sum);
}

// Target log-posterior
double l_target(NumericVector loc, double R, double sigma,
                NumericMatrix Ob, NumericMatrix Os,
                NumericMatrix Xb, NumericMatrix Xs,
                NumericMatrix Yb, NumericMatrix Ys) {

  std::vector<double> First_term;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 4; ++j) {
      NumericVector xbj = Xb(i, _);
      NumericVector locj = loc[Range(2 * j, 2 * j + 1)];
      double d = norm2(xbj, locj);
      double exp_term = std::exp(-d * d / (2 * R * R));
      double val = std::pow(exp_term, Ob(j, i)) *
                   std::pow(1 - exp_term, 1 - Ob(j, i));
      First_term.push_back(val);
    }
  }

  std::vector<double> Second_term;
  for (int i = 0; i < 3; ++i) {
    for (int j = i + 1; j < 4; ++j) {
      NumericVector loci = loc[Range(2 * i, 2 * i + 1)];
      NumericVector locj = loc[Range(2 * j, 2 * j + 1)];
      double d = norm2(loci, locj);
      double exp_term = std::exp(-d * d / (2 * R * R));
      double val = std::pow(exp_term, Os(i, j)) *
                   std::pow(1 - exp_term, 1 - Os(i, j));
      Second_term.push_back(val);
    }
  }

  std::vector<double> First_obs_term;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 4; ++j) {
      NumericVector xbj = Xb(i, _);
      NumericVector locj = loc[Range(2 * j, 2 * j + 1)];
      double d = norm2(xbj, locj);
      double obs = Yb(j, i);
      double dens = R::dnorm(obs, d, sigma, false);
      First_obs_term.push_back(std::pow(dens, Ob(j, i)));
    }
  }

  std::vector<double> Second_obs_term;
  for (int i = 0; i < 3; ++i) {
    for (int j = i + 1; j < 4; ++j) {
      NumericVector loci = loc[Range(2 * i, 2 * i + 1)];
      NumericVector locj = loc[Range(2 * j, 2 * j + 1)];
      double d = norm2(loci, locj);
      double obs = Ys(i, j);
      double dens = R::dnorm(obs, d, sigma, false);
      Second_obs_term.push_back(std::pow(dens, Os(i, j)));
    }
  }

  double log_lik = 0.0;
  for (double x : First_term) log_lik += std::log(x + 1e-100);
  for (double x : Second_term) log_lik += std::log(x + 1e-100);
  for (double x : First_obs_term) log_lik += std::log(x + 1e-100);
  for (double x : Second_obs_term) log_lik += std::log(x + 1e-100);

  double prior = 0.0;
  for (int i = 0; i < loc.size(); ++i) {
    prior += R::dnorm(loc[i], 0.0, 10.0, true);
  }

  return log_lik + prior;
}

// Two-coin Barker acceptance
// [[Rcpp::export]]
int two_coin_algorithm(double c_xy, NumericVector mu, NumericVector mu_p,
                       int loc_number, double scale,
                       NumericMatrix Ob, NumericMatrix Os,
                       NumericMatrix Xb, NumericMatrix Xs,
                       NumericMatrix Yb, NumericMatrix Ys) {
  double eps = 1e-100;
  NumericVector M = clone(mu);

  while (true) {
    int C1 = R::rbinom(1, 1.0 / (1.0 + c_xy));
    if (C1 == 1) {
      M = clone(mu);
      for (int i = 0; i < 2; ++i) {
        M[2 * (loc_number - 1) + i] += R::rnorm(0.0, scale);
      }
      double log_mu = l_target(mu, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys);
      double log_M = l_target(M, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys);
      double py = std::min(1.0, (std::exp(log_mu) + eps) / (std::exp(log_M) + eps));
      int C2 = R::rbinom(1, py);
      if (C2 == 1) return 1;
    } else {
      M = clone(mu_p);
      for (int i = 0; i < 2; ++i) {
        M[2 * (loc_number - 1) + i] += R::rnorm(0.0, scale);
      }
      double log_mup = l_target(mu_p, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys);
      double log_M = l_target(M, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys);
      double px = std::min(1.0, (std::exp(log_mup) + eps) / (std::exp(log_M) + eps));
      int C2 = R::rbinom(1, px);
      if (C2 == 1) return 0;
    }
  }
}
')


ram_2coin_barker_kernel = function(current.location, loc.number, scale){
  eps <- 10^(-100)
  x.c <- current.location 
  log.x.c.den <- l_target(x.c, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  x.c.den <- exp(log.x.c.den)
  


  # downhill
  x.p <- x.c
  x.p[(2 * loc.number - 1) : (2 * loc.number)] <- x.p[(2 * loc.number - 1) : (2 * loc.number)] + 
                                                   rnorm(2, 0, scale)
  log.x.p.den <- l_target(x.p, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  x.p.den <- exp(log.x.p.den)

  c.xy <- exp(log.x.c.den - log.x.p.den)
  
  
  decision = two_coin_algorithm(c.xy, x.c, x.p, loc.number, scale, Ob, Os, Xb, Xs, Yb, Ys)

  if(decision == 1){
    x.c = x.p
  }

  return(x.c)
}

MHwG.RAM.2coin.Barker <- function(initial.loc, jump.scale, 
             Ob, Os, Xb, Xs, Yb, Ys, n.sample = 10, n.burn = 0){

  n.total <- n.sample + n.burn
  out <- matrix(NA, nrow = n.total, ncol = 8)
  loc.t <- initial.loc
  
  for (i in 1 : n.total) {
    for (j in 1 : 4) {
      TEMP <- ram_2coin_barker_kernel(loc.t, j, jump.scale[j])
      loc.t <- TEMP
    }
    out[i, ] <- loc.t
  }
  
  list(x = out)
}

### Calculating ESJD 
Rcpp::cppFunction('
Rcpp::NumericVector esjd(Rcpp::NumericMatrix chain) {
  int T = chain.nrow();
  int d = chain.ncol();

  Rcpp::NumericVector indiv_sq_dists(d, 0.0);
  Rcpp::NumericVector row_sq_dists(T - 1, 0.0);

  for (int t = 0; t < T - 1; ++t) {
    for (int j = 0; j < d; ++j) {
      double diff = chain(t + 1, j) - chain(t, j);
      double sq = diff * diff;
      indiv_sq_dists[j] += sq;
      row_sq_dists[t] += sq;
    }
  }

  for (int j = 0; j < d; ++j) {
    indiv_sq_dists[j] /= (T - 1);
  }

  double total_esjd = std::accumulate(row_sq_dists.begin(), row_sq_dists.end(), 0.0) / (T - 1);

  Rcpp::NumericVector output(d + 1);
  for (int j = 0; j < d; ++j) {
    output[j] = indiv_sq_dists[j];
  }
  output[d] = total_esjd;

  return output;
}
')



component_ess = function(chain){
  comp_ess = rep(0, 4)
  for(i in 1:4){
    minichain = chain[, c(2*i -1, 2*i)]
    cov_est = mcse.multi(minichain, method = "bm", r =1)$cov
    comp_ess[i] = multiESS(minichain, cov_est)
  }
  return(comp_ess)
}

# Comp_ess = function(chain, mat){
#   out = rep(0, dim(chain)[2])
#   for(i in 1:dim(chain)[2]){
#     out[i] = dim(chain)[1]*var(chain[, i])/mat[i, i]
#   }
#   return(out)
# }



 # chain <- matrix(rnorm(1000 * 8), ncol = 8)
 # component_ess(chain)
# mat <- diag(rep(1, 8))

# Comp_ess_cpp(chain, mat)

# Comp_ess(chain, mat)


# ?var

# chain <- matrix(rnorm(1000 * 5), nrow = 1000, ncol = 5)

# # Call the C++ function
# esjd(chain)
