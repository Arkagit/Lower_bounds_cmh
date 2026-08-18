library(pracma)

# \phi(n) = \log (n)

design_mat_poor_cond <- function(n, d) {
  if (n <= 1) stop("n must be > 1 so that log(n) > 0")
  if (d < 2) stop("d must be at least 2")
  
  Sigma <- diag(c(rep(1, d - 1), 1 / log(n)))
  
  Z <- matrix(rnorm(n * d), nrow = n, ncol = d)
  
  Sigma_half <- diag(c(rep(1, d - 1), 1 / sqrt(log(n))))
  
  W <- Z %*% Sigma_half
  
  return(list(W = W, Sigma = Sigma))
}

design_mat_good_cond <- function(n, d) {
  if (n <= 1) stop("n must be > 1 so that log(n) > 0")
  if (d < 2) stop("d must be at least 2")
  
  Sigma <- diag(c(rep(1, d - 1), 1 / n))
  
  Z <- matrix(rnorm(n * d), nrow = n, ncol = d)
  
  Sigma_half <- diag(c(rep(1, d - 1), 1 / sqrt(n)))
  
  W <- Z %*% Sigma_half
  
  return(list(W = W, Sigma = Sigma))
}

VM_lower_bound = function(n, d, h, b, zeta_star, xi_star, W){
  WtW  = t(W)%*%W
  ev <- eigen(WtW, symmetric = TRUE, only.values = TRUE)$values

  lambda_max <- max(ev)
  lambda_min <- min(ev)
  kappa_WtW <- lambda_max / lambda_min

  c_1 = lambda_max/n
  c_2 = kappa_WtW/log(n)
  tilde_c = (c_1 / c_2)^d

  tilde_Y = W%*%xi_star + rnorm(n, 0, zeta_star)
  m_xi = solve(WtW + diag(rep(1, d)))%*%t(W)%*%tilde_Y
  brack = ((Norm(tilde_Y - W%*%m_xi, 2) + Norm(m_xi, 2))/2 + b)/n

  vm_lb = 1 - kappa_WtW*(log(n)/(h*n))^(d/2)*(brack)^(d/2)

  return(vm_lb)
}
