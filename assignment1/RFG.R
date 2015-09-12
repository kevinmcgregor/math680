#Random function generator from Friedman 2001

#Main random function generator
RFG = function(N=100,p=10) {
  #Input data
  x = generate_x(N, rep(0,p), diag(p))
  
  
}

#Generates observed variables x from multivariate normal distribution.
# N specifies the sample size
# mu is a vector of means for the multivar normal distribution
# sigma is a square matrix of the same dimension as mu
generate_x = function(N, mu, sigma) {
  if (N <= 0) stop("N must be positive")
  if (length(mu)!=NROW(sigma) | length(mu)!=NCOL(sigma)) stop("mu and sigma are not of the same dimension")
  require(MASS)
  return(mvrnorm(N, mu, sigma))
}

#Generates random error term
# N is the sample size
generate_re = function(N) {
  if (N <= 0) stop("N must be positive")
  return(rnorm(N))
}



