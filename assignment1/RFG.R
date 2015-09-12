#Random function generator from Friedman 2001

#Main random function generator
RFG = function(N=100,p=10,l=20) {
  #Input data
  x = generate_x(N, rep(0,p), diag(p))
  error = generate_re(N)
  
  #Creating the random function
  a = generate_a(l)
  pl = generate_pl(l)
  
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

#Generate the coefficients a1...al in the sum term
# l is the number of terms in the sum... defaults to 20
generate_a = function(l=20) {
  return(runif(l,-1,1))
}

#Generate r
generate_r = function() {
  return(rexp(1,0.5))
}

#Generate the p_l
# p is the number of predictors
# l is the number of terms in the sum... defaults to 20
generate_pl = function(p,l=20) {
  rvec = replicate(l, generate_r())
  pl = pmin(floor(1.5+rvec),p)
  return(pl)
}

