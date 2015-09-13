#Random function generator from Friedman 2001
#Kevin McGregor - MATH680 HW1 - Due Sept 16th, 2015

#Main function for random function generator
#Returns data frame with first column y and remaining columns x_1 through x_p
RFG = function(N=100,p=10,l=20) {
  if (!is.numeric(N) | !is.numeric(p) | !is.numeric(l)) stop("N, p, and l must be numeric")
  if (N <= 0 | N != round(N)) stop("N must be a positive integer")
  if (p <= 0 | p != round(p)) stop("p must be a positive integer")
  if (l <= 0 | l != round(l)) stop("l must be a positive integer")
  #Input data
  x = generate_x(N, rep(0,p), diag(p))
  error = generate_re(N)
  
  #Creating the random function
  a = generate_a(l)
  pl = generate_pl(p,l)
  gl = generate_gl(x,pl)
  fx = gl %*% a
  
  y = fx + error
  return(data.frame(y,x))
}

#Generates observed variables x from multivariate normal distribution.
# N specifies the sample size
# mu is a vector of means for the multivar normal distribution
# sigma is a square matrix of the same dimension as mu
generate_x = function(N, mu, sigma) {
  require(MASS)
  return(mvrnorm(N, mu, sigma))
}

#Generates random error term
# N is the sample size
generate_re = function(N) {
  return(rnorm(N))
}

#Generate the coefficients a1...al in the sum term
# l is the number of terms in the overall sum
generate_a = function(l) {
  return(runif(l,-1,1))
}

#Generate r
generate_r = function() {
  return(rexp(1,0.5))
}

#Generate the p_l
# p is the number of predictors
# l is the number of terms in the overall sum
generate_pl = function(p,l) {
  rvec = replicate(l, generate_r())
  pl = pmin(floor(1.5+rvec),p)
  return(pl)
}

#Generate the z_l
# x is the input data
# pl is the number of (permuted) columns of x to include in each zl
generate_zl = function(x, pl) {
  l = length(pl)
  p = NCOL(x)
  
  zl = list()
  perm=NULL
  for (i in 1:l) {
    perm = sample(p)
    zl[i] = list(x[,perm[1:pl[i]]])
  }
  return(zl)
}

#Generate mu_l
# N is the number of samples
# pl is the number of (permuted) columns of x included in each zl
generate_mul = function(N,pl) {
  l = length(pl)
  mul = list()
  for (i in 1:l) {
    #mu_l is generated using the same distribution as each row of x
    mul[i] = list(generate_x(N,rep(0,pl[i]),diag(pl[i])))
  }
  return(mul)
}

#Generates D_l
# pl is the number of (permuted) columns of x included in each zl
# returns a list of length l
generate_Dl = function(pl) {
  l = length(pl)
  
  Dl = list()
  for (i in 1:l) {
    if (pl[i]>1){
      Dl[i] = list(diag(runif(pl[i],0.1,2)))
    } else {
      Dl[i] = list(runif(1,0.1,2))
    }
  }
  return(Dl)
}

#Generate Vl
# pl is the number of (permuted) columns of x included in each zl
# returns a list of length l
generate_Vl = function(pl) {
  l = length(pl)
  Dl = generate_Dl(pl)
  
  Vl = list()
  curQ = NULL
  curDl = NULL
  for (i in 1:l) {
    #Random orthonormal matrix
    curQ = genQ(pl[i])
    #Also need to get Dl to proper form for matrix mult
    curDl = matrix(unlist(Dl[i]),pl[i],pl[i])
    Vl[i] = list(curQ %*% curDl %*% t(curQ))
  }
  
  return(Vl)
}

#Generate gl
# x is the input data
# pl is the number of (permuted) columns of x included in each zl
# returns an Nxl matrix
generate_gl = function(x,pl) {
  N = NROW(x)
  l = length(pl)
  zl = generate_zl(x,pl)
  mul = generate_mul(N,pl)
  Vl = generate_Vl(pl)
  
  curZl = NULL
  curmul = NULL
  curVl = NULL
  gl = matrix(0, N, l)
  for (i in 1:l) {
    curZl = matrix(unlist(zl[i]),ncol=pl[i],byrow=TRUE)
    curmul = matrix(unlist(mul[i]),ncol=pl[i],byrow=TRUE)
    curVl = matrix(unlist(Vl[i]),ncol=pl[i],byrow=TRUE)
    for (j in 1:N) {
      gl[j,i] = exp(-1/2 * t(curZl[j,]-curmul[j,]) %*% curVl %*% (curZl[j,]-curmul[j,]) )
    }
  }
  return(gl)
}



