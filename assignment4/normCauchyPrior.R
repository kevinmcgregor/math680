#A function to draw n independent samples from the posterior dist coming from a normal dist whose
#mean has a cauchy prior.

# n is the number of samples to generate
# x is the observed datum
normCauchyPrior = function(n=1, x) {
  if (n<1 | n!=floor(n)) stop("n must be a positive integer")
  if (!is.numeric(x)) stop("x must be numeric")
  
  samp = rep(0,n)
  n.trial = 0
  n.accept = 0
  u = 0
  while(n.accept < n) {
    n.trial = n.trial + 1
    #Get an observation from the trial distribution (Cauchy(0,1))
    trialobs = rcauchy(1,0,1)
    u = runif(1)
    #If condition met for trial... accept observation
    if (u < exp(-1/2*(trialobs-x)^2) ) {
      n.accept = n.accept + 1
      samp[n.accept] = trialobs
    }
  }
  
  return(list(samp=samp, n.trial=n.trial))
}

#Another function for importance sampling using the same trial density

importNormCauchy = function(n, x) {
  if (n<1 | n!=floor(n)) stop("n must be a positive integer")
  if (!is.numeric(x)) stop("x must be numeric")
  
  Z = rcauchy(n,0,1)
  
  return(sum(Z*exp(-(Z-x)^2))/sum(exp(-1/2*(Z-x)^2)))
}






