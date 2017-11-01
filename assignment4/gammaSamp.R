# A function to generate observations from Gamma(alpha, beta) using rejection sampling,
# assuming alpha>1
#returns a list with "samp" being the sample from a gamma distribution, and "n.trial" being the number
#of trial samples generated before the sample size was reached.
gammaSamp = function(n=1, alpha=2, beta=1) {
  if (n<1 | n!=floor(n)) stop("n must be a positive integer")
  if (alpha<=1) stop("alpha must be greater than 1")
  if (beta<=0) stop("beta must be greater than 0")
  
  #Find the floor of alpha
  fla = floor(alpha)
  
  
  #Defining c for rejection sampling
  c = beta^(alpha-fla)*gamma(fla)/gamma(alpha)*((alpha-fla)/(beta-1))^(alpha-fla)*exp(fla-alpha)
  
  #Vector to hold final sample
  samp = rep(0,n)
  n.trial = 0
  n.accept = 0
  u = 0
  while(n.accept < n) {
    n.trial = n.trial + 1
    #Get an observation from the trial distribution
    trialobs = getTrialObs(alpha, fla)
    u = runif(1)
    #If condition met for trial... accept observation
    if (u < dgamma(trialobs,alpha,beta)/(c*dgamma(trialobs,fla,1))) {
      n.accept = n.accept + 1
      samp[n.accept] = trialobs
    }
  }
  
  #Multiply all observations by 1/beta to get gamma(alpha,beta)
  samp = samp/beta
  
  return(list(samp=samp, n.trial=n.trial))
}

#Function to get a trial observation
getTrialObs = function(alpha, fla) {
  #Generate an observation from an exponential dist using inverse cdf
  r = runif(fla, 0, 1)
  #Sum to get gamma(floor(alpha),1) RV
  gamobs = sum(-log(1-r))
  
  return(gamobs)
}



