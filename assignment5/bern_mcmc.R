#Function to do MCMC sampling of posterior distribution (A5 Q1)

bern_mcmc = function(X,y,nsamp,Ti,a0,b0,vb,quiet=FALSE) {

  n = length(X)
  p = NCOL(X[[1]])
  
  #Data structures to hold parameter samples
  beta.post = matrix(0,nsamp,p)
  u.post = matrix(0,nsamp,n)
  vu.post = rep(0,nsamp)
  
  #Initialize parameters
  beta = rep(0,p)
  u = rep(1,n)
  vu = 1
  
  #Variance in normal distribution (which is the trial distribution)
  trial_var = 0.1
  
  beta.new = rep(0,p)
  u.new = rep(0,n)
  vu.new = 0
  diff_loglik = 0
  log_acc_prob = 0
  unif_val = 0
  #Get samples
  for (i in 1:nsamp) {
    if (!quiet & i%%1000==0) cat(i, "\n")
    #Get new values of params from trial distribution
    beta.new = rnorm(p, beta, sqrt(trial_var))
    u.new = rnorm(n, u, sqrt(trial_var))
    vu.new = rnorm(1, vu, sqrt(trial_var))
    
    #Difference in log-likelihood
    diff_loglik = getLogLike(X,y,Ti,beta.new,u.new,vu.new,vb,a0,b0) - getLogLike(X,y,Ti,beta,u,vu,vb,a0,b0)
    #Get log of acceptance prob
    log_acc_prob = min(0, diff_loglik)
    
    #Decide to accept/reject
    unif_val = runif(1)
    if (log(unif_val)<log_acc_prob) {
      #Accept... use these values in next iteration
      beta = beta.new
      u = u.new
      vu = vu.new 
    } 
    #Saving current sample
    beta.post[i,] = beta
    u.post[i,] = u
    vu.post[i] = vu
    
  }
  
  return(list(beta.post=beta.post,u.post=u.post,vu.post=vu.post))
}

#Inverse logit function
i.logit = function(x){exp(x)/(1+exp(x))}

#Helper function to calculatethe log likelihood
getLogLike = function(X,y,Ti,beta,u,vu,vb,a0,b0) {
  n = length(X)
  
  loglike = 0
  
  #If vu is negative, set log likelihood to minus infinity
  if (vu<0) return(-Inf)
  
  #Pre-calculating the Bernoulli part since the calculation is complicated
  bernpart = 0
  for (i in 1:n) {
    for (t in 1:Ti[i]) {
      bernpart = bernpart + y[[i]][t]*log(i.logit( X[[i]][t,]%*%beta+u[i] )) +
          (1-y[[i]][t])*log( 1/(1+exp(X[[i]][t,]%*%beta+u[i])) )
    }
  }
  
  loglike = bernpart - 1/(2*vb)*sum(beta^2) - n/2*log(vu) - 1/(2*vu)*sum(u^2) - (a0+1)*log(vu) - b0/vu
                        
  
  return(loglike)
}

