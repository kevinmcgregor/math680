#Function to generate samples from Bayesian fused ridge regression (using Gibbs sampling)

bayesFused = function(X,y,nsamp,w,av,bv,al,bl) {
  require(MCMCpack)
  require(MASS)
  
  n = NROW(X)
  p = NCOL(X)
  
  #Building the matrix to contrast coefficients
  H = matrix(0,p,p)
  for (i in 1:p) {
    if (i==1) {
      H[i,] = c(1, -1, rep(0, p-2))
    } else if (i==p) {
      H[i,] = c(rep(0, p-2), -1, 1)
    } else {
      H[i,] = c(rep(0, i-2), -1, 2, -1, rep(0, p-i-1))
    }
  }
  
  #Initialization
  beta = rep(1,p)
  v = 1
  lambda = 1
  
  #Storage for posterior samples
  beta.post = matrix(0, nsamp, p)
  v.post = rep(0, nsamp)
  lambda.post = rep(0, nsamp)
  
  for (k in 1:nsamp) {
    #Mean of B posterior
    mu.b = qr.solve(crossprod(X)+lambda*w*diag(p)+lambda*(1-w)*H, crossprod(X,y))
    
    #Sample v given lambda and y
    v.post[k] = rinvgamma(1, av+n/2, bv + 0.5*(crossprod(y)-crossprod(mu.b,(crossprod(X)+lambda*w*diag(p)+
                                                                              lambda*(1-w)*H))%*%mu.b))
    v = v.post[k]
    
    #Sample beta given lambda, v, and y
    beta.post[k,] = mvrnorm(1, mu.b, v*qr.solve(crossprod(X)+lambda*w*diag(p)+lambda*(1-w)*H))
    beta = beta.post[k,]
    
    #Sample lambda given beta, v, and y
    lambda.post[k] = rgamma(1, p/2+al, 1/(2*v)*crossprod(beta,w*diag(p)+(1-w)*H)%*%beta + bl)
    
  }
  
  return(list(beta.post=beta.post, v.post=v.post, lambda.post=lambda.post))
}





