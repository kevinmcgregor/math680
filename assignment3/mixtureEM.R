#EM algorithm for Gaussian mixture model

#Some code borrowed from notes
mixtureEM = function(X, C, tol=1e-8, maxit=100) {
  p=NCOL(X)
  n=NROW(X)
  
  #Inital param values
  pi.list = rep(1/C, C)
  mu.mat=matrix(rnorm(p*C), nrow=p, ncol=C) + apply(X, 2, mean)%*%matrix(1, nrow=1, ncol=C)
  Sigma=vector(length=C, mode="list")
  for(y in 1:C){Sigma[[y]]=diag(p)}
  
  objfnval=-Inf
  
  k=0
  iterating=TRUE
  while(iterating)
  {
    k=k+1
    
    #########################################################################
    ##  Compute the n by C matrix P.mat who's (i,j)th entry 
    ##  is the current iterate's estimate of P(Y=j|X=x_i)
    logPhinum=matrix(NA, nrow=n, ncol=C)
    for (y in 1:C)
    {
      Xcentered=scale(X, scale=FALSE, center=mu.mat[,y])
      qf=-0.5*apply(t(Xcentered)*qr.solve(Sigma[[y]], t(Xcentered)), 2, sum)
      logPhinum[,y]=qf-0.5*determinant(Sigma[[y]], logarithm=TRUE)$mod[1]-0.5*p*log(2*pi)+log(pi.list[y])
    }
    
    newobjfnval=sum(log(apply(exp(logPhinum), 1, sum)))  
    ## control numerical stability by adjusting on the log scale:
    ##   subtract the row maximum from each element in the row 
    logPhinum = logPhinum - apply(logPhinum, 1, max)%*%matrix(1, nrow=1, ncol=C)
    Phinum=exp(logPhinum)
    Phiden=apply(Phinum, 1, sum)
    P.mat=Phinum/Phiden
    ########################################################################
        
    if( ((newobjfnval - objfnval) < tol)  | (k > maxit) )
      iterating=FALSE
    
    objfnval=newobjfnval 
    
    for(y in 1:C)
    {
      ## update the pi's
      sumprob=sum(P.mat[,y])
      pi.list[y]=sumprob/n
      
      ## update the mu's
      weight.matrix=diag(P.mat[,y]/sumprob)
      mu.mat[,y] = apply(weight.matrix %*% X, 2, sum) 
      
      ## update the Sigma's
      Xcentered=scale(X, scale=FALSE, center=mu.mat[,y])
      Sigma[[y]] =  crossprod(Xcentered, weight.matrix%*%Xcentered)
    }
  }
  return(list(mu.mat=mu.mat, Sigma=Sigma, pi.list=pi.list, P.mat=P.mat, total.iterations=k))
  
}


#Gaussian mixture model where we assume sigma is constant over all categories
mixtureEM_csig = function(X, C, tol=1e-8, maxit=100) {
  p=NCOL(X)
  n=NROW(X)
  
  #Inital param values
  pi.list = rep(1/C, C)
  mu.mat=matrix(rnorm(p*C), nrow=p, ncol=C) + apply(X, 2, mean)%*%matrix(1, nrow=1, ncol=C)
  
  #Here sigma is constant over different categories
  Sigma=diag(p)
  
  objfnval=-Inf
  
  k=0
  iterating=TRUE
  while(iterating)
  {
    k=k+1
    
    #########################################################################
    ##  Compute the n by C matrix P.mat who's (i,j)th entry 
    ##  is the current iterate's estimate of P(Y=j|X=x_i)
    logPhinum=matrix(NA, nrow=n, ncol=C)
    for (y in 1:C)
    {
      Xcentered=scale(X, scale=FALSE, center=mu.mat[,y])
      qf=-0.5*apply(t(Xcentered)*qr.solve(Sigma, t(Xcentered)), 2, sum)
      logPhinum[,y]=qf-0.5*determinant(Sigma, logarithm=TRUE)$mod[1]-0.5*p*log(2*pi)+log(pi.list[y])
    }
    
    newobjfnval=sum(log(apply(exp(logPhinum), 1, sum)))  
    ## control numerical stability by adjusting on the log scale:
    ##   subtract the row maximum from each element in the row 
    logPhinum = logPhinum - apply(logPhinum, 1, max)%*%matrix(1, nrow=1, ncol=C)
    Phinum=exp(logPhinum)
    Phiden=apply(Phinum, 1, sum)
    P.mat=Phinum/Phiden
    ########################################################################
    
    if( ((newobjfnval - objfnval) < tol)  | (k > maxit) )
      iterating=FALSE
    
    objfnval=newobjfnval 
    
    
    Sigsum = matrix(0,p,p)
    sumprobsum = 0
    for(y in 1:C)
    {
      ## update the pi's
      sumprob=sum(P.mat[,y])
      pi.list[y]=sumprob/n
      
      ## update the mu's
      weight.matrix=diag(P.mat[,y]/sumprob)
      mu.mat[,y] = apply(weight.matrix %*% X, 2, sum) 
      
      ## update Sigma
      Xcentered=scale(X, scale=FALSE, center=mu.mat[,y])
      sig.weight.matrix = diag(P.mat[,y])
      Sigsum = Sigsum + crossprod(Xcentered, sig.weight.matrix%*%Xcentered)
      sumprobsum = sumprobsum + sumprob
    }
    
    Sigma = Sigsum/sumprobsum
  }
  return(list(mu.mat=mu.mat, Sigma=Sigma, pi.list=pi.list, P.mat=P.mat, total.iterations=k))
  
}



