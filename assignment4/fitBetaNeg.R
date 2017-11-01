#Function to do least squares with the last beta forced to be negative

#X is the design matrix (nxp)
#y is the response (nx1 vector or matrix)
fitBetaNeg = function(X, y) {
  if (!is.matrix(X)) stop("X must be a matrix")
  if (!is.vector(y) & !is.matrix(y)) stop("y must be a vector or matrix")
  if (NROW(X)!=length(y)) stop("Number of observations in X and y differs")
  
  n = NROW(X)
  p = NCOL(X)
  
  #Get the value of u in the KKT conditions
  R = c(rep(0,p-1),1)
  u = 2*t(R)%*%qr.solve(crossprod(X))%*%crossprod(X,y)/(crossprod(R,qr.solve(crossprod(X)))%*%R)
  
  #If u is negative, then the OLS estimate of beta_p is already negative so we don't need to do any
  #constraining.  So just set u=0
  u = ifelse(u<0, 0, u)
  
  #Final constrained estimate of beta
  beta.bar = qr.solve(crossprod(X), crossprod(X,y)-u/2*R)
  
  return(beta.bar)
}

