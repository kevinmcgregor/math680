#Function to do bridge-penalized logistic regression using MM algorithm

bridgeLogistic = function(y,X,n.list,lam=1,alpha=1.5,tol=1e-8,maxit=100) {
  if (!is.matrix(X)) stop("X must be a matrix")
  if (!(is.vector(y) | is.matrix(y)) | !is.numeric(y)) stop("y must be a numeric vector or matrix")
  if (!is.numeric(lam) | lam <= 0) stop("lam must be positive numeric")
  if (!is.numeric(alpha) | alpha <= 1 | alpha >=2) stop("alpha must be numeric and in (1,2), exclusive")
  if (NROW(X)!=NROW(y)) stop("Number of observations differs between X and y")
  
  
  y = as.vector(y)
  
  p = ncol(X)
  n = length(y)
  
  count=1
  beta.cur = rep(1,p)
  beta.new = rep(0,p)
  while (sum(abs(beta.new-beta.cur)) > tol & count<=maxit) {
    beta.cur = beta.new
    mvec = abs(beta.cur)^(alpha-2)
    M = matrix(0,p,p)
    diag(M) = mvec
    
    
    
    count=count+1
  }
  
  return(list(b=beta.new, total.iterations=count-1))
}





