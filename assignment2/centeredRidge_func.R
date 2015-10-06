#Function to do centered ridge regression

# A function to do centered ridge regression for a given value of lambda
# X is the design matrix (must be a matrix)
# y is the response vector (must be a vector) 
# lam is lambda (must be numeric and positive)
centeredRidge = function(X,y,lam) {
  if (!is.matrix(X)) stop("X must be a matrix")
  if (!(is.vector(y) | is.matrix(y)) | !is.numeric(y)) stop("y must be a numeric vector or matrix")
  if (lam<=0 | !is.numeric(lam)) stop("lam must be numeric and positive")
  if (NROW(X)!=NROW(y)) stop("Number of observations differs between X and y")
  
  y=as.vector(y)
  
  #Centering X and y
  Xc = scale(X[,-1], scale = FALSE)
  yc = scale(y, scale = FALSE)
  ybar = mean(y)
  Xbar = apply(Xc,2,mean)
  
  p = NCOL(Xc)
  n = length(y)
  b1 = 0
  b = rep(0,p-1)
  q = min(n-1,p)
  
  #Do a singular value decomposition
  sv = svd(Xc, nu=n, nv=p)
  djj = sv$d[1:q]
  M = matrix(0,p,n)
  diag(M)[1:q] = djj/(djj^2+lam)
  u=sv$u; v=sv$v
  
  #calculate slopes and intercept
  b = v%*%M%*%t(u)%*%y
  b1 = ybar - t(Xbar)%*%b
  
  return(list(b1=b1,b=b))
}