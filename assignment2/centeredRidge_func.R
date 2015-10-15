#Function to do centered ridge regression

# A function to do centered ridge regression for a given value of lambda
# X is the design matrix (must be a matrix)
# y is the response vector (must be a vector) 
# lam is a vector containing values of lambda (must be numeric and positive)
# SVD is an optional variable to take in a predefined (full) singular value decomposition to save time
#   in a simulation with multiple replications
# n.orig is optional, and contains the original sample size (which would be reduced in cross validation)
# returns a list with two elements.  The first element is a vector of intercepts for each lambda,
# and the second element is a slopes matrix with rows corresponding to different values of lambda, and
# columns corresponding to the columns of X
centeredRidge = function(X,y,lam,SVD=NULL,n.orig=NULL) {
  if (!is.matrix(X)) stop("X must be a matrix")
  if (!(is.vector(y) | is.matrix(y)) | !is.numeric(y)) stop("y must be a numeric vector or matrix")
  if (!is.vector(lam) | !is.numeric(lam)) stop("lam must be a numeric vector")
  if (NROW(X)!=NROW(y)) stop("Number of observations differs between X and y")

  
  y=as.vector(y)
  
  #Centering X and y
  Xc = scale(X[,-1], scale = FALSE)
  yc = scale(y, scale = FALSE)
  ybar = mean(y)
  Xbar = colMeans(Xc) #apply(Xc,2,mean)
  
  p = NCOL(Xc)
  n = length(y)
  q = min(n-1,p)
  
  #Do a singular value decomposition (or use predefined one if it was included as a parameter)
  if (is.null(SVD)) {
    sv = svd(Xc, nu=n, nv=p)
    n.orig = n
  } else {
    sv = SVD
  }
  
  djj = sv$d[1:q]
  u=sv$u; v=sv$v
  
  #Add on respective values of lambda
  nlam = length(lam)
  M = vector("list", nlam)
  b1 = rep(0, nlam)
  b = matrix(0, p, nlam)
  for (i in 1:nlam) {
    M[[i]] = matrix(0,p,n.orig)
    diag(M[[i]])[1:q] = djj/(djj^2+lam[i])
    #calculate slopes and intercept
    b[,i] = v%*%M[[i]]%*%crossprod(u,y)  #t(u)%*%y
    b1[i] = ybar - crossprod(Xbar,b[,i])
  }
  
  return(list(b1=b1,b=b))
}