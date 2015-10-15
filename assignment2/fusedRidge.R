#Function to do fused ridge regression

fusedRidge = function(X, y, lam1, lam2) {
  if (!is.matrix(X)) stop("X must be a matrix")
  if (!(is.vector(y) | is.matrix(y)) | !is.numeric(y)) stop("y must be a numeric vector or matrix")
  if (!is.numeric(lam1) | !is.numeric(lam2) | lam1 <= 0 | lam2 <= 0) stop("lam1 and lam2 must be positive numeric")
  if (NROW(X)!=NROW(y)) stop("Number of observations differs between X and y")
  
  y = as.vector(y)
  
  Xc = scale(X[,-1], scale=FALSE)
  yc = scale(y, scale=FALSE)
  ybar = mean(y)
  xbar = colMeans(Xc)
  
  p = NCOL(Xc)
  
  #Building the matrix to contrast coefficients
  M = matrix(0,p,p)
  for (i in 1:p) {
    if (i==1) {
      M[i,] = c(1, -1, rep(0, p-2))
    } else if (i==p) {
      M[i,] = c(rep(0, p-2), -1, 1)
    } else {
      M[i,] = c(rep(0, i-2), -1, 2, -1, rep(0, p-i-1))
    }
  }
  
  #Return closed-form estimates for fused-ridge regression
  b = qr.solve(crossprod(Xc)+lam1*diag(p)+lam2*M)%*%crossprod(Xc,yc)
  b1 = ybar - crossprod(xbar, b)
  
  return(list(b1=b1, b=b))
  
}

