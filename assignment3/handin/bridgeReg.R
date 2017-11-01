# Function to do bridge regression.  Uses MM algorithm

#A function to do bridge regression
# y is a vector for the response
# X is a matrix for the 
bridgeReg = function(y,X,lam=1,alpha=1.5,tol=1e-10,maxit=1000) {
  if (!is.matrix(X)) stop("X must be a matrix")
  if (!(is.vector(y) | is.matrix(y)) | !is.numeric(y)) stop("y must be a numeric vector or matrix")
  if (!is.numeric(lam) | lam <= 0) stop("lam must be positive numeric")
  if (alpha <= 1 | alpha >=2) stop("alpha must be in (1,2), exclusive")
  if (NROW(X)!=NROW(y)) stop("Number of observations differs between X and y")
  
  y = as.vector(y)
  
  #Centering
  Xc = scale(X[,-1], scale = FALSE)
  yc = scale(y, scale = FALSE)
  ybar = mean(y)
  Xbar = colMeans(Xc)
  
  p = NCOL(Xc)
  n = NROW(Xc)
  
  count=1
  beta.cur = rep(0,p)
  beta.new = rep(1,p)
  mvec = rep(0,p)
  while (sum(abs(beta.new-beta.cur)) > tol & count<=maxit) {
    beta.cur = beta.new
    mvec = abs(beta.cur)^(alpha-2)
    M = matrix(0,p,p)
    diag(M) = mvec
    #Mimimizing majorization function
    beta.new = qr.solve(crossprod(Xc)+lam*M)%*%crossprod(Xc,y)
    count=count+1
  }
  
  #Calculating intercept
  b1 = ybar - crossprod(Xbar, beta.new)
  return(list(b1=b1,b=beta.new,total.iterations=count-1))
}
