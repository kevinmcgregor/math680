#Function to to cross-validation for fused ridge regression

cvfr = function(X, y, lam1.vec, lam2.vec, K) {
  n = length(y)
  if (!(K %in% 2:n)) stop("K must be one of 2, 3,..., n")
  if (!is.matrix(X)) stop("X must be a matrix")
  if (!(is.vector(y) | is.matrix(y)) | !is.numeric(y)) stop("y must be a numeric vector or matrix")
  if (!is.vector(lam1.vec) | !is.numeric(lam1.vec)) stop("lam1.vec must be a numeric vector")
  if (!is.vector(lam2.vec) | !is.numeric(lam2.vec)) stop("lam2.vec must be a numeric vector")
  if (any(lam1.vec<=0|lam2.vec<=0)) stop("All elements of lam1.vec and lam2.vec must be positive")
  if (NROW(X)!=NROW(y)) stop("Number of observations differs between X and y")
  
  #create K folds
  groups = makeFolds(n,K)
  
  cv.error = matrix(0, length(lam1.vec), length(lam2.vec))
  tot_err_k = NULL
  train = NULL
  beta = NULL
  count=1
  tot_err_k = matrix(0, length(lam1.vec), length(lam2.vec))
  #loop over different values of lam1 and lam2
  for (i in 1:length(lam1.vec)) {
    for (j in 1:length(lam2.vec)) {
      for (k in 1:K) {
          #fit on training data
          train = fusedRidge(X[groups!=k,], y[groups!=k], lam1.vec[i], lam2.vec[j])
          
          #fit on test data
          beta = c(train$b1, t(train$b))
          tot_err_k[i,j] = tot_err_k[i,j] + sum((y[groups==k]-X[groups==k,]%*%beta)^2)
      }
    }
  }
  
  #Calculate the cross validation error
  cv.error = tot_err_k / n    
  
  #Find lambdas giving minimum error
  ind.min = which(cv.error == min(cv.error), arr.ind = TRUE)
  best.lam1 = lam1.vec[ind.min[1]]
  best.lam2 = lam2.vec[ind.min[2]]

  finalfit = fusedRidge(X, y, best.lam1, best.lam2)
  
  return(list(m=finalfit$b1, b=finalfit$b, best.lam1=best.lam1, best.lam2=best.lam2, cv.error=cv.error))
}

# A function to randomly make K folds.  Returns a vector of size n specifying
# which group each observation was put into.
makeFolds = function(n,K) {
  
  #Divide into K groups of equal size, and any remainder gets randomly sorted into one of the K groups
  init_groupsize = n%/%K
  num_remain = n%%K
  fold = rep(1:K, each=init_groupsize)
  fold = c(fold, sample(1:K, num_remain))
  #fold = NULL
  #for (i in 1:K) {
  #  fold = c(fold, rep(i, init_groupsize))
  #}
  #fold = c(fold, sample(1:K, num_remain))
  
  #Randomly permute the group choices
  perm = sample(1:n)
  
  return(fold[perm])
}


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

