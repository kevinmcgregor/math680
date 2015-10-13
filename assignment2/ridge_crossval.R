source("~/Documents/mcgill/math680/assignment2/centeredRidge_func.R")

# Funtion to do cross validation to get optimal lambda in ridge regression
# X is the design matrix
# y is the reponse vector
# lam.vec is a vector of lambda values to iterate over
# K is the number of folds in the cross validation
ridgeCrossval = function(X, y, lam.vec, K) {
  n = length(y)
  if (!(K %in% 2:n)) stop("K must be one of 2, 3,..., n")
  if (!is.matrix(X)) stop("X must be a matrix")
  if (!(is.vector(y) | is.matrix(y)) | !is.numeric(y)) stop("y must be a numeric vector or matrix")
  if (!is.vector(lam.vec) | !is.numeric(lam.vec)) stop("lam.vec must be a numeric vector")
  if (any(lam.vec<=0)) stop("All elements of lam.vec must be positive")
  if (NROW(X)!=NROW(y)) stop("Number of observations differs between X and y")
  
  groups = makeFolds(n,K)
  cv.error = rep(0, length(lam.vec))
  tot_err_k = NULL
  train = NULL
  beta = NULL
  count=1
  tot_err_k = rep(0, length(lam.vec))
  for (j in 1:K) {
    train = centeredRidge(X[groups!=j,], y[groups!=j], lam.vec)
    beta = cbind(train$b1, t(train$b))
    for (i in 1:length(lam.vec)) {
      tot_err_k[i] = tot_err_k[i] + sum((y[groups==j]-X[groups==j,]%*%beta[i,])^2)
    }
  }
  
  cv.error = tot_err_k / n    
  
  best.lam = lam.vec[which.min(cv.error)]
  finalfit = centeredRidge(X, y, best.lam)
  
  return(list(b1=finalfit$b1, b=finalfit$b, best.lam=best.lam, cv.error=cv.error))
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


#Testing
#set.seed(342)
#n=500; p=40; sigma=2
#X=cbind(1, matrix(rnorm(n*p), nrow=n, ncol=p))
#Xtilde = scale(X[,-1], scale=FALSE)
#beta.star=c(1,runif(p,-5,5))
#y=X%*%beta.star + sigma*rnorm(n)
#lamvec = seq(0.001,5,0.05)

#tmp = ridgeCrossval(X, y, lamvec, K=5)
#plot(lamvec, tmp$cv.error, type = "l")



