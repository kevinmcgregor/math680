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
  
  #Do a reduced singular value decomposition
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


#Testing out the function (borrowed from Notes 2)
set.seed(680)
n=50; p=10; lambda=1; sigma=1
X=cbind(1, matrix(rnorm(n*p), nrow=n, ncol=p))
Xtilde = scale(X[,-1], scale=FALSE)
beta.star=c(1,runif(p,-5,5))
y=X%*%beta.star + sigma*rnorm(n)
centeredRidge(X,y,lambda)

#Trying replications and comparing to theoretical variance of estimator
nrep = 1000
betarep = matrix(0,p,nrep)
ycur = NULL
for (i in 1:nrep) {
  ycur=X%*%beta.star + sigma*rnorm(n)
  betarep[,i] = centeredRidge(X,ycur,lambda)$b
}

betaexpect = apply(betarep,1,mean)
theor_betaexpect = qr.solve(crossprod(Xtilde, Xtilde)+lambda*diag(p))%*%crossprod(Xtilde, Xtilde)%*%beta.star[-1]

betavar = (n-1)/n * cov(t(betarep))
theor_betavar = sigma*qr.solve(crossprod(Xtilde, Xtilde)+lambda*diag(p))%*%t(Xtilde)%*%
                    t(qr.solve(crossprod(Xtilde, Xtilde)+lambda*diag(p))%*%t(Xtilde))
