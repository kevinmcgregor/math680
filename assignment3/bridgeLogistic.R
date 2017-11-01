#Function to do bridge penalized logisitic regression

bridgeLogistic = function(X,y,n.list,lam=1,alpha=1.5,tol=1e-8,maxit=100) {
  if (!is.matrix(X)) stop("X must be a matrix")
  if (!(is.vector(y) | is.matrix(y)) | !is.numeric(y)) stop("y must be a numeric vector or matrix")
  if (!is.numeric(lam) | lam <= 0) stop("lam must be positive numeric")
  if (!is.numeric(alpha) | alpha <= 1 | alpha >=2) stop("alpha must be numeric and in (1,2), exclusive")
  if (NROW(X)!=NROW(y)) stop("Number of observations differs between X and y")
  
  y = as.vector(y)
  
  
  p = ncol(X)
  n = length(y)
  
  count=1
  beta.cur = rep(0,p)
  beta.new = rep(0.5,p)
  while (sum(abs(beta.new-beta.cur)) > tol & count<=maxit) {
    beta.cur = beta.new
  
    m = abs(beta.cur)^(alpha-2)
    lam.m = lam*m
    lam.M = diag(c(0,lam.m[-1]))
    X.t.n.list.y=crossprod(X, n.list*y)
    pi.t=ilogit(as.numeric(X%*%beta.cur))
    #W=diag(n.list*pi.t*(1-pi.t))
    #W = diag(max(n.list)*rep(1/4, n))
    
    LHS = max(n.list)*1/4*crossprod(X,X) + lam.M
    RHS = X.t.n.list.y-crossprod(X, n.list*pi.t)
    
    beta.new = qr.solve(LHS, RHS)
    
    count=count+1
  }
  
  cat(sum(abs(beta.new-beta.cur)))
  
  return(list(b=beta.new, total.iterations=count-1))
}

#ilogit=function(u) { require(boot); return(inv.logit(u)) }
ilogit=function(u) { return(exp(u)/(1+exp(u))) }

#Minimizes the current majorization function using Newton's method
min_major = function(X, y, n.list, lam, alpha, beta) {
  n = NROW(X)
  
  maj.tol = 1e-5
  maj.maxit = 10
  
  #M = matrix(0,p,p)
  #diag(M) = abs(beta)^(alpha-2)
  m = abs(beta)^(alpha-2)
  #m=c(0, rep(1, p-1))
  lam.m = lam*m
  lam.M = diag(lam.m)
  X.t.n.list.y=crossprod(X, n.list*y)
  
  add=1
  count = 1
  mb.cur = rep(1,p)
  mb.new = rep(0,p)
  while (sum(abs(add)) > maj.tol & count<=maj.maxit) {
    mb.cur = mb.new
    #eta = ilogit(X%*%mb.cur)
    #W = matrix(0,n,n)
    #diag(W) = n.list*eta*(1-eta)
  
    #gradient = crossprod(X, n.list*(eta-y)) #+ lam*M%*%mb.cur
    #hessian = crossprod(X,W%*%X) #+ lam*M
    
    pi.t=ilogit(as.numeric(X%*%mb.cur))
    W=diag(n.list*pi.t*(1-pi.t))
    
    
    minusGrad=X.t.n.list.y-crossprod(X, n.list*pi.t) - lam.m*mb.cur
    Hess=crossprod(X,W%*%X)+lam.M
    add=qr.solve(Hess, minusGrad)
    mb.new=mb.cur+add
    
    #mb.new = mb.cur + qr.solve(hessian, -gradient)
    
    count=count+1
  }
  
  return(mb.new)
}


#function to calculate the gradient
grad = function(X, y, n.list, lam, alpha, beta) {
  eta = ilogit(X%*%beta)
  return( crossprod(X, n.list*(eta-y)) + lam*sign(beta)*abs(beta)^(alpha-1) )
}

#function to calculate the Hessian matrix
hess = function(X, y, n.list, lam, alpha, beta) {
  n = NROW(X)
  p = NCOL(X)
  eta = ilogit(X%*%beta)
  W = matrix(0,n,n)
  diag(W) = n.list*eta*(1-eta)
  M = matrix(0,p,p)
  diag(M) = abs(beta)^(alpha-2)
  
  return( crossprod(X,W%*%X) + lam*(alpha-1)*M )
}




