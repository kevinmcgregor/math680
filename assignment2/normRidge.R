#A function to do normal likelihood ridge estimator.


normRidge = function(X,y,lam, tol=1e-12) {
  
  Xc = scale(X[,-1], scale = FALSE)
  yc = scale(y, scale = FALSE)
  ybar = mean(y)
  Xbar = colMeans(Xc)
  
  p = NCOL(Xc)
  n = NROW(Xc)
  
  beta.new = rep(1,p)
  beta.cur = rep(0,p)
  sigma2.new = 1
  sigma2.cur = 0
  while (sum((beta.new-beta.cur)^2) > tol & (sigma2.new-sigma2.cur)^2 > tol) {
    beta.cur=beta.new; sigma2.cur=sigma2.new
    #Update beta holding sigma constant
    beta.new = beta.cur + changeBeta(beta.cur,Xc,yc,lam,sigma2.cur)
    #Update sigma^2 holding beta constant at new beta
    sigma2.new = as.numeric(sigma2.cur + changeSigma(beta.new,Xc,yc,lam,sigma2.cur))
  }
  
  #UPDATE INTERCEPT!!!
  b1 = ybar - crossprod(Xbar, beta.new)
  
  return(list(b1=b1, b=beta.new, sigma2=sigma2.new))
}


#Change in beta according to Newton's method
changeBeta = function(beta, X, y, lam, sigma2) {
  grad = (crossprod(X)%*%beta - crossprod(X, y))/sigma2 + lambda*beta
  hess = crossprod(X)/sigma2 + lambda
  
  return(-1*solve(hess)%*%grad)
}

#Change in sigma^2 according to Newton's method
changeSigma = function(beta, X, y, lam, sigma2) {
  n = length(y)
  grad = n/sigma2 - (crossprod(y)-2*crossprod(X%*%beta,y)+crossprod(X%*%beta,X)%*%beta)/(sigma2)^2
  hess = -1*n/(sigma2)^2 + 2*(crossprod(y)-2*crossprod(X%*%beta,y)+crossprod(X%*%beta,X)%*%beta)/(sigma2^3)
  
  return(-1*solve(hess)%*%grad)
}

