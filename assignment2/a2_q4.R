# Doing simulations for assignment 2 question 4
source("~/Documents/mcgill/math680/assignment2/ridge_crossval.R")

#Setting simulation parameters
n=100; p1=50; p2=1000; theta1=0.5; theta2=0.9

sigma.star = 0.5

beta1 = rnorm(p1, 0, sqrt(sigma.star))
beta2 = rnorm(p2, 0, sqrt(sigma.star))

#Covariance matrices for generating X
Sigma1 = matrix(0, p1-1, p1-1)
for (i in 1:(p1-1)) {
  for (j in 1:(p1-1)) {
    Sigma1[i,j] = theta1^abs(i-j)
  }
}

Sigma2 = matrix(0, p1-1, p1-1)
for (i in 1:(p1-1)) {
  for (j in 1:(p1-1)) {
    Sigma2[i,j] = theta2^abs(i-j)
  }
}

Sigma3 = matrix(0, p2-1, p2-1)
for (i in 1:(p2-1)) {
  for (j in 1:(p2-1)) {
    Sigma3[i,j] = theta1^abs(i-j)
  }
}

Sigma4 = matrix(0, p2-1, p2-1)
for (i in 1:(p2-1)) {
  for (j in 1:(p2-1)) {
    Sigma4[i,j] = theta2^abs(i-j)
  }
}

#Use covariance matrices in multivar normal dist to create X
Z1=matrix(rnorm(n*(p1-1)), nrow=n, ncol=p1-1)
eo1=eigen(Sigma1, symmetric=TRUE)
Sigma1.sqrt=eo1$vec %*% diag(eo1$val^0.5)%*%t(eo1$vec)
Xtilde1 = matrix(0,n,p1-1) + Z1%*%Sigma1.sqrt
X1 = cbind(1,Xtilde1)

Z2=matrix(rnorm(n*(p1-1)), nrow=n, ncol=p1-1)
eo2=eigen(Sigma2, symmetric=TRUE)
Sigma2.sqrt=eo2$vec %*% diag(eo2$val^0.5)%*%t(eo2$vec)
Xtilde2 = matrix(0,n,p1-1) + Z2%*%Sigma2.sqrt
X2 = cbind(1,Xtilde2)

Z3=matrix(rnorm(n*(p2-1)), nrow=n, ncol=p2-1)
eo3=eigen(Sigma3, symmetric=TRUE)
Sigma3.sqrt=eo3$vec %*% diag(eo3$val^0.5)%*%t(eo3$vec)
Xtilde3 = matrix(0,n,p2-1) + Z3%*%Sigma3.sqrt
X3 = cbind(1,Xtilde3)

Z4=matrix(rnorm(n*(p2-1)), nrow=n, ncol=p2-1)
eo4=eigen(Sigma4, symmetric=TRUE)
Sigma4.sqrt=eo4$vec %*% diag(eo4$val^0.5)%*%t(eo4$vec)
Xtilde4 = matrix(0,n,p2-1) + Z4%*%Sigma4.sqrt
X4 = cbind(1,Xtilde4)




ridgeSim = function(X, beta, lam.vec, sigma.star, numrep=200) {

  n = NROW(x)
  
  k5=k10=kn=NULL
  beta_ols=beta_k5=beta_k10=beta_kn=NULL
  loss_ols_b=loss_ols_xb=rep(0,numrep)
  loss_k5_b=loss_k5_xb=rep(0,numrep)
  loss_k10_b=loss_k10_xb=rep(0,numrep)
  loss_kn_b=loss_kn_xb=rep(0,numrep)
  for (i in 1:numrep) {
    #Generate y
    y = Xbeta + rnorm(n, 0, sqrt(sigma.star))
    
    k5 = ridgeCrossval(X, y, lam.vec, 5)
    k10 = ridgeCrossval(X, y, lam.vec, 10)
    kn = ridgeCrossval(X, y, lam.vec, n)
    
    beta_ols = lm.fit(X, y)$coefficients
    beta_k5 = c(k5$b1, k5$b)
    beta_k10 = c(k10$b1, k10$b)
    beta_kn = c(kn$b1, kn$b)
    
    loss_ols_b[i] = sum((beta_ols-beta)^2)
    loss_k5_b[i] = sum((beta_k5-beta)^2)
    loss_k10_b[i] = sum((beta_k10-beta)^2)
    loss_kn_b[i] = sum((beta_kn-beta)^2)
    
    loss_ols_xb[i] = sum((X%*%beta_ols-Xbeta)^2)
    loss_k5_xb[i] = sum((X%*%beta_k5-Xbeta)^2)
    loss_k10_xb[i] = sum((X%*%beta_k10-Xbeta)^2)
    loss_kn_xb[i] = sum((X%*%beta_kn-Xbeta)^2)
    
  }
  
  #Means and standard errors of the different loss functions
  avg_loss = c(mean(loss_ols_b),mean(loss_k5_b),mean(loss_k10_b),mean(loss_kn_b),
               mean(loss_ols_xb), mean(loss_k5_xb), mean(loss_k10_xb), mean(loss_kn_xb))
  se_loss = c(sd(loss_ols_b),sd(loss_k5_b),sd(loss_k10_b),sd(loss_kn_b),
              sd(loss_ols_xb), sd(loss_k5_xb), sd(loss_k10_xb), sd(loss_kn_xb))
  
  final = data.frame(avg_loss=avg_loss, se_loss=se_loss)
  rownames(final) = c("OLS_B","K5_B","K10_B","Kn_B","OLS_XB","K5_XB","K10_XB","Kn_XB")
  
  return(list(beta_ols,final))
}

# A function to generate y and calculate losses. Main function to be included in parallel step  
getLosses = function(X, beta, lam.vec, sigma.star) {
  
  Xbeta = X%*%beta
  #Generate y
  y = Xbeta + rnorm(n, 0, sqrt(sigma.star))
  
  k5 = ridgeCrossval(X, y, lam.vec, 5)
  k10 = ridgeCrossval(X, y, lam.vec, 10)
  kn = ridgeCrossval(X, y, lam.vec, n)
  
  beta_ols = lm.fit(X, y)$coefficients
  beta_k5 = c(k5$b1, k5$b)
  beta_k10 = c(k10$b1, k10$b)
  beta_kn = c(kn$b1, kn$b)
  
  loss_ols_b = sum((beta_ols-beta)^2)
  loss_k5_b = sum((beta_k5-beta)^2)
  loss_k10_b = sum((beta_k10-beta)^2)
  loss_kn_b = sum((beta_kn-beta)^2)
  
  loss_ols_xb = sum((X%*%beta_ols-Xbeta)^2)
  loss_k5_xb = sum((X%*%beta_k5-Xbeta)^2)
  loss_k10_xb = sum((X%*%beta_k10-Xbeta)^2)
  loss_kn_xb = sum((X%*%beta_kn-Xbeta)^2)
  
  return(c(loss_ols_b, loss_k5_b, loss_k10_b, loss_kn_b, loss_ols_xb, loss_k5_xb, loss_k10_xb, loss_kn_xb))
}



lambda = 10^(seq(-8,8,0.5))
test = ridgeSim(X1,beta1,lambda,sigma.star)
test2 = ridgeSim(X2,beta1,lambda,sigma.star)
test3 = ridgeSim(X3,beta2,lambda,sigma.star)
test4 = ridgeSim(X4,beta2,lambda,sigma.star)

save(test, file="~/Documents/mcgill/math680/assignment2/dat4a.RData")
save(test2, file="~/Documents/mcgill/math680/assignment2/dat4a2.RData")
save(test3, file="~/Documents/mcgill/math680/assignment2/dat4a3.RData")
save(test4, file="~/Documents/mcgill/math680/assignment2/dat4a4.RData")

#ty = X1%*%beta1 + rnorm(100,0,0.5)
#tyfit = lm.fit(X1, ty)


