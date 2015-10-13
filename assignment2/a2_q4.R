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




ridgeSim = function(X, beta, lam.vec, sigma.star, numrep=200, n.cores=1) {
  require(parallel)
  
  n = NROW(X)
  Xbeta = X%*%beta
  y = replicate(numrep, Xbeta + rnorm(n, 0, sqrt(sigma.star)), simplify=FALSE)
  
  losses = mclapply(y, FUN=getLosses, Xmat=X, beta=beta, Xbeta=Xbeta, lam.vec=lam.vec, sigma.star=sigma.star, mc.cores = n.cores)
  
  loss_ols_b = sapply(losses, function(x){x[1]})
  loss_k5_b = sapply(losses, function(x){x[2]})
  loss_k10_b = sapply(losses, function(x){x[3]})
  loss_kn_b = sapply(losses, function(x){x[4]})
  loss_ols_xb = sapply(losses, function(x){x[5]})
  loss_k5_xb = sapply(losses, function(x){x[6]})
  loss_k10_xb = sapply(losses, function(x){x[7]})
  loss_kn_xb = sapply(losses, function(x){x[8]})
  
  #Means and standard errors of the different loss functions
  avg_loss = c(mean(loss_ols_b),mean(loss_k5_b),mean(loss_k10_b),mean(loss_kn_b),
               mean(loss_ols_xb), mean(loss_k5_xb), mean(loss_k10_xb), mean(loss_kn_xb))
  se_loss = c(sd(loss_ols_b),sd(loss_k5_b),sd(loss_k10_b),sd(loss_kn_b),
              sd(loss_ols_xb), sd(loss_k5_xb), sd(loss_k10_xb), sd(loss_kn_xb))
  
  final = data.frame(avg_loss=avg_loss, se_loss=se_loss)
  rownames(final) = c("OLS_B","K5_B","K10_B","Kn_B","OLS_XB","K5_XB","K10_XB","Kn_XB")
  
  return(final)
}

# A function to do the fitting and calculate losses. Main function to be included in parallel step  
getLosses = function(y, Xmat, beta, Xbeta, lam.vec, sigma.star) {
  
  y=unlist(y)
  n = NROW(Xmat)
  
  
  k5 = ridgeCrossval(Xmat, y, lam.vec, 5)
  k10 = ridgeCrossval(Xmat, y, lam.vec, 10)
  kn = ridgeCrossval(Xmat, y, lam.vec, n)
  
  beta_ols = lm.fit(Xmat, y)$coefficients
  beta_k5 = c(k5$b1, k5$b)
  beta_k10 = c(k10$b1, k10$b)
  beta_kn = c(kn$b1, kn$b)
  
  loss_ols_b = sum((beta_ols-beta)^2)
  loss_k5_b = sum((beta_k5-beta)^2)
  loss_k10_b = sum((beta_k10-beta)^2)
  loss_kn_b = sum((beta_kn-beta)^2)
  
  loss_ols_xb = sum((Xmat%*%beta_ols-Xbeta)^2)
  loss_k5_xb = sum((Xmat%*%beta_k5-Xbeta)^2)
  loss_k10_xb = sum((Xmat%*%beta_k10-Xbeta)^2)
  loss_kn_xb = sum((Xmat%*%beta_kn-Xbeta)^2)
  
  return(c(loss_ols_b, loss_k5_b, loss_k10_b, loss_kn_b, loss_ols_xb, loss_k5_xb, loss_k10_xb, loss_kn_xb))
}


#4:51
lambda = 10^(seq(-8,8,0.5))
test = ridgeSim(X1,beta1,lambda,sigma.star,n.cores=4)
test2 = ridgeSim(X2,beta1,lambda,sigma.star,n.cores=4)
test3 = ridgeSim(X3,beta2,lambda,sigma.star,n.cores=4)
test4 = ridgeSim(X4,beta2,lambda,sigma.star,n.cores=4)

save(test, file="~/Documents/mcgill/math680/assignment2/dat4a.RData")
save(test2, file="~/Documents/mcgill/math680/assignment2/dat4a2.RData")
save(test3, file="~/Documents/mcgill/math680/assignment2/dat4a3.RData")
save(test4, file="~/Documents/mcgill/math680/assignment2/dat4a4.RData")

#ty = X1%*%beta1 + rnorm(100,0,0.5)
#tyfit = lm.fit(X1, ty)


