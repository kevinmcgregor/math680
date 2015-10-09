ridgeSim = function(X, beta, lam.vec, sigma.star, numrep=200) {
  
  n = NROW(x)
  #Calculating true value of X times Beta
  Xbeta = X%*%beta
  
  y=rep(0,n)
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

