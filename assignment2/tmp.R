getLosses = function(y, X, beta, Xbeta, lam.vec, sigma.star) {
  
  y=unlist(y)
  n = NROW(X)
  
  cat(y)
  
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

