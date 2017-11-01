#Testing bayesian fused ridge regression function
source("~/Documents/mcgill/math680/assignment5/bayesFused.R")
source("~/Documents/mcgill/math680/assignment5/fusedRidge_crossval.R")

#Quick test to make sure the function runs (not part of the assignment)

#set.seed(680)
#n=50; p=3; sigma=1
#X=matrix(rnorm(n*p), nrow=n, ncol=p)
#beta.star=runif(p,-5,5)
#y=X%*%beta.star + sigma*rnorm(n)

#br_fit = bayesFused(X,y,10000,0.5,1,1,1,1)

##########################################################
##########################################################
#Assignment Question 2(b)(ii)
#Generating data
set.seed(9387)
n = 100; p = 5; v.star = 2
beta.star = c(1,1,0.75,0.5,0.5)
Sigma = matrix(0,p,p)
for (i in 1:p) {
  for (j in 1:p) {
    Sigma[i,j] = 0.7^abs(i-j)
  }
}
#Sampling from normal distribution using inverse prob transform, and eigendecomp
Z1=matrix(qnorm(runif(n*p)), nrow=n, ncol=p)
eo1=eigen(Sigma, symmetric=TRUE)
Sigma.sqrt=eo1$vec %*% diag(eo1$val^0.5)%*%t(eo1$vec)
X = matrix(0,n,p) + Z1%*%Sigma.sqrt

y = qnorm(runif(n), X%*%beta.star, v.star)

#Doing 5-fold cross validation
lam = 10^(seq(-8,8,0.5))
cv_fit = cvfr(cbind(0,X), y, lam, lam, 5)

#Now doing Bayesian version
w = cv_fit$best.lam1/(cv_fit$best.lam1+cv_fit$best.lam2)
al = 2
bl = 1/(cv_fit$best.lam1+cv_fit$best.lam2)
av = 2
bv = 3/100*crossprod(X%*%cv_fit$b-y)

bayes_fit = bayesFused(X,y,10000,w,av,bv,al,bl)

save(bayes_fit, file="~/Documents/mcgill/math680/assignment5/dat_a5_q2bii.RData")

#99% credible intervals for beta
bcred = matrix(0,5,2); colnames(bcred) = c("lower","upper")
bcred[1,] = quantile(bayes_fit$beta.post[,1], c(0.005,0.995))
bcred[2,] = quantile(bayes_fit$beta.post[,2], c(0.005,0.995))
bcred[3,] = quantile(bayes_fit$beta.post[,3], c(0.005,0.995))
bcred[4,] = quantile(bayes_fit$beta.post[,4], c(0.005,0.995))
bcred[5,] = quantile(bayes_fit$beta.post[,5], c(0.005,0.995))

#Estimates and credible intervals for beta
beta_est = cbind(beta.star, colMeans(bayes_fit$beta.post),bcred); colnames(beta_est) = c("Truth","Estimate","Lower","Upper")
library(xtable)
xtable(beta_est, digits = 4)

#Autocorrelation plot
plot(acf(rowMeans(bayes_fit$beta.post)), main="Autocorrelation of mean of beta")

#Plotting trace
par(mfrow=c(2,2), mar=rep(2,4))
plot(cumsum(bayes_fit$beta.post[,1])/(1:length(bayes_fit$beta.post[,1])), type="l")
plot(cumsum(bayes_fit$beta.post[,2])/(1:length(bayes_fit$beta.post[,2])), type="l")
plot(cumsum(bayes_fit$beta.post[,3])/(1:length(bayes_fit$beta.post[,3])), type="l")
plot(cumsum(bayes_fit$beta.post[,4])/(1:length(bayes_fit$beta.post[,4])), type="l")
plot(cumsum(bayes_fit$beta.post[,5])/(1:length(bayes_fit$beta.post[,5])), type="l")

