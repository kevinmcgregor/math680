#Testing out and running simulation on, fusedRidge

source("~/Documents/mcgill/math680/assignment2/fusedRidge.R")

set.seed(680)
n=50; p=10; lam1=1; lam2=1; sigma=1
X=cbind(1, matrix(rnorm(n*p), nrow=n, ncol=p))
Xtilde = scale(X[,-1], scale=FALSE)
beta.star=c(1,runif(p,-5,5))
y=X%*%beta.star + sigma*rnorm(n)
fusedtmp = fusedRidge(X, y, lam1, lam2)

#Running simulation to estimate mean and variance of estimates
nrep = 1000
betarep = matrix(0,p,nrep)
ycur = NULL
for (i in 1:nrep) {
  ycur=X%*%beta.star + sigma*rnorm(n)
  betarep[,i] = fusedRidge(X,ycur,lam1, lam2)$b
}

#Getting M
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

betaexpect = apply(betarep,1,mean)
theor_betaexpect = qr.solve(crossprod(Xtilde, Xtilde)+lam1*diag(p)+lam2*M)%*%
            crossprod(Xtilde, Xtilde)%*%beta.star[-1]

betavar = (n-1)/n * cov(t(betarep))
theor_betavar = sigma*qr.solve(crossprod(Xtilde, Xtilde)+lam1*diag(p)+lam2*M)%*%t(Xtilde)%*%
  t(qr.solve(crossprod(Xtilde, Xtilde)+lam1*diag(p)+lam2*M)%*%t(Xtilde))



